"""
qPCR probe and amplicon design
"""
import re

# LIBS
import seqfold
import itertools
import multiprocessing

# varVAMP
from varvamp.scripts import config
from varvamp.scripts import primers
from varvamp.scripts import reporting


def choose_probe_direction(seq):
    """
    choose the direction of the probe so
    that the probe always has more c than g
    """
    c_count = seq.count("c")
    g_count = seq.count("g")
    direction = "-+"
    if c_count < g_count:
        direction = "-"
    if c_count > g_count:
        direction = "+"

    return direction


def ambiguous_ends(amb_seq):
    """
    determine if kmers have ambiguous ends
    """
    if amb_seq[:1] in config.NUCS and amb_seq[-1] in config.NUCS:
        return False
    else:
        return True


def filter_probe_direction_dependent(seq):
    """
    filter probes for their direction
    """
    return (
            (seq[:1] != "g")  # no 5' g as this leads to quenching
            and (primers.calc_hairpin(seq).tm <= config.PRIMER_HAIRPIN)
            and (config.QPROBE_GC_END[0] <= primers.calc_end_gc(seq) <= config.QPROBE_GC_END[1])
    )


def get_qpcr_probes(kmers, ambiguous_consensus, alignment_cleaned):
    """
    find potential qPCR probes
    """
    probe_candidates = {}
    probe_idx = 0

    for kmer in kmers:
        # filter probe for base params
        if not primers.filter_kmer_direction_independent(kmer[0], config.QPROBE_TMP, config.QPROBE_GC_RANGE,
                                                         config.QPROBE_SIZES):
            continue
        # do not allow ambiguous chars at both ends
        if ambiguous_ends(ambiguous_consensus[kmer[1]:kmer[2]]):
            continue
        # calc penalties analogous to primer search
        base_penalty = primers.calc_base_penalty(kmer[0], config.QPROBE_TMP, config.QPROBE_GC_RANGE,
                                                 config.QPROBE_SIZES)
        per_base_mismatches = primers.calc_per_base_mismatches(
            kmer,
            alignment_cleaned,
            ambiguous_consensus
        )
        permutation_penalty = primers.calc_permutation_penalty(
            ambiguous_consensus[kmer[1]:kmer[2]]
        )
        # determine the direction with more cytosine or set both if 50 %
        direction = choose_probe_direction(kmer[0])
        # create probe dictionary
        if "+" in direction:
            if filter_probe_direction_dependent(kmer[0]):
                probe_name = f"PROBE_{probe_idx}_LEFT"
                three_prime_penalty = primers.calc_3_prime_penalty("+", per_base_mismatches)
                probe_candidates[probe_name] = [kmer[0], kmer[1], kmer[2],
                                                base_penalty + permutation_penalty + three_prime_penalty,
                                                per_base_mismatches, direction]
                probe_idx += 1
        if "-" in direction:
            if filter_probe_direction_dependent(primers.rev_complement(kmer[0])):
                probe_name = f"PROBE_{probe_idx}_RIGHT"
                three_prime_penalty = primers.calc_3_prime_penalty("-", per_base_mismatches)
                probe_candidates[probe_name] = [primers.rev_complement(kmer[0]), kmer[1], kmer[2],
                                                base_penalty + permutation_penalty + three_prime_penalty,
                                                per_base_mismatches, direction]
                probe_idx += 1
    # sort by penalty
    probe_candidates = dict(sorted(probe_candidates.items(), key=lambda x: x[1][3]))

    return probe_candidates


def flanking_primer_subset(primer_list, primer_direction, probe):
    """
    subset for primers flanking the probe and sort by penalty
    """
    subset = []

    if primer_direction == "+":
        window_start = probe[1] - config.QAMPLICON_LENGTH[1] + len(probe[0])
        window_stop = probe[1]
    else:
        window_start = probe[2]
        window_stop = probe[2] + config.QAMPLICON_LENGTH[1] - len(probe[0])
    for primer in primer_list:
        if window_start < primer[1] and primer[2] < window_stop:
            subset.append(primer)

    # sort by penalty and start if same penalty
    subset.sort(key=lambda x: (x[3], x[1]))

    return subset


def hardfilter_amplicon(majority_consensus, left_primer, right_primer):
    """
    hardfilter possible amplicon for length, gc and if there are Ns
    (indicating deletions in amp) or n (no folding analysis possible)
    """
    amplicon_length = right_primer[2] - left_primer[1]
    amplicon_seq = majority_consensus[left_primer[1]:right_primer[2]]
    # check for the right amplicon length
    return (
            (config.QAMPLICON_LENGTH[0] <= amplicon_length <= config.QAMPLICON_LENGTH[1])
            and (config.QAMPLICON_GC[0] <= primers.calc_gc(amplicon_seq) <= config.QAMPLICON_GC[1])
            and "NN" not in amplicon_seq
    )


def forms_dimer_or_overhangs(right_primer, left_primer, probe, ambiguous_consensus):
    """
    checks if primers cause dimers and if combinations of primers/probe including all permutations form dimers
    """

    forms_structure = False

    # first check if there are dimers between the two flanking primers
    if primers.is_dimer(left_primer[0], right_primer[0]):
        return True
    # for the probe check all permutations and possible overhangs to ensure
    # that none of the primers could cause unspecific probe binding.
    # first get all permutations
    probe_per = reporting.get_permutations(ambiguous_consensus[probe[1]:probe[2]])
    left_per = reporting.get_permutations(ambiguous_consensus[left_primer[1]:left_primer[2]])
    right_per = reporting.get_permutations(ambiguous_consensus[right_primer[1]:right_primer[2]])
    # then check all permutations
    for combination in [(probe_per, left_per), (probe_per, right_per)]:
        for oligo1 in combination[0]:
            for oligo2 in combination[1]:
                if primers.is_dimer(oligo1, oligo2):
                    forms_structure = True
                    break
            # break all loops because we found an unwanted structure in one of the permutations
            # (either dimer formation or a too long overlap at the ends of the primer)
            if forms_structure:
                break
        if forms_structure:
            break

    return forms_structure


def assess_amplicons(left_subset, right_subset, qpcr_probes, probe, majority_consensus, ambiguous_consensus):
    """
    assess if a potential amplicon is a qPCR scheme for a specific probe and return the one with the lowest penalty
    """

    primer_combinations = ()
    amplicon_found = False

    # consider a combination of flanking primers if ...
    for left_primer in left_subset:
        for right_primer in right_subset:
            # ... the amplicon is large enough and is in gc range, ...
            if not hardfilter_amplicon(majority_consensus, left_primer, right_primer):
                continue
            # ... the probe is close enough to the primer on the same strand
            if "LEFT" in probe:
                if not qpcr_probes[probe][1] in range(
                        left_primer[2] + config.QPROBE_DISTANCE[0],
                        left_primer[2] + config.QPROBE_DISTANCE[1]
                ):
                    continue
            elif "RIGHT" in probe:
                if not right_primer[1] in range(
                        qpcr_probes[probe][2] + config.QPROBE_DISTANCE[0],
                        qpcr_probes[probe][2] + config.QPROBE_DISTANCE[1]

                ):
                    continue
            # ... the primer temps do not differ too much, ...
            primer_temps = (primers.calc_temp(right_primer[0]), primers.calc_temp(left_primer[0]))
            if abs(primer_temps[0] - primer_temps[1]) > config.QPRIMER_DIFF:
                continue
            # ... the probe has a higher temp than the primers and ...
            probe_temp = primers.calc_temp(qpcr_probes[probe][0])
            if not all(
                    [config.QPROBE_TEMP_DIFF[0] <= probe_temp - x <= config.QPROBE_TEMP_DIFF[1] for x in primer_temps]):
                continue
            # .... all combination of oligos do not form dimers or overhangs.
            if forms_dimer_or_overhangs(right_primer, left_primer, qpcr_probes[probe], ambiguous_consensus):
                continue
            # append to list and break as this is the primer combi
            # with the lowest penalty (primers are sorted by penalty)
            amplicon_found = True
            break
        # break also the outer loop
        if amplicon_found:
            primer_combinations = (left_primer, right_primer)
            break

    return primer_combinations


def find_qcr_schemes(qpcr_probes, left_primer_candidates, right_primer_candidates, majority_consensus,
                     ambiguous_consensus):
    """
    this finds the final qPCR schemes. it slices for primers flanking a probe and
    test all left/right combinations whether they are potential amplicons. as primers
    are sorted by penalty, only the very first match is considered as this has the
    lowest penalty. however, probes are overlapping and there is a high chance that
    left and right primers are found multiple times. to consider only one primer-probe
    combination the probes are also sorted by penalty. therefore, if a primer
    combination has been found already the optimal probe was already selected and
    there is no need to consider this primer probe combination.
    """

    qpcr_scheme_candidates = []
    found_amplicons = []
    amplicon_nr = -1

    for probe in qpcr_probes:
        left_subset = flanking_primer_subset(left_primer_candidates, "+", qpcr_probes[probe])
        right_subset = flanking_primer_subset(right_primer_candidates, "-", qpcr_probes[probe])
        # consider if there are primers flanking the probe ...
        if not left_subset or not right_subset:
            continue
        primer_combination = assess_amplicons(left_subset, right_subset, qpcr_probes, probe, majority_consensus,
                                              ambiguous_consensus)
        # ... a combi has been found, ...
        if not primer_combination:
            continue
        # ...and this combi is not already present for a probe with a better penalty.
        if primer_combination in found_amplicons:
            continue
        # populate the primer dictionary:
        amplicon_nr += 1
        found_amplicons.append(primer_combination)
        qpcr_scheme_candidates.append(
            {
                "id": f"AMPLICON_{amplicon_nr}",
                "penalty": qpcr_probes[probe][3] + primer_combination[0][3] + primer_combination[1][3],
                "PROBE": qpcr_probes[probe],
                "LEFT": primer_combination[0],
                "RIGHT": primer_combination[1]
            }
        )
    # and again sort by total penalty (left + right + probe)
    return qpcr_scheme_candidates


def process_single_amplicon_deltaG(amplicon, majority_consensus):
    """
    Process a single amplicon to test its deltaG and apply filtering.
    This function will be called concurrently by multiple threads.
    """
    seq = majority_consensus[amplicon["LEFT"][1]:amplicon["RIGHT"][2]]
    seq = seq.replace("N", "")
    seq = seq.replace("n", "")
    # check if the amplicon overlaps with an amplicon that was previously
    # found and had a high enough deltaG
    min_temp = min((primers.calc_temp(amplicon["LEFT"][0]),
                    primers.calc_temp(amplicon["RIGHT"][0])))
    # calculate deltaG at the minimal primer temp
    amplicon["deltaG"] = seqfold.dg(seq, min_temp)

    return amplicon


def test_amplicon_deltaG_parallel(qpcr_schemes_candidates, majority_consensus, n_to_test, deltaG_cutoff, n_threads):
    """
    Test all amplicon deltaGs for the top n hits at the lowest primer temperature
    and filters if they fall below the cutoff. Multiple processes are used
    for processing amplicons in parallel.
    """
    final_amplicons = []

    # Create a pool of processes to handle the concurrent processing
    with multiprocessing.Pool(processes=n_threads) as pool:
        # Create a list of the first n amplicon tuples for processing
        # The list is sorted first on whether offset targets were predicted for the amplicon,
        # then by penalty. This ensures that amplicons with offset targets are always considered last
        amplicons = itertools.islice(
            sorted(qpcr_schemes_candidates, key=lambda x: (x.get("offset_targets", False), x["penalty"])),
            n_to_test
        )
        # process amplicons concurrently
        results = pool.starmap(process_single_amplicon_deltaG, [(amp, majority_consensus) for amp in amplicons])
        # Process the results
        retained_ranges = []
        for amp in results:
            # check if the amplicon overlaps with an amplicon that was previously
            # found and had a high enough deltaG
            if amp["deltaG"] <= deltaG_cutoff:
                continue
            amp_range = range(amp["LEFT"][1], amp["RIGHT"][2])
            overlaps_retained = False
            for r in retained_ranges:
                if amp_range.start < r.stop and r.start < amp_range.stop:
                    overlaps_retained = True
                    break
            if not overlaps_retained:
                final_amplicons.append(amp)
                retained_ranges.append(amp_range)

    return final_amplicons
