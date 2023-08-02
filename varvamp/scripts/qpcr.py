"""
qPCR probe and amplicon design
"""

# LIBS
import seqfold
import itertools
import multiprocessing

# varVAMP
from varvamp.scripts import config
from varvamp.scripts import primers

def choose_probe_direction(seq):
    """
    choose the direction of the probe so
    that the probe always has more c than g
    """
    c_count = seq.count("c")
    g_count = seq.count("g")
    if c_count < g_count:
        direction = "-"
    if c_count == g_count:
        direction = "-+"
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
    return(
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
        if not primers.filter_kmer_direction_independent(kmer[0], config.QPROBE_TMP, config.QPROBE_GC_RANGE, config.QPROBE_SIZES):
            continue
        # do not allow ambiguous chars at both ends
        if ambiguous_ends(ambiguous_consensus[kmer[1]:kmer[2]]):
            continue
        # calc penalties analogous to primer search
        base_penalty = primers.calc_base_penalty(kmer[0], config.QPROBE_TMP, config.QPROBE_GC_RANGE, config.QPROBE_SIZES)
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
                probe_name = f"PROBE_{probe_idx}_FW"
                three_prime_penalty = primers.calc_3_prime_penalty("+", per_base_mismatches)
                probe_candidates[probe_name] = [kmer[0], kmer[1], kmer[2], base_penalty + permutation_penalty + three_prime_penalty, per_base_mismatches, direction]
                probe_idx += 1
        if "-" in direction:
            if filter_probe_direction_dependent(primers.rev_complement(kmer[0])):
                probe_name = f"PROBE_{probe_idx}_RW"
                three_prime_penalty = primers.calc_3_prime_penalty("-", per_base_mismatches)
                probe_candidates[probe_name] = [primers.rev_complement(kmer[0]), kmer[1], kmer[2], base_penalty + permutation_penalty + three_prime_penalty, per_base_mismatches, direction]
                probe_idx += 1
    # sort by score
    probe_candidates = dict(sorted(probe_candidates.items(), key=lambda x: x[1][3]))

    return probe_candidates


def flanking_primer_subset(primer_list, primer_direction, probe):
    """
    subset for primers flanking the probe and sort by score
    """
    subset = []

    if primer_direction == "+":
        window_start = probe[1] - config.QAMPLICON_LENGTH[1] + len(probe[0])
        window_stop = probe[1]
    elif primer_direction == "-":
        window_start = probe[2]
        window_stop = probe[2] + config.QAMPLICON_LENGTH[1] - len(probe[0])
    for primer in primer_list:
        if window_start < primer[1] and primer[2] < window_stop:
            subset.append(primer)

    # sort by score
    subset.sort(key=lambda x: x[3])

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
        and not "NN" in amplicon_seq
    )


def forms_dimer(right_primer, left_primer, probe):
    """
    checks if combinations of primers/probe form dimers
    """

    forms_dimer = False

    for combination in [(right_primer[0], left_primer[0]), (right_primer[0], probe[0]), (left_primer[0], probe[0])]:
        if primers.calc_dimer(combination[0], combination[1]).tm > config.PRIMER_MAX_DIMER_TMP:
            forms_dimer = True
            break

    return forms_dimer


def assess_amplicons(left_subset, right_subset, qpcr_probes, probe, majority_consensus):
    """
    assess if a potential amplicon is a qPCR scheme for a specific probe and return the best scoring
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
            if "FW" in probe:
                if not qpcr_probes[probe][1] in range(
                            left_primer[2] + config.QPROBE_DISTANCE[0],
                            left_primer[2] + config.QPROBE_DISTANCE[1] + 1
                ):
                    continue
            elif "RW" in probe:
                if not right_primer[1] in range(
                            qpcr_probes[probe][2] + config.QPROBE_DISTANCE[0],
                            qpcr_probes[probe][2] + config.QPROBE_DISTANCE[1] + 1

                ):
                    continue
            # ... the primer temps do not differ too much, ...
            primer_temps = (primers.calc_temp(right_primer[0]), primers.calc_temp(left_primer[0]))
            if abs(primer_temps[0] - primer_temps[1]) > config.QPRIMER_DIFF:
                continue
            # ... the probe has a higher temp than the primers and ...
            probe_temp = primers.calc_temp(qpcr_probes[probe][0])
            if not all([config.QPROBE_TEMP_DIFF[0] <= probe_temp-x <= config.QPROBE_TEMP_DIFF[1] for x in primer_temps]):
                continue
            # .... all combination of oligos do not form dimers.
            if forms_dimer(right_primer, left_primer, qpcr_probes[probe]):
                continue
            # append to list and break as this is the lowest scoring primer combi (primers are sorted by score)
            amplicon_found = True
            break
        # break also the outer loop
        if amplicon_found:
            primer_combinations = (left_primer, right_primer)
            break

    return primer_combinations


def find_qcr_schemes(qpcr_probes, left_primer_candidates, right_primer_candidates, majority_consensus):
    """
    this finds the final qPCR schemes. it slices for primers flanking a probe and
    test all left/right combinations whether they are potential amplicons. as primers
    are sorted by score, only the very first match is considered as this has the
    lowest score. however, probes are overlapping and there is a high chance that
    left and right primers are found multiple times. to consider only one primer-
    probe combination the probes are also sorted by score. therefore, if a primer
    combination has been found already the optimal probe was already selected and
    there is no need to consider this primer probe combination.
    """

    qpcr_scheme_candidates = {}
    found_amplicons = []
    amplicon_nr = -1

    for probe in qpcr_probes:
        left_subset = flanking_primer_subset(left_primer_candidates, "+", qpcr_probes[probe])
        right_subset = flanking_primer_subset(right_primer_candidates, "-", qpcr_probes[probe])
        # consider if there are primers flanking the probe ...
        if not left_subset or not right_subset:
            continue
        primer_combination = assess_amplicons(left_subset, right_subset, qpcr_probes, probe, majority_consensus)
        # ... a combi has been found, ...
        if not primer_combination:
            continue
        # ...and this combi is not already present for a probe with a better score.
        if primer_combination in found_amplicons:
            continue
        # populate the primer dictionary:
        amplicon_nr += 1
        found_amplicons.append(primer_combination)
        qpcr_scheme_candidates[f"AMPLICON_{amplicon_nr}"] = {
            "score": qpcr_probes[probe][3]+primer_combination[0][3]+primer_combination[1][3],
            "probe": qpcr_probes[probe],
            "left": primer_combination[0],
            "right": primer_combination[1]
        }
    # and again sort by total score (left + right + probe)
    qpcr_scheme_candidates = dict(sorted(qpcr_scheme_candidates.items(), key=lambda x: x[1]["score"]))

    return qpcr_scheme_candidates


def process_single_amplicon_deltaG(amplicon, majority_consensus):
    """
    Process a single amplicon to test its deltaG and apply filtering.
    This function will be called concurrently by multiple threads.
    """
    name, data = amplicon
    start = data["left"][1]
    stop = data["right"][2]
    seq = majority_consensus[start:stop]
    seq = seq.replace("N", "")
    seq = seq.replace("n", "")
    amp_positions = list(range(start, stop + 1))
    # check if the amplicon overlaps with an amplicon that was previously
    # found and had a high enough deltaG
    min_temp = min((primers.calc_temp(data["left"][0]),
                    primers.calc_temp(data["right"][0])))
    # calculate deltaG at the minimal primer temp
    deltaG = seqfold.dg(seq, min_temp)

    return deltaG, amp_positions, name


def test_amplicon_deltaG_parallel(qpcr_schemes_candidates, majority_consensus, n_to_test, deltaG_cutoff):
    """
    Test all amplicon deltaGs for the top n hits at the lowest primer temperature
    and filters if they fall below the cutoff. Multiple processes are used
    for processing amplicons in parallel.
    """
    final_schemes = {}
    passed_counter = 0  # counter for re-naming amplicons that passed deltaG cutoff
    amplicon_set = set()

    # Create a pool of processes to handle the concurrent processing
    with multiprocessing.Pool() as pool:
        # Create a list of the first n amplicon tuples for processing
        amplicons = itertools.islice(qpcr_schemes_candidates.items(), n_to_test)
        # process amplicons concurrently
        results = pool.starmap(process_single_amplicon_deltaG, [(amp, majority_consensus) for amp in amplicons])

        # Process the results
        for deltaG, amp_positions, amp_name in results:
            # check if the amplicon overlaps with an amplicon that was previously
            # found and had a high enough deltaG
            if any(x in amp_positions for x in amplicon_set):
                continue
            # and if this passes cutoff make a dict entry and do not allow further
            # amplicons in that region (they will have a lower score)
            if deltaG > deltaG_cutoff:
                new_name = f"QPCR_SCHEME_{passed_counter}"
                final_schemes[new_name] = qpcr_schemes_candidates[amp_name]
                final_schemes[new_name]["deltaG"] = deltaG
                amplicon_set.update(amp_positions)
                passed_counter += 1

    return final_schemes

#TODO: add option for number of threads (currently all are used)