# LIBS
import seqfold

# varVAMP
from varvamp.scripts import config
from varvamp.scripts import primers

# overwrite config variables to use primer functions
config.PRIMER_TMP = config.QPROBE_TMP
config.PRIMER_GC_RANGE = config.QPROBE_GC_RANGE
config.PRIMER_SIZES = config.QPROBE_SIZES
config.PRIMER_GC_CLAMP = config.QPROBE_GC_CLAMP
config.PRIMER_MAX_GC_END = config.QPROBE_MAX_GC_END


def choose_probe_direction(seq):
    """
    choose the direction of the probe so
    that the probe always has more c than g
    """
    c_count = seq.count("c")
    g_count = seq.count("g")
    if c_count < g_count:
        return "-"
    elif c_count == g_count:
        return "-+"
    elif c_count > g_count:
        return "+"


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
        and (config.PRIMER_GC_CLAMP <= primers.calc_end_gc(seq) <= config.PRIMER_MAX_GC_END)
    )


def get_qpcr_probes(kmers, ambiguous_consensus, alignment_cleaned):
    """
    find potential qPCR probes
    """
    probe_candidates = {}
    probe_idx = 0

    for kmer in kmers:
        # filter probe for base params
        if not primers.filter_kmer_direction_independent(kmer[0]):
            continue
        # do not allow ambiguous chars at both ends
        if ambiguous_ends(ambiguous_consensus[kmer[1]:kmer[2]]):
            continue
        # calc penalties analogous to primer search
        base_penalty = primers.calc_base_penalty(kmer[0])
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
                probe_candidates[probe_name] = [kmer[0], kmer[1], kmer[2], base_penalty + permutation_penalty + three_prime_penalty, per_base_mismatches]
                probe_idx += 1
        if "-" in direction:
            if filter_probe_direction_dependent(primers.rev_complement(kmer[0])):
                probe_name = f"PROBE_{probe_idx}_RW"
                three_prime_penalty = primers.calc_3_prime_penalty("-", per_base_mismatches)
                probe_candidates[probe_name] = [primers.rev_complement(kmer[0]), kmer[1], kmer[2], base_penalty + permutation_penalty + three_prime_penalty, per_base_mismatches]
                probe_idx += 1

    return probe_candidates


def find_best_compatible_probe(qpcr_probes, left_primer, right_primer, primer_temps):
    """
    find the lowest internal probe for a given amplicon
    """

    best_probe = ("none", float("inf"))

    for probe in qpcr_probes:
        probe_temp = primers.calc_temp(qpcr_probes[probe][0])
        # check if probe is close enough to the primer on the same strand
        # and does not overlap with the primer on the other strand
        if "FW" in probe:
            if not all(
                (
                    qpcr_probes[probe][1] in range(
                        left_primer[2] + config.QPROBE_DISTANCE[0],
                        left_primer[2] + config.QPROBE_DISTANCE[1] + 1
                    ),
                    qpcr_probes[probe][2] < right_primer[1]
                )
            ):
                continue
        elif "RW" in probe:
            if not all(
                (
                    right_primer[1] in range(
                        qpcr_probes[probe][2] + config.QPROBE_DISTANCE[0],
                        qpcr_probes[probe][2] + config.QPROBE_DISTANCE[1] + 1
                    ),
                    left_primer[2] < qpcr_probes[probe][1]
                )
            ):
                continue
        # check if the probe is at least 5-10°C higher in temperature than the primers
        if all([5 <= probe_temp-x <= 10 for x in primer_temps]):
            continue
        # check if the probe forms dimers with the primers
        if any(
            (
                primers.calc_dimer(right_primer[0], qpcr_probes[probe][0]).tm > config.PRIMER_MAX_DIMER_TMP,
                primers.calc_dimer(left_primer[0], qpcr_probes[probe][0]).tm > config.PRIMER_MAX_DIMER_TMP
            )
        ):
            continue
        # determine for each amplicon the lowest scoring probe
        if qpcr_probes[probe][1] < best_probe[1]:
            best_probe = (probe, qpcr_probes[probe][1])

    return best_probe


def populate_qPCR_dictionary(qpcr_amplicons, qpcr_probes, best_probe, left_primer, right_primer):
    """
    populate qPCR dictionary and update if a amplicon with a better score is found
    """
    if not best_probe[0] in qpcr_amplicons:
        qpcr_amplicons[best_probe[0]] = {
            "score": qpcr_probes[best_probe[0]][3] + left_primer[3] + right_primer[3],
            "probe": qpcr_probes[best_probe[0]],
            "left": left_primer,
            "right": right_primer
        }
    else:
        if qpcr_amplicons[best_probe[0]]["score"] > qpcr_probes[best_probe[0]][3] + left_primer[3] + right_primer[3]:
            qpcr_amplicons[best_probe[0]] = {
                "score": qpcr_probes[best_probe[0]][3] + left_primer[3] + right_primer[3],
                "probe": qpcr_probes[best_probe[0]],
                "left": left_primer,
                "right": right_primer
            }


def sort_qPCR_dictionary(qpcr_amplicons):
    """
    sort dictionary by cummulative score
    """
    sorted_qpcr_amplicons = {}
    for key in sorted(qpcr_amplicons, key=lambda x: qpcr_amplicons[x]["score"]):
        sorted_qpcr_amplicons[key] = qpcr_amplicons[key]

    return sorted_qpcr_amplicons


def find_qPCR_schemes(qpcr_probes, left_primer_candidates, right_primer_candidates):
    """
    for each qPCR probe find the best possible amplicon primers
    """

    qpcr_amplicons = {}

    for left_primer in left_primer_candidates:
        for right_primer in right_primer_candidates:
            amplicon_length = right_primer[2] - left_primer[1]
            # check for the right amplicon length
            if not config.QAMPLICON_LENGTH[0] <= amplicon_length <= config.QAMPLICON_LENGTH[1]:
                continue
            primer_temps = (primers.calc_temp(right_primer[0]), primers.calc_temp(left_primer[0]))
            # check if temperature diff is lower than 2°C
            if abs(primer_temps[0] - primer_temps[1]) > 2:
                continue
            # check if primers form dimer
            if primers.calc_dimer(left_primer[0], right_primer[0]).tm > config.PRIMER_MAX_DIMER_TMP:
                continue
            # check for internal oligos and choose the one with the lowest score
            best_probe = find_best_compatible_probe(qpcr_probes, left_primer, right_primer, primer_temps)
            if best_probe[1] == float("inf"):
                continue
            # create dic for probe and the best amplicon
            populate_qPCR_dictionary(qpcr_amplicons, qpcr_probes, best_probe, left_primer, right_primer)

    return sort_qPCR_dictionary(qpcr_amplicons)
