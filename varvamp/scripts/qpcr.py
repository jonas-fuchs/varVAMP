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
        # determine if it is in the amplicon and not overlapping with amplicon primers
        if not all((qpcr_probes[probe][1] > left_primer[2], qpcr_probes[probe][2] < right_primer[1])):
            continue
        # check if the probe is at least 5-10°C higher in temperature than the primers
        probe_temp = primers.calc_temp(qpcr_probes[probe][0])
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
