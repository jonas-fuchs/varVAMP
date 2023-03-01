from Bio.Seq import Seq
import primer3 as p3

from scr import config


def calc_gc(primer):
    """
    calculate the gc of a sequence
    """
    return 100*(primer.count("g")+primer.count("G")+primer.count("c")+primer.count("C"))/len(primer)

def calc_temp(primer):
    """
    calculate the melting temperature
    """
    return p3.calcTm(
            primer,
            mv_conc = config.MV_CONC,
            dv_conc = config.DV_CONC,
            dntp_conc = config.DNTP_CONC,
            dna_conc = config.DNA_CONC
        )

def calc_hairpin(primer):
    """
    calculates hairpins
    """
    return p3.calcHairpin(
            primer,
            mv_conc = config.MV_CONC,
            dv_conc = config.DV_CONC,
            dntp_conc = config.DNTP_CONC,
            dna_conc = config.DNA_CONC
        )

def calc_max_polyx(primer):
    """
    calculate maximum polyx in primers
    """
    previous_nuc = primer[0]
    counter = 1
    max_polyx = 1
    for nuc in primer[1:]:
        if nuc == previous_nuc:
            counter += 1
        else:
            counter = 1
            previous_nuc = nuc
        if counter > max_polyx:
            max_polyx = counter
    return(max_polyx)

def calc_max_dinuc_repeats(primer):
    """
    calculate the amount of repeating dinucleotides
    """
    for s in [primer, primer[1:]]:
        previous_dinuc = s[0:2]
        max_dinuc = 0
        counter = 0
        for i in range(2,len(s),2):
            if s[i:i+2] == previous_dinuc:
                counter += 1
            else:
                if counter > max_dinuc:
                    max_dinuc = counter
                counter = 0
                previous_dinuc = s[i:i+2]
    return max_dinuc

def three_prime_ambiguous(amb_primer):
    """
    determine if a sequence contains an ambiguous char at the 3'prime
    """
    len_3_prime = config.MIN_3_WITHOUT_AMB

    if len_3_prime != 0:
        for nuc in amb_primer[len(amb_primer)-len_3_prime:]:
            if nuc not in config.nucs:
                is_amb = True
                break
            else:
                is_amb = False
    else:
        is_amb = False

    return is_amb

def calc_base_penalty(primer):
    """
    Calculate intrinsic primer penalty.
    """
    penalty = 0

    tm = calc_temp(primer)
    gc = calc_gc(primer)
    size = len(primer)

    # TEMP penalty
    if tm > config.PRIMER_TMP[2]:
        penalty += config.PRIMER_TM_PENALTY*(
            tm - config.PRIMER_TMP[2]
            )
    if tm < config.PRIMER_TMP[2]:
        penalty += config.PRIMER_TM_PENALTY*(
            config.PRIMER_TMP[2] - tm
            )
    # GC penalty
    if gc > config.PRIMER_GC_RANGE[2]:
        penalty += config.PRIMER_GC_PENALTY*(
            gc - config.PRIMER_GC_RANGE[2]
            )
    if gc < config.PRIMER_GC_RANGE[2]:
        penalty += config.PRIMER_GC_PENALTY*(
            config.PRIMER_GC_RANGE[2] - gc
        )
    # SIZE penalty
    if size > config.PRIMER_SIZES[2]:
        penalty += config.PRIMER_SIZE_PENALTY*(
            size - config.PRIMER_SIZES[2]
        )
    if size < config.PRIMER_SIZES[2]:
        penalty += config.PRIMER_SIZE_PENALTY * (
            config.PRIMER_SIZES[2] - size
        )

    return penalty

def rev_complement(seq):
    """
    reverse complement a sequence
    """
    return(str(Seq(seq).reverse_complement()))

def hardfilter_primers(primer):
    """
    hard filter primers for temperature, gc content,
    poly x and dinucleotide repeats.
    """
    return(
        (config.PRIMER_TMP[0] <= calc_temp(primer) <= config.PRIMER_TMP[1]) and
        (config.PRIMER_GC_RANGE[0] <= calc_gc(primer) <= config.PRIMER_GC_RANGE[1]) and
        (calc_max_polyx(primer) <= config.MAX_POLYX) and
        (calc_max_dinuc_repeats(primer) <= config.MAX_DINUC_REPEATS)
    )

def filter_primer_direction_specific(direction, primer, ambiguous_consensus):
    """
    filter for 3'ambiguous and hairpin - this differs
    depending on the direction of the primer.
    """
    if direction == "LEFT":
        amb_kmer = ambiguous_consensus[primer[1]:primer[2]]
        hairpin_tm = calc_hairpin(primer[0]).tm
    elif direction == "RIGHT":
        amb_kmer = rev_complement(ambiguous_consensus[primer[1]:primer[2]])
        hairpin_tm = calc_hairpin(rev_complement(primer[0])).tm
    if hairpin_tm <= config.PRIMER_HAIRPIN:
        if not three_prime_ambiguous(amb_kmer):
            return primer

def find_lowest_scoring(direction, hardfiltered_kmers):
    """
    sort the primers by start and base penalty.
    for primers that have the same start, retain only
    the lowest scoring and now give a fixed number.
    """
    candidates = []

    # sort for start of the primer and score
    if direction == "LEFT":
        start = 1   # start index of the forward primer
    elif direction == "RIGHT":
        start = 2   # stop index asit is the start of the reverse primer
    hardfiltered_kmers.sort(key = lambda x:(x[start], x[3]))

    # pick the lowest scoring primer for primers that have the same start
    primer_index = 1
    prev_start = -1
    for primer in hardfiltered_kmers:
        if primer[start] != prev_start:
            primer.append(direction+"_"+str(primer_index))
            candidates.append(primer)
            prev_start = primer[start]
            primer_index += 1

    return candidates

def primer_per_base_mismatch(primer, alignment):
    """
    calculate for a given primer with [seq, start, stop]
    percent mismatch per primer pos with the alignment.
    considers if primer or sequences have an amb nuc.
    """
    primer_per_base_mismatch = len(primer[0])*[0]

    for sequence in alignment:
        # slice each sequence of the alignment for the primer
        # start and stop positions
        seq_slice = sequence[1][primer[1]:primer[2]]
        for idx, slice_nuc in enumerate(seq_slice):
            # find the respective nuc to that of the slice
            current_primer_pos = primer[0][idx]
            if slice_nuc != current_primer_pos:
                # check if the slice nucleotide is an amb pos
                if slice_nuc in config.ambig_nucs:
                    # check if the primer has an amb pos
                    if current_primer_pos in config.ambig_nucs:
                        slice_nuc_set = set(config.ambig_nucs[slice_nuc])
                        pri_set = set(config.ambig_nucs[current_primer_pos])
                        # check if these sets have no overlap
                        # -> mismatch
                        if len(slice_nuc_set.intersection(pri_set)) == 0:
                            primer_per_base_mismatch[idx] += 1
                    # if no amb pos is in primer then check if primer nuc
                    # is part of the amb slice nuc
                    elif current_primer_pos not in config.ambig_nucs[slice_nuc]:
                            primer_per_base_mismatch[idx] += 1
                # check if primer has an amb pos but the current
                # slice_nuc is not part of this amb nucleotide
                elif current_primer_pos in config.ambig_nucs:
                        if slice_nuc not in config.ambig_nucs[current_primer_pos]:
                            primer_per_base_mismatch[idx] += 1
                # mismatch
                else:
                    primer_per_base_mismatch[idx] += 1

    # gives a percent mismatch over all positions of the primer from 5' to 3'
    primer_per_base_mismatch = [round(x/len(alignment),2) for x in primer_per_base_mismatch]

    return primer_per_base_mismatch

def penalty_3_prime(direction, primer):
    """
    calculate the penalty for mismatches at the 3' end.
    the more mismatches are closer to the 3' end of the primer,
    the higher the penalty.
    """
    penalty = 0

    for i in range(0,len(config.PRIMER_3_PENALTY)):
        if direction == "RIGHT":
            penalty += primer[5][i]*config.PRIMER_3_PENALTY[i]
        elif direction == "LEFT":
            penalty += primer[5][len(primer[0])-1-i]*config.PRIMER_3_PENALTY[i]

    return(penalty)

def find_primers(kmers, ambiguous_consensus, alignment):
    """
    hardfilter kmers and further filter
    for potential primers
    """

    hardfiltered_left_kmers = []
    hardfiltered_right_kmers = []

    for kmer in kmers:
        if hardfilter_primers(kmer[0]):
            base_penalty = calc_base_penalty(kmer[0])
            if base_penalty <= config.PRIMER_MAX_BASE_PENALTY:
                if filter_primer_direction_specific(
                    "LEFT",
                    kmer,
                    ambiguous_consensus
                ):
                    hardfiltered_left_kmers.append(
                        [kmer[0], kmer[1], kmer[2], base_penalty]
                    )
                if filter_primer_direction_specific(
                    "RIGHT",
                    kmer,
                    ambiguous_consensus
                ):
                    hardfiltered_right_kmers.append(
                        [kmer[0], kmer[1], kmer[2], base_penalty]
                    )

    # filter kmers and complement kmers for possible primers
    left_primer_candidates = find_lowest_scoring("LEFT", hardfiltered_left_kmers)
    right_primer_candidates = find_lowest_scoring("RIGHT", hardfiltered_right_kmers)

    # now calculate the mismatches for each position in the primer.
    # based on this score calculate the 3' penalty and add to base penalty.
    for direction, primer_candidates in [("LEFT", left_primer_candidates), ("RIGHT", right_primer_candidates)]:
        for primer in primer_candidates:
            primer.append(primer_per_base_mismatch(primer, alignment))
            primer[3] = primer[3] + penalty_3_prime(direction, primer)

    return left_primer_candidates, right_primer_candidates
