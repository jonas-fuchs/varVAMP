#!/usr/bin/env python
"""
            INFO
------------------------------
This contains the main workflow.

           LICENCE
-------------------------------
todo

         EXPLANATIONS
-------------------------------
This contains the definitions for evaluating if a kmer is a potential
primer (the core is highly similar to primalscheme!). Importantly, the user
can specify if the primer can contain ambiguous characters at the 3' end.
If passing kmers have the same start only the kmer with the min base
penalty is retained. Further the ambigous version of the primer is checked
for mismatches against each sequence in the alignment. This allows to calculate
the penalty score for mismatches at the 3' end. Also the number of permutations
is calculated and multiplied by the permutation penalty.
"""

# BUILT-INS
import itertools

# LIBS
from Bio.Seq import Seq
import primer3 as p3

# varVAMP
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
    return p3.calc_tm(
            primer.upper(),
            mv_conc = config.MV_CONC,
            dv_conc = config.DV_CONC,
            dntp_conc = config.DNTP_CONC,
            dna_conc = config.DNA_CONC
        )

def calc_hairpin(primer):
    """
    calculates hairpins
    """
    return p3.calc_hairpin(
            primer.upper(),
            mv_conc = config.MV_CONC,
            dv_conc = config.DV_CONC,
            dntp_conc = config.DNTP_CONC,
            dna_conc = config.DNA_CONC
        )


def calc_dimer(primer1, primer2):
    """
    Calculate the heterodimerization thermodynamics of two DNA sequences.
    Return primer3 thermo object.
    """
    return p3.calc_heterodimer(
        primer1.upper(),
        primer2.upper(),
        mv_conc=config.MV_CONC,
        dv_conc=config.DV_CONC,
        dna_conc=config.DNA_CONC,
        dntp_conc=config.DNTP_CONC,
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
    poly x, dinucleotide repeats and homodimerization.
    """
    return(
        (config.PRIMER_TMP[0] <= calc_temp(primer) <= config.PRIMER_TMP[1]) and
        (config.PRIMER_GC_RANGE[0] <= calc_gc(primer) <= config.PRIMER_GC_RANGE[1]) and
        (calc_max_polyx(primer) <= config.MAX_POLYX) and
        (calc_max_dinuc_repeats(primer) <= config.MAX_DINUC_REPEATS) and
        (calc_dimer(primer, primer).tm <= config.MAX_DIMER_TMP)
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
        start = 2   # stop index as it is the start of the reverse primer
    hardfiltered_kmers.sort(key = lambda x:(x[start], x[3]))

    # pick the lowest scoring primer for primers that have the same start
    prev_start = -1
    for primer in hardfiltered_kmers:
        if primer[start] != prev_start:
            candidates.append(primer)
            prev_start = primer[start]

    return candidates


def get_permutation_penalty(ambiguous_primer):
    # get all permutations of a primer with ambiguous nucleotides and
    # multiply with permutation penalty
    permutations = 0

    for nuc in ambiguous_primer:
        if nuc in config.ambig_nucs:
            n = len(config.ambig_nucs[nuc])
            if permutations != 0:
                permutations = permutations*n
            else:
                permutations = n

    return permutations*config.PRIMER_PERMUTATION_PENALTY


def primer_per_base_mismatch(primer, alignment, ambiguous_consensus):
    """
    calculate for a given primer with [seq, start, stop]
    percent mismatch per primer pos with the alignment.
    considers if primer or sequences have an amb nuc.
    """
    primer_per_base_mismatch = len(primer[0])*[0]
    ambigous_primer = ambiguous_consensus[primer[1]:primer[2]]


    for sequence in alignment:
        # slice each sequence of the alignment for the primer
        # start and stop positions
        seq_slice = sequence[1][primer[1]:primer[2]]
        for idx, slice_nuc in enumerate(seq_slice):
            # find the respective nuc to that of the slice
            current_primer_pos = ambigous_primer[idx]
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

def get_penalty_3_prime(direction, primer):
    """
    calculate the penalty for mismatches at the 3' end.
    the more mismatches are closer to the 3' end of the primer,
    the higher the penalty.
    """
    penalty = 0
    if config.PRIMER_3_PENALTY:
        for i in range(0,len(config.PRIMER_3_PENALTY)):
            if direction == "RIGHT":
                penalty += primer[4][i]*config.PRIMER_3_PENALTY[i]
            elif direction == "LEFT":
                penalty += primer[4][len(primer[0])-1-i]*config.PRIMER_3_PENALTY[i]

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
    # based on this score calculate the 3' penalty and n ambiguous chars.
    for direction, primer_candidates in [("LEFT", left_primer_candidates), ("RIGHT", right_primer_candidates)]:
        for primer in primer_candidates:
            primer.append(primer_per_base_mismatch(
                primer,
                alignment,
                ambiguous_consensus
            ))
            primer[3] += get_penalty_3_prime(direction, primer) #add 3 prime penalty to base penalty
            primer[3] += get_permutation_penalty( #also add permutation penalty
                ambiguous_consensus[primer[1]:primer[2]]
            )

    return left_primer_candidates, right_primer_candidates


def create_primer_dictionary(primer_candidates, direction):
    """
    creates a primer dictionary from primer list
    """
    primer_dict = {}
    primer_idx = 0

    for primer in primer_candidates:
        primer_name = direction + "_" + str(primer_idx)
        primer_dict[primer_name] = primer
        primer_idx += 1

    return primer_dict

def find_best_primers(left_primer_candidates, right_primer_candidates):
    """
    Primer candidates might be overlapping. Here all primers are found within a
    window that is defined by the start of the first primer of the window and its
    stop + maximum primer length. From this list the best scoring primer is choosen.
    This reduces the complexity of the later calculated amplicon graph while
    retaining well scoring primers within a reasonable range.
    """
    for direction, primer_candidates in [("LEFT", left_primer_candidates), ("RIGHT", right_primer_candidates)]:
        # ini lists
        primers_to_retain = []
        primers_temp = []
        writeable = False
        # set the first search window in relation to the first primer
        best_scoring = [0,primer_candidates[0][3]]
        search_window = [primer_candidates[0][1], primer_candidates[0][2]+config.PRIMER_SIZES[1]]

        for idx, primer in enumerate(primer_candidates):
            # append primers if their start is within the current window
            if primer[1] in range(search_window[0], search_window[1]):
                primers_temp.append([idx, primer[3]])
                writeable = False
                # if the end of the primer list is reached, make the appended
                # primers writeable
                if idx == len(primer_candidates)-1:
                    writeable = True
            else:
                writeable = True

            if writeable:
                # check in the temporary list of primers for the one with
                # the best score...
                for primer_temp in primers_temp:
                    if primer_temp[1] < best_scoring[1]:
                        best_scoring = primer_temp
                # ... and append to final list
                primers_to_retain.append(best_scoring[0])
                # if the end has not been reached, initialize a new window
                if idx != len(primer_candidates)-1:
                    search_window = [primer[1], primer[2]+config.PRIMER_SIZES[1]]
                    best_scoring = [idx, primer[3]]
                    primers_temp = [[idx, primer[3]]]

        # subset the primer candidate list with the remembered indices
        subset = [primer_candidates[i] for i in primers_to_retain]
        # and put these in the right lists
        if direction == "LEFT":
            left_primer_candidates = create_primer_dictionary(subset, direction)
        if direction == "RIGHT":
            right_primer_candidates = create_primer_dictionary(subset, direction)

    return left_primer_candidates, right_primer_candidates
