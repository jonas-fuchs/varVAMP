"""
primer creation and evaluation
"""

# LIBS
from Bio.Seq import Seq
import primer3 as p3

# varVAMP
from varvamp.scripts import config


def calc_gc(seq):
    """
    calculate the gc of a sequence
    """
    return 100*(seq.count("g")+seq.count("c"))/len(seq)


def calc_temp(seq):
    """
    calculate the melting temperature
    """
    return p3.calc_tm(
            seq.upper(),
            mv_conc=config.PCR_MV_CONC,
            dv_conc=config.PCR_DV_CONC,
            dntp_conc=config.PCR_DNTP_CONC,
            dna_conc=config.PCR_DNA_CONC
        )


def calc_hairpin(seq):
    """
    calculates hairpins
    """
    return p3.calc_hairpin(
            seq.upper(),
            mv_conc=config.PCR_MV_CONC,
            dv_conc=config.PCR_DV_CONC,
            dntp_conc=config.PCR_DNTP_CONC,
            dna_conc=config.PCR_DNA_CONC
        )


def calc_dimer(seq1, seq2):
    """
    Calculate the heterodimerization thermodynamics of two DNA sequences.
    Return primer3 thermo object.
    """
    return p3.calc_heterodimer(
        seq1.upper(),
        seq2.upper(),
        mv_conc=config.PCR_MV_CONC,
        dv_conc=config.PCR_DV_CONC,
        dna_conc=config.PCR_DNA_CONC,
        dntp_conc=config.PCR_DNTP_CONC,
    )


def calc_max_polyx(seq):
    """
    calculate maximum polyx of a seq
    """
    previous_nuc = seq[0]
    counter = 0
    max_polyx = 0
    for nuc in seq[1:]:
        if nuc == previous_nuc:
            counter += 1
        else:
            counter = 0
            previous_nuc = nuc
        if counter > max_polyx:
            max_polyx = counter
    return(max_polyx)


def calc_max_dinuc_repeats(seq):
    """
    calculate the amount of repeating
    dinucleotides in a sequence
    """
    for s in [seq, seq[1:]]:
        previous_dinuc = s[0:2]
        max_dinuc = 0
        counter = 0
        for i in range(2, len(s), 2):
            if s[i:i+2] == previous_dinuc:
                counter += 1
            else:
                if counter > max_dinuc:
                    max_dinuc = counter
                counter = 0
                previous_dinuc = s[i:i+2]
    return max_dinuc


def calc_base_penalty(seq):
    """
    Calculate intrinsic primer penalty.
    """
    penalty = 0

    tm = calc_temp(seq)
    gc = calc_gc(seq)
    size = len(seq)

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


def filter_kmer_direction_independent(kmer):
    """
    filter kmer for temperature, gc content,
    poly x, dinucleotide repeats and homodimerization
    """
    return(
        (config.PRIMER_TMP[0] <= calc_temp(kmer) <= config.PRIMER_TMP[1])
        and (config.PRIMER_GC_RANGE[0] <= calc_gc(kmer) <= config.PRIMER_GC_RANGE[1])
        and (calc_max_polyx(kmer) <= config.PRIMER_MAX_POLYX)
        and (calc_max_dinuc_repeats(kmer) <= config.PRIMER_MAX_DINUC_REPEATS)
        and (calc_dimer(kmer, kmer).tm <= config.PRIMER_MAX_DIMER_TMP)
        and (calc_dimer(kmer, kmer[-5:]).tm <= config.PRIMER_MAX_DIMER_TMP_3_PRIME)
    )


def calc_end_gc(seq):
    """
    check how many gc nucleotides
    are within the last 5 bases of
    the 3' end
    """
    return seq[-5:].count('g') + seq[-5:].count('c')


def gc_clamp_present(seq):
    """
    checks if a gc clamp is present
    """
    if config.PRIMER_GC_CLAMP > 0:
        for nuc in seq[-config.PRIMER_GC_CLAMP:]:
            if nuc in "cg":
                clamp_present = True
            else:
                clamp_present = False
                break
    else:
        clamp_present = True

    return clamp_present


def is_three_prime_ambiguous(amb_seq):
    """
    determine if a sequence contains an ambiguous char at the 3'prime
    """
    len_3_prime = config.PRIMER_MIN_3_WITHOUT_AMB

    if len_3_prime != 0:
        for nuc in amb_seq[len(amb_seq)-len_3_prime:]:
            if nuc not in config.nucs:
                is_amb = True
                break
            else:
                is_amb = False
    else:
        is_amb = False

    return is_amb


def rev_complement(seq):
    """
    reverse complement a sequence
    """
    return(str(Seq(seq).reverse_complement()))


def filter_kmer_direction_dependend(direction, kmer, ambiguous_consensus):
    """
    filter for 3'ambiguous, hairpin temp and end GC content.
    this differs depending on the direction of the kmer.
    """
    # get the correct kmer to test
    if direction == "LEFT":
        kmer_seq = kmer[0]
        amb_kmer_seq = ambiguous_consensus[kmer[1]:kmer[2]]
    elif direction == "RIGHT":
        kmer_seq = rev_complement(kmer[0])
        amb_kmer_seq = rev_complement(ambiguous_consensus[kmer[1]:kmer[2]])
    # filter kmer
    return(
        (calc_hairpin(kmer_seq).tm <= config.PRIMER_HAIRPIN)
        and (calc_end_gc(kmer_seq) <= config.PRIMER_MAX_GC_END)
        and gc_clamp_present(kmer_seq)
        and not is_three_prime_ambiguous(amb_kmer_seq)
    )


def hardfilter_kmers(kmers, ambiguous_consensus):
    """
    hardfilter kmers based on their seq and direction
    """

    left_primer_candidates = []
    right_primer_candidates = []

    for kmer in kmers:
        if filter_kmer_direction_independent(kmer[0]):
            base_penalty = calc_base_penalty(kmer[0])
            if base_penalty <= config.PRIMER_MAX_BASE_PENALTY:
                for direction in ["LEFT", "RIGHT"]:
                    if filter_kmer_direction_dependend(direction, kmer, ambiguous_consensus):
                        if direction == "LEFT":
                            left_primer_candidates.append(
                                [kmer[0], kmer[1], kmer[2], base_penalty]
                            )
                        if direction == "RIGHT":
                            right_primer_candidates.append(
                                [kmer[0], kmer[1], kmer[2], base_penalty]
                            )

    return left_primer_candidates, right_primer_candidates


def filter_for_lowest_scoring(direction, primer_candidates):
    """
    sort the primers by start and base penalty.
    for primers that have the same start, retain only
    the lowest scoring.
    """
    candidates = []

    # sort for start of the primer and score
    if direction == "LEFT":
        start = 1   # start index of the forward primer
    elif direction == "RIGHT":
        start = 2   # stop index as it is the start of the reverse primer
    primer_candidates.sort(key=lambda x: (x[start], x[3]))

    # pick the lowest scoring primer for primers that have the same start
    prev_start = -1
    for primer in primer_candidates:
        if primer[start] != prev_start:
            candidates.append(primer)
            prev_start = primer[start]

    return candidates


def calc_permutation_penalty(amb_seq):
    """
    get all permutations of a primer with ambiguous
    nucleotides and multiply with permutation penalty
    """
    permutations = 0

    for nuc in amb_seq:
        if nuc in config.ambig_nucs:
            n = len(config.ambig_nucs[nuc])
            if permutations != 0:
                permutations = permutations*n
            else:
                permutations = n

    return permutations*config.PRIMER_PERMUTATION_PENALTY


def get_per_base_mismatches(primer, alignment, ambiguous_consensus):
    """
    calculate for a given primer with [seq, start, stop]
    the percent mismatch per primer pos with the alignment.
    considers if primer or sequences have an amb nuc. returns
    a list of percent mismatches for each primer position.
    """
    # ini list
    mismatches = len(primer[0])*[0]
    # get primer with ambiguous nucs
    ambigous_primer = ambiguous_consensus[primer[1]:primer[2]]
    # test it against all sequences in the alignment
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
                            mismatches[idx] += 1
                    # if no amb pos is in primer then check if primer nuc
                    # is part of the amb slice nuc
                    elif current_primer_pos not in config.ambig_nucs[slice_nuc]:
                        mismatches[idx] += 1
                # check if primer has an amb pos but the current
                # slice_nuc is not part of this amb nucleotide
                elif current_primer_pos in config.ambig_nucs:
                    if slice_nuc not in config.ambig_nucs[current_primer_pos]:
                        mismatches[idx] += 1
                # mismatch
                else:
                    mismatches[idx] += 1

    # gives a percent mismatch over all positions of the primer from 5' to 3'
    mismatches = [round(x/len(alignment), 2) for x in mismatches]

    return mismatches


def calc_3_prime_penalty(direction, mismatches):
    """
    calculate the penalty for mismatches at the 3' end.
    the more mismatches are closer to the 3' end of the primer,
    the higher the penalty. uses the previously calculated
    mismatch list.
    """
    if config.PRIMER_3_PENALTY:
        if direction == "RIGHT":
            penalty = sum([m * p for m, p in zip(mismatches[0:len(config.PRIMER_3_PENALTY)], config.PRIMER_3_PENALTY)])
        elif direction == "LEFT":
            penalty = sum([m * p for m, p in zip(mismatches[::-1][0:len(config.PRIMER_3_PENALTY)], config.PRIMER_3_PENALTY)])
    else:
        penalty = 0

    return(penalty)


def find_primers(kmers, ambiguous_consensus, alignment):
    """
    filter kmers direction specific, get the lowest scoring and append penalties
    --> potential primers
    """

    # filter kmers direction specific
    left_primer_candidates, right_primer_candidates = hardfilter_kmers(kmers, ambiguous_consensus)

    for direction, primer_candidates in [("LEFT", left_primer_candidates), ("RIGHT", right_primer_candidates)]:
        # check if some kmers passed
        if primer_candidates:
            # filter for lowest scoring with the same start
            primer_candidates = filter_for_lowest_scoring(direction, primer_candidates)
            # append scores
            for primer in primer_candidates:
                primer.append(
                    get_per_base_mismatches(
                        primer,
                        alignment,
                        ambiguous_consensus
                    )
                )
            primer[3] += calc_3_prime_penalty(direction, primer[4])  # add 3 prime penalty to base penalty
            primer[3] += calc_permutation_penalty(  # also add permutation penalty
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
        best_scoring = [0, primer_candidates[0][3]]
        search_window = [primer_candidates[0][1], primer_candidates[0][2]+config.PRIMER_SIZES[1]]

        for idx, primer in enumerate(primer_candidates):
            # append primers if their start is within the current window
            if primer[1] in range(search_window[0], search_window[1]):
                # remember the index and the score
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
