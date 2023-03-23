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


def calc_per_base_mismatches(kmer, alignment, ambiguous_consensus):
    """
    calculate for a given kmer with [seq, start, stop]
    the percent mismatch per kmer pos with the alignment.
    considers if kmer or aln sequences have an amb nuc. returns
    a list of percent mismatches for each kmer position.
    """
    # ini list
    mismatches = len(kmer[0])*[0]
    # get kmer with ambiguous nucs
    amb_kmer = ambiguous_consensus[kmer[1]:kmer[2]]
    # test it against all sequences in the alignment
    for sequence in alignment:
        # slice each sequence of the alignment for the kmer
        # start and stop positions
        seq_slice = sequence[1][kmer[1]:kmer[2]]
        for idx, slice_nuc in enumerate(seq_slice):
            # find the respective nuc to that of the slice
            current_kmer_pos = amb_kmer[idx]
            if slice_nuc == current_kmer_pos:
                continue
            # check if the slice nucleotide is an amb pos
            if slice_nuc in config.ambig_nucs:
                # check if the kmer has an amb pos
                if current_kmer_pos in config.ambig_nucs:
                    slice_nuc_set = set(config.ambig_nucs[slice_nuc])
                    pri_set = set(config.ambig_nucs[current_kmer_pos])
                    # check if these sets have no overlap
                    # -> mismatch
                    if len(slice_nuc_set.intersection(pri_set)) == 0:
                        mismatches[idx] += 1
                # if no amb pos is in kmer then check if kmer nuc
                # is part of the amb slice nuc
                elif current_kmer_pos not in config.ambig_nucs[slice_nuc]:
                    mismatches[idx] += 1
            # check if kmer has an amb pos but the current
            # slice_nuc is not part of this amb nucleotide
            elif current_kmer_pos in config.ambig_nucs:
                if slice_nuc not in config.ambig_nucs[current_kmer_pos]:
                    mismatches[idx] += 1
            # mismatch
            else:
                mismatches[idx] += 1

    # gives a percent mismatch over all positions of the kmer from 5' to 3'
    mismatches = [round(x/len(alignment), 2) for x in mismatches]

    return mismatches


def calc_3_prime_penalty(direction, mismatches):
    """
    calculate the penalty for mismatches at the 3' end.
    the more mismatches are closer to the 3' end of the kmer,
    the higher the penalty. uses the previously calculated
    mismatch list.
    """
    if config.PRIMER_3_PENALTY:
        if direction == "-":
            penalty = sum([m * p for m, p in zip(mismatches[0:len(config.PRIMER_3_PENALTY)], config.PRIMER_3_PENALTY)])
        elif direction == "+":
            penalty = sum([m * p for m, p in zip(mismatches[::-1][0:len(config.PRIMER_3_PENALTY)], config.PRIMER_3_PENALTY)])
    else:
        penalty = 0

    return(penalty)


def filter_kmer_direction_independent(seq):
    """
    filter kmer for temperature, gc content,
    poly x, dinucleotide repeats and homodimerization
    """
    return(
        (config.PRIMER_TMP[0] <= calc_temp(seq) <= config.PRIMER_TMP[1])
        and (config.PRIMER_GC_RANGE[0] <= calc_gc(seq) <= config.PRIMER_GC_RANGE[1])
        and (calc_max_polyx(seq) <= config.PRIMER_MAX_POLYX)
        and (calc_max_dinuc_repeats(seq) <= config.PRIMER_MAX_DINUC_REPEATS)
        and (calc_base_penalty(seq) <= config.PRIMER_MAX_BASE_PENALTY)
        and (calc_dimer(seq, seq).tm <= config.PRIMER_MAX_DIMER_TMP)
    )


def filter_kmer_direction_dependend(direction, kmer, ambiguous_consensus):
    """
    filter for 3'ambiguous, hairpin temp and end GC content.
    this differs depending on the direction of the kmer.
    """
    # get the correct amb kmer to test
    if direction == "+":
        kmer_seq = kmer[0]
        amb_kmer_seq = ambiguous_consensus[kmer[1]:kmer[2]]
    elif direction == "-":
        kmer_seq = rev_complement(kmer[0])
        amb_kmer_seq = rev_complement(ambiguous_consensus[kmer[1]:kmer[2]])
    # filter kmer
    return(
        (calc_hairpin(kmer_seq).tm <= config.PRIMER_HAIRPIN)
        and (calc_end_gc(kmer_seq) <= config.PRIMER_MAX_GC_END)
        and gc_clamp_present(kmer_seq)
        and not is_three_prime_ambiguous(amb_kmer_seq)
    )


def find_primers(kmers, ambiguous_consensus, alignment):
    """
    filter kmers direction specific and append penalties
    --> potential primers
    """
    left_primer_candidates = []
    right_primer_candidates = []

    for kmer in kmers:
        # filter kmers based on their direction independend stats
        if not filter_kmer_direction_independent(kmer[0]):
            continue
        # calc base penalty
        base_penalty = calc_base_penalty(kmer[0])
        # calcualte per base mismatches
        per_base_mismatches = calc_per_base_mismatches(
                                kmer,
                                alignment,
                                ambiguous_consensus
                            )
        # calculate permutation penealty
        permutation_penalty = calc_permutation_penalty(
                                ambiguous_consensus[kmer[1]:kmer[2]]
                            )
        # now check direction specific
        for direction in ["+", "-"]:
            # check if kmer passes direction filter
            if not filter_kmer_direction_dependend(direction, kmer, ambiguous_consensus):
                continue
            # calculate the 3' penalty
            three_prime_penalty = calc_3_prime_penalty(
                                    direction,
                                    per_base_mismatches
                                )
            # add all penalties
            primer_penalty = base_penalty + permutation_penalty + three_prime_penalty
            # sort into lists
            if direction == "+":
                left_primer_candidates.append(
                    [kmer[0], kmer[1], kmer[2], primer_penalty, per_base_mismatches]
                )
            if direction == "-":
                right_primer_candidates.append(
                    [rev_complement(kmer[0]), kmer[1], kmer[2], primer_penalty, per_base_mismatches]
                )

    return left_primer_candidates, right_primer_candidates


def create_primer_dictionary(primer_candidates, direction):
    """
    creates a primer dictionary from primer list
    """
    primer_dict = {}
    primer_idx = 0

    for primer in primer_candidates:
        if direction == "+":
            direction_name = "LEFT"
        elif direction == "-":
            direction_name = "RIGHT"
        primer_name = direction_name + "_" + str(primer_idx)
        primer_dict[primer_name] = primer
        primer_idx += 1

    return primer_dict


def find_best_primers(left_primer_candidates, right_primer_candidates):
    """
    Primer candidates are likely overlapping. Here, the list of primers
    is sorted for the best to worst scoring. Then, the next best scoring
    is retained if it does not have any nucleotides that have already
    been covered by a better scoring primer candidate. This significantly
    reduces the amount of primers while retaining the best scoring ones.
    """
    all_primers = {}

    for direction, primer_candidates in [("+", left_primer_candidates), ("-", right_primer_candidates)]:
        # sort the primers for the best scoring
        primer_candidates.sort(key=lambda x: x[3])
        # ini everything with the top scoring primer
        to_retain = [primer_candidates[0]]
        primer_ranges = list(range(primer_candidates[0][1], primer_candidates[0][2]+1))
        primer_set = set(primer_ranges)

        for primer in primer_candidates:
            primer_positions = list(range(primer[1], primer[2]+1))
            # check if none of the nucleotides of the next primer
            # are already covered by a better primer
            if not any(x in primer_positions for x in primer_set):
                # update the primer set
                primer_set.update(primer_positions)
                # append this primer as it is well scoring and not overlapping
                # with another already retained primer
                to_retain.append(primer)

        # sort by start
        to_retain.sort(key=lambda x: x[1])
        # create dict
        all_primers[direction] = create_primer_dictionary(to_retain, direction)

    # and create a dict
    return all_primers
