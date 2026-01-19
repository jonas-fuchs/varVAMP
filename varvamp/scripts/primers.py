"""
primer creation and evaluation
"""

# BUILTIN
from itertools import chain

# LIBS
from Bio.Seq import Seq
import primer3 as p3

# varVAMP
from varvamp.scripts import config, reporting


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


def calc_dimer(seq1, seq2, structure=False):
    """
    Calculate the hetero-dimerization thermodynamics of two DNA sequences.
    Return primer3 thermo object.
    """
    return p3.calc_heterodimer(
        seq1.upper(),
        seq2.upper(),
        mv_conc=config.PCR_MV_CONC,
        dv_conc=config.PCR_DV_CONC,
        dna_conc=config.PCR_DNA_CONC,
        dntp_conc=config.PCR_DNTP_CONC,
        output_structure=structure
    )


def calc_max_polyx(seq):
    """
    calculate maximum polyx of a seq
    """
    previous_nuc, counter, max_polyx = seq[0], 1, 1

    for nuc in seq[1:]:
        if nuc == previous_nuc:
            counter += 1
        else:
            counter = 1
            previous_nuc = nuc
        if counter > max_polyx:
            max_polyx = counter

    return max_polyx


def calc_max_dinuc_repeats(seq):
    """
    calculate maximum polyxy of a seq
    """
    for s in [seq, seq[1:]]:
        previous_dinuc = s[0:2]
        max_dinuc = 1
        counter = 1
        for i in range(2, len(s), 2):
            if s[i:i+2] == previous_dinuc:
                counter += 1
            else:
                if counter > max_dinuc:
                    max_dinuc = counter
                counter = 1
                previous_dinuc = s[i:i+2]

    return max_dinuc


def calc_end_gc(seq):
    """
    check how many gc nucleotides
    are within the last 5 bases of
    the 3' end
    """
    return seq[-5:].count('g') + seq[-5:].count('c')


def is_three_prime_ambiguous(amb_seq):
    """
    determine if a sequence contains an ambiguous char at the 3'prime
    """
    len_3_prime, is_amb = config.PRIMER_MIN_3_WITHOUT_AMB, False

    if len_3_prime != 0:
        for nuc in amb_seq[len(amb_seq)-len_3_prime:]:
            if nuc not in config.NUCS:
                is_amb = True
                break

    return is_amb


def rev_complement(seq):
    """
    reverse complement a sequence
    """
    return str(Seq(seq).reverse_complement())


def calc_permutation_penalty(amb_seq):
    """
    get all permutations of a primer with ambiguous
    nucleotides and multiply with permutation penalty
    """
    permutations = 0

    for nuc in amb_seq:
        if nuc in config.AMBIG_NUCS:
            n = len(config.AMBIG_NUCS[nuc])
            if permutations != 0:
                permutations = permutations*n
            else:
                permutations = n

    return permutations*config.PRIMER_PERMUTATION_PENALTY


def calc_base_penalty(seq, primer_temps, gc_range, primer_sizes):
    """
    Calculate intrinsic primer penalty.
    """
    penalty = 0

    tm = calc_temp(seq)
    gc = calc_gc(seq)
    size = len(seq)

    # TEMP penalty
    if tm > primer_temps[2]:
        penalty += config.PRIMER_TM_PENALTY*(
            tm - primer_temps[2]
            )
    if tm < primer_temps[2]:
        penalty += config.PRIMER_TM_PENALTY*(
            primer_temps[2] - tm
            )
    # GC penalty
    if gc > gc_range[2]:
        penalty += config.PRIMER_GC_PENALTY*(
            gc - gc_range[2]
            )
    if gc < gc_range[2]:
        penalty += config.PRIMER_GC_PENALTY*(
            gc_range[2] - gc
        )
    # SIZE penalty
    if size > primer_sizes[2]:
        penalty += config.PRIMER_SIZE_PENALTY*(
            size - primer_sizes[2]
        )
    if size < primer_sizes[2]:
        penalty += config.PRIMER_SIZE_PENALTY * (
            primer_sizes[2] - size
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
            if slice_nuc in config.AMBIG_NUCS:
                # check if the kmer has an amb pos
                if current_kmer_pos in config.AMBIG_NUCS:
                    slice_nuc_set = set(config.AMBIG_NUCS[slice_nuc])
                    pri_set = set(config.AMBIG_NUCS[current_kmer_pos])
                    # check if these sets have no overlap
                    # -> mismatch
                    if len(slice_nuc_set.intersection(pri_set)) == 0:
                        mismatches[idx] += 1
                # if no amb pos is in kmer then check if kmer nuc
                # is part of the amb slice nuc
                elif current_kmer_pos not in config.AMBIG_NUCS[slice_nuc]:
                    mismatches[idx] += 1
            # check if kmer has an amb pos but the current
            # slice_nuc is not part of this amb nucleotide
            elif current_kmer_pos in config.AMBIG_NUCS:
                if slice_nuc not in config.AMBIG_NUCS[current_kmer_pos]:
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

    penalty = 0

    if config.PRIMER_3_PENALTY:
        if direction == "-":
            penalty = sum([m * p for m, p in zip(mismatches[0:len(config.PRIMER_3_PENALTY)], config.PRIMER_3_PENALTY)])
        elif direction == "+":
            penalty = sum([m * p for m, p in zip(mismatches[::-1][0:len(config.PRIMER_3_PENALTY)], config.PRIMER_3_PENALTY)])

    return penalty


def filter_kmer_direction_independent(seq, primer_temps=config.PRIMER_TMP, gc_range=config.PRIMER_GC_RANGE, primer_sizes=config.PRIMER_SIZES):
    """
    filter kmer for temperature, gc content,
    poly x, dinucleotide repeats and homodimerization
    """
    return(
        (primer_temps[0] <= calc_temp(seq) <= primer_temps[1])
        and (gc_range[0] <= calc_gc(seq) <= gc_range[1])
        and (calc_max_polyx(seq) <= config.PRIMER_MAX_POLYX)
        and (calc_max_dinuc_repeats(seq) <= config.PRIMER_MAX_DINUC_REPEATS)
        and (calc_base_penalty(seq, primer_temps, gc_range, primer_sizes) <= config.PRIMER_MAX_BASE_PENALTY)
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
    return (
        (calc_hairpin(kmer_seq).tm <= config.PRIMER_HAIRPIN)
        and (config.PRIMER_GC_END[0] <= calc_end_gc(kmer_seq) <= config.PRIMER_GC_END[1])
        and not is_three_prime_ambiguous(amb_kmer_seq)
    )


def parse_primer_fasta(fasta_path):
    """
    Parse a primer FASTA file and return a list of sequences.
    """
    sequences = []
    current_seq = []

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    seq = ''.join(current_seq).lower()
                    if len(current_seq) <= 40:
                        sequences.append(reporting.get_permutations(seq))
                current_seq = []
            else:
                current_seq.append(line)
        # read last sequence
        if current_seq:
            seq = ''.join(current_seq).lower()
            if len(current_seq) <= 40:
                sequences.append(reporting.get_permutations(seq))

    return list(chain.from_iterable(sequences))


def check_dimer_with_sequences(primer_seq, external_sequences):
    """
    Check if a primer forms dimers with any of the provided sequences.
    Considers all permutations if primer contains degenerate bases.
    Returns True if a dimer is formed above the threshold.
    """


    for seq in external_sequences:
        if calc_dimer(primer_seq, seq).tm > config.PRIMER_MAX_DIMER_TMP:
            return True

    return False


def filter_non_dimer_candidates(primer_candidates, external_sequences):
    """
    Filter out primer candidates that form dimers with external sequences.
    """
    filtered = []
    for primer in primer_candidates:
        if not check_dimer_with_sequences(primer[0], external_sequences):
            filtered.append(primer)

    return filtered


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
        base_penalty = calc_base_penalty(kmer[0],config.PRIMER_TMP, config.PRIMER_GC_RANGE, config.PRIMER_SIZES)
        # calculate per base mismatches
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
        primer_name = f"{direction_name}_{primer_idx}"
        primer_dict[primer_name] = primer
        primer_idx += 1

    return primer_dict


def find_best_primers(left_primer_candidates, right_primer_candidates):
    """
    Primer candidates are likely overlapping. Here, the list of primers
    is sorted for the lowest to highest penalty. Then, the next lowest
    is retained if it does not have any nucleotides that have already
    been covered by the middle third of a better primer candidate.
    This significantly reduces the amount of primers while retaining
    the ones with the lowest penalties.

    Example:
    -------- (penalty 1) 1
        ------------- (penalty 1) 2
            ------------ (penalty 0.8) 3
                      ----------(penalty 0.9) 4

    --> primer 3 would be retained and primer 2 excluded, primer 4 and 1
    will however be considered in the next set of overlapping primers.

    """
    all_primers = {}

    for direction, primer_candidates in [("+", left_primer_candidates), ("-", right_primer_candidates)]:
        # sort the primers by penalty, and if same penalty by start
        primer_candidates.sort(key=lambda x: (x[3], x[1]))
        # ini everything with the primer with the lowest penalty
        to_retain = [primer_candidates[0]]
        primer_ranges = list(range(primer_candidates[0][1], primer_candidates[0][2]))
        primer_set = set(primer_ranges)

        for primer in primer_candidates:
            # get the thirds of the primer, only consider the middle
            thirds_len = int((primer[2] - primer[1])/3)
            primer_positions = list(range(primer[1] + thirds_len, primer[2] - thirds_len))
            # check if none of the nucleotides of the next primer
            # are already covered by a better primer
            if not any(x in primer_positions for x in primer_set):
                # update the primer set
                primer_set.update(primer_positions)
                # append this primer as it has a low penalty and is not overlapping
                # with another already retained primer
                to_retain.append(primer)

        # sort by start
        to_retain.sort(key=lambda x: x[1])
        # create dict
        all_primers[direction] = create_primer_dictionary(to_retain, direction)

    # and create a dict
    return all_primers
