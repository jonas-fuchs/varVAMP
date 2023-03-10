"""
consensus creation
"""

# BUILT-INS
import collections

# varVAMP
from varvamp.scripts import config


def determine_nucleotide_counts(alignment, idx):
    """
    count the number of each nucleotides at
    an idx of the alignment. return sorted dic.
    handels ambiguous nucleotides in sequences.
    also handels gaps.
    """
    nucleotide_list = []

    # get all nucleotides
    for sequence in alignment:
        nucleotide_list.append(sequence[1][idx])
    # count occurences of nucleotides
    counter = dict(collections.Counter(nucleotide_list))
    # get permutations of an ambiguous nucleotide
    to_delete = []
    temp_dict = {}
    for nucleotide in counter:
        if nucleotide in config.ambig_nucs:
            to_delete.append(nucleotide)
            permutations = config.ambig_nucs[nucleotide]
            adjusted_freq = 1/len(permutations)
            for permutation in permutations:
                if permutation in temp_dict:
                    temp_dict[permutation] += adjusted_freq
                else:
                    temp_dict[permutation] = adjusted_freq
        if nucleotide == "-":
            to_delete.append(nucleotide)

    # drop ambiguous entrys and add adjusted freqs to
    if to_delete:
        for i in to_delete:
            counter.pop(i)
        for nucleotide in temp_dict:
            if nucleotide in counter:
                counter[nucleotide] += temp_dict[nucleotide]
            else:
                counter[nucleotide] = temp_dict[nucleotide]

    return dict(sorted(counter.items(), key=lambda x: x[1], reverse=True))


def get_consensus_nucleotides(nucleotide_counts, consensus_cutoff):
    """
    get a list of nucleotides for the consensus seq
    """
    n = 0

    consensus_nucleotides = []
    for nuc in nucleotide_counts:
        n += nucleotide_counts[nuc]
        consensus_nucleotides.append(nuc)
        if n >= consensus_cutoff:
            break

    return consensus_nucleotides


def get_ambiguous_char(nucleotides):
    """
    get ambiguous char from a list of nucleotides
    """
    for ambiguous, permutations in config.ambig_nucs.items():
        if set(permutations) == set(nucleotides):
            return ambiguous


def create_consensus(alignment, threshold):
    """
    build a majority sequence and a sequence that
    has ambiguous chars as determined by the freq
    threshold.
    """

    # ini the consensus seq
    ambiguous_consensus = str()
    majority_consensus = str()

    # define consensus cut-off
    consensus_cutoff = len(alignment)*threshold
    # define length of the consensus from the first seq in alignment
    length_consensus = len(alignment[0][1])

    # built consensus sequences
    for idx in range(length_consensus):
        nucleotide_counts = determine_nucleotide_counts(alignment, idx)
        consensus_nucleotide = get_consensus_nucleotides(
            nucleotide_counts,
            consensus_cutoff
        )
        if len(consensus_nucleotide) > 1:
            amb_consensus_nucleotide = get_ambiguous_char(consensus_nucleotide)
            ambiguous_consensus = ambiguous_consensus + amb_consensus_nucleotide
        else:
            ambiguous_consensus = ambiguous_consensus + consensus_nucleotide[0]

        majority_consensus = majority_consensus + consensus_nucleotide[0]

    return majority_consensus, ambiguous_consensus
