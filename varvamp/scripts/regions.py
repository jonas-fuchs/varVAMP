"""
finding and digesting potential primer regions.
"""

# varVAMP
from varvamp.scripts import config


def find_regions(consensus_amb, allowed_ambiguous):
    """
    finds primer regions as specified by a
    certain amount of ambiguous bases in a given
    sequence length
    """

    current_window = []
    primer_regions = []
    last_amb = 0

    # append one N so the last window is closed
    seq = str(consensus_amb) + "N"

    for idx, nuc in enumerate(seq):
        # open a new window if none is there and we have a nucleotide
        if nuc in config.NUCS and not current_window:
            current_window = [idx, 0]
            all_previous_ambiguous_pos = []
        # if we have a non-nucleotide, check if its time to writ the window
        elif nuc not in config.NUCS and current_window:
            # are the current and last ambig char sufficiently far?
            if idx - last_amb >= config.PRIMER_SIZES[0] and nuc != "N":
                current_window[1] = 1
                all_previous_ambiguous_pos = [idx]
            else:
                current_window[1] += 1
                all_previous_ambiguous_pos.append(idx)
            # check if there were too many previous amb char in subwindow
            # or if N is reached -> writeable
            if current_window[1] > allowed_ambiguous or nuc == "N":
                # check if the writable window has a sufficient length.
                if idx-current_window[0] >= config.PRIMER_SIZES[0]:
                    primer_regions.append([current_window[0], idx])
                # open new window if N is reached
                if nuc == "N":
                    current_window = []
                # or reset to a previous ambig char
                else:
                    current_window[0] = all_previous_ambiguous_pos[0]+1
                    current_window[1] = current_window[1] - 1
                    all_previous_ambiguous_pos.pop(0)
            last_amb = idx

    return primer_regions


def mean(primer_regions, consensus):
    """
    calculate the percentage of regions
    that cover the consensus sequence
    """
    covered_set = set()
    for region in primer_regions:
        pos_list = list(range(region[0], region[1]))
        [covered_set.add(x) for x in pos_list]
    return round(len(covered_set)/len(consensus)*100, 1)


def digest_seq(seq, kmer_size):
    """
    digest the sequence into kmers
    """
    return[[seq[i:i+kmer_size], i, i+len(seq[i:i+kmer_size])] for i in range(len(seq)-kmer_size+1)]


def produce_kmers(primer_regions, consensus, sizes=config.PRIMER_SIZES):
    """
    produce kmers for all primer regions
    """
    kmers = []

    for region in primer_regions:
        sliced_seq = consensus[region[0]:region[1]]
        for kmer_size in range(sizes[0], sizes[1]+1):
            kmers_temp = digest_seq(sliced_seq, kmer_size)
            # adjust the start and stop position of the kmers
            for kmer_temp in kmers_temp:
                kmer_temp[1] = kmer_temp[1]+region[0]
                kmer_temp[2] = kmer_temp[2]+region[0]
                # check if kmer is already in list (overlapping regions)
                if kmer_temp in kmers:
                    continue
                kmers.append(kmer_temp)

    return kmers
