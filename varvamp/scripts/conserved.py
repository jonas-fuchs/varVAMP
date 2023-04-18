"""
finding and digesting conserved regions.
"""

# varVAMP
from varvamp.scripts import config


def find_regions(consensus_amb, allowed_ambiguous):
    """
    finds conserved regions as specified by a
    certain amount of ambiguous bases in a given
    sequence length
    """
    # init the variables
    current_window = []
    writable = False
    in_ambiguous_region = True
    last_amb = 0
    conserved_regions = []

    # append enough Ns to ensure that last window is closed
    seq = str(consensus_amb) + allowed_ambiguous*'N'

    for idx, nuc in enumerate(seq):
        if in_ambiguous_region and nuc in config.NUCS:
            in_ambiguous_region = False
            # just entered a new stretch of non-ambiguous bases
            # may be time to open a new window
            if not current_window:
                # track window start and ambiguous chars so far
                current_window = [idx, 0]
                all_previous_ambiguous_pos = []
            continue
        if nuc not in config.NUCS:
            if current_window:
                in_ambiguous_region = True
                amb_to_amb_len = idx - last_amb
                all_previous_ambiguous_pos.append(idx)
                # check if there were too many previous amb char in subwindow
                # or if N is reached -> writeable
                if current_window[1] == allowed_ambiguous or nuc == "N":
                    writable = True
                    # are the current and last ambig char sufficiently far?
                    if amb_to_amb_len >= config.PRIMER_SIZES[0] and nuc != "N":
                        current_window[1] = 0
                        all_previous_ambiguous_pos = [idx]
                        writable = False

                current_window[1] += 1

                if writable:
                    writable = False
                    # check if the writable window has a sufficient length.
                    window_length = idx-current_window[0]
                    if window_length >= config.PRIMER_SIZES[0]:
                        conserved_regions.append([current_window[0], idx])
                    # reset the window start to the last ambig char
                    if nuc != "N":
                        current_window[0] = all_previous_ambiguous_pos[0]+1
                        current_window[1] = current_window[1] - 1
                        all_previous_ambiguous_pos.pop(0)
                    # or open a new window
                    else:
                        current_window = []

            last_amb = idx

    return conserved_regions


def mean(conserved_regions, consensus):
    """
    calculate the percentage of regions
    that are conserved
    """
    covered_set = set()
    for region in conserved_regions:
        pos_list = list(range(region[0], region[1]))
        [covered_set.add(x) for x in pos_list]
    return round(len(covered_set)/len(consensus)*100, 1)


def digest_seq(seq, kmer_size):
    """
    digest the sequence into kmers
    """
    return[[seq[i:i+kmer_size], i, i+len(seq[i:i+kmer_size])] for i in range(len(seq)-kmer_size+1)]


def produce_kmers(conserved_regions, consensus, sizes=config.PRIMER_SIZES):
    """
    produce kmers for all conserved regions
    """
    kmers = []

    for region in conserved_regions:
        sliced_seq = consensus[region[0]:region[1]]
        for kmer_size in range(sizes[0], sizes[1]+1):
            kmers_temp = digest_seq(sliced_seq, kmer_size)
            # adjust the start and stop position of the kmers
            for kmer_temp in kmers_temp:
                kmer_temp[1] = kmer_temp[1]+region[0]
                kmer_temp[2] = kmer_temp[2]+region[0]
                if kmer_temp in kmers:
                    continue
                kmers.append(kmer_temp)

    return kmers
