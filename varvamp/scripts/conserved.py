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

    seq = str(consensus_amb) + 2*'N'
    for idx, nuc in enumerate(seq):
        if in_ambiguous_region and nuc in config.nucs:
            in_ambiguous_region = False
            # just entered a new stretch of non-ambiguous bases
            # may be time to open a new window
            if not current_window:
                current_window = [idx, 0]
                amb_pos = []
                # create new window if none is there. First element
                # keeps track of start of the window, second element is
                # a counter that resets if two ambiguous chars are longer
                # than specified apart and last one counts all ambiguous
                # chars. also track all amb chars after a window has opened
            continue
        if nuc not in config.nucs:
            if current_window:
                in_ambiguous_region = True
                amb_to_amb_len = idx - last_amb
                if nuc != "N":
                    # track previous amb pos only if current pos is not a N as this
                    # region is witeable
                    amb_pos.append(idx)
                if current_window[1] >= allowed_ambiguous or nuc == "N":
                    # check if there were too many previous amb char in subwindow
                    # and make it writable. Always make it writeable if N is
                    # reached
                    writable = True
                    if amb_to_amb_len >= config.PRIMER_SIZES[0] and nuc != "N":
                        # check if the last amb is sufficiently far, if yes keep
                        # window open and set amb counter to 0, reset also the
                        # list of amb positions and track only the current pos
                        current_window[1] = 0
                        writable = False
                        amb_pos = [idx]

                current_window[1] += 1

                if writable:
                    writable = False
                    window_length = idx-current_window[0]
                    if window_length >= config.PRIMER_SIZES[0]:
                        # check if the writable window has a sufficient length.
                        conserved_regions.append([current_window[0], idx])
                        # reset the window and the list of amb positions
                        # after it was written
                        current_window = []
                    elif nuc == "N":
                        # if nuc was a N and region was not written also open a
                        # new window
                        current_window = []
                    else:
                        # else set the start pos to the next amb pos and
                        # check again if the new window matches the criteria
                        current_window[0] = amb_pos[0]+1
                        current_window[1] = current_window[1]-1
                        amb_pos.pop(0)
            last_amb = idx

    return conserved_regions


def mean(conserved_regions, consensus):
    """
    calculate the percentage of regions
    that are conserved
    """
    sum = 0
    for region in conserved_regions:
        sum += region[1]-region[0]
    return round(sum/len(consensus)*100, 1)


def digest_seq(seq, kmer_size):
    """
    digest the sequence into kmers
    """
    return[[seq[i:i+kmer_size], i, i+len(seq[i:i+kmer_size])] for i in range(len(seq)-kmer_size+1)]


def produce_kmers(conserved_regions, consensus):
    """
    produce kmers for all conserved regions
    """
    kmers = []

    for region in conserved_regions:
        sliced_seq = consensus[region[0]:region[1]]
        for kmer_size in range(config.PRIMER_SIZES[0], config.PRIMER_SIZES[1]+1):
            kmers_temp = digest_seq(sliced_seq, kmer_size)
            # adjust the start and stop position of the kmers
            for kmer_temp in kmers_temp:
                kmer_temp[1] = kmer_temp[1]+region[0]
                kmer_temp[2] = kmer_temp[2]+region[0]
            kmers += kmers_temp

    return kmers
