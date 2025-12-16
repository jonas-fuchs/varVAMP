"""
estimate varVAMP threshold and max n of ambiguous bases if none are given
"""

# libs
import numpy as np

# varVAMP
from varvamp.scripts import config


def calculate_frequencies(preprocessed_alignment):
    """
    calculate the max nucleotide freq at each pos
    """
    all_freqs = []

    for i in range(0, len(preprocessed_alignment[0][1])):
        nuc_dict = {
            "a": 0,
            "c": 0,
            "t": 0,
            "g": 0,
        }
        for seq in preprocessed_alignment:
            if seq[1][i] in nuc_dict:
                nuc_dict[seq[1][i]] += 1
        highest_freq = max(nuc_dict.values())
        all_freqs.append(highest_freq/len(preprocessed_alignment))

    return all_freqs


def calculate_distances(all_freqs, threshold):
    """
    calc the distance for each nuc freq to the prior
    nuc freq that fell below the cutoff
    """

    current_dis = 0
    previous_dis = 0
    all_dis = []

    for idx, freq in enumerate(all_freqs):
        if freq < threshold:
            current_dis = 0
        current_dis += 1
        if current_dis <= previous_dis or idx == len(all_freqs)-1:
            all_dis.append(previous_dis)

        previous_dis = current_dis

    return all_dis


def get_parameters(preprocessed_alignment, args, log_file):
    """
    give an estimate for number of ambiguous chars and/or threshold
    writes to log file
    """
    # set coverage to max
    args.threshold = 0.1
    # read in the alignment and calc freqs
    frequencies = calculate_frequencies(preprocessed_alignment)
    # write to log
    with open(log_file, 'a') as f:
        print(f"AUTOMATIC THRESHOLD SELECTION\n", file=f)
        print(f"-t\tmaximum non-covered region", file=f)
        # calculate distance between passing potential primer regions
        while args.threshold < 1:
            distances = calculate_distances(frequencies, args.threshold)
            max_distance_between_passing, previous_passing = 0, 0
            for idx, dis in enumerate(distances):
                if sum(distances[idx:idx + 1 + args.n_ambig]) >= config.PRIMER_SIZES[1]:
                    current_dis = idx - previous_passing
                    if max_distance_between_passing < current_dis:
                        max_distance_between_passing = current_dis
                    previous_passing = idx
            # consider the threshold in which the distance between passing regions is smaller than the optimal primer
            if max_distance_between_passing < args.max_length:
                # write each iteration to log
                print(round(args.threshold, 2), max_distance_between_passing, sep="\t", file=f)
                args.threshold += 0.005
            # or reset to the param of the prior iteration
            else:
                args.threshold -= 0.005
                break
        print(f"Automatic parameter selection set -t {round(args.threshold, 3)} at -a {args.n_ambig}.", file=f)

    return round(args.threshold, 2)
