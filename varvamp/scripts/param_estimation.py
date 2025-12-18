"""
estimate varVAMP threshold of ambiguous bases if none are given
"""

# libs
import numpy as np

# varVAMP
from varvamp.scripts import config


def calculate_frequencies(preprocessed_alignment):
    """
    calculate individual frequencies for a, t, c, g, and gaps at each position
    returns a 2D numpy array where each column is a position and rows are [a, c, t, g, -]
    """
    # convert alignment to numpy array
    alignment_array = np.array([list(seq.lower()) for _, seq in preprocessed_alignment])
    num_sequences, sequence_length = alignment_array.shape

    # calculate frequencies
    frequencies = np.zeros((5, sequence_length))
    nucleotides = ['a', 'c', 't', 'g', '-']

    for i in range(sequence_length):
        column = alignment_array[:, i]
        for nuc_idx, nuc in enumerate(nucleotides):
            frequencies[nuc_idx, i] = np.count_nonzero(column == nuc) / num_sequences

    return np.max(frequencies[0:4, :], axis=0), frequencies[4, :]


def calculate_distances(highest_freq, deletion_freq, threshold):
    """
    calc the distance for each nuc freq to the prior
    nuc freq that fell below the cutoff
    """

    current_dis = 0
    previous_dis = 0
    all_dis = []

    for idx, (freq, del_freq) in enumerate(zip(highest_freq, deletion_freq)):
        if freq < threshold:
            current_dis = 0
        # ignore gaps -> they will be excluded later and might skew distance calc
        if del_freq > 1-threshold:
            continue
        current_dis += 1
        if current_dis <= previous_dis or idx == len(highest_freq)-1:
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
    highest_frequencies, deletion_frequency = calculate_frequencies(preprocessed_alignment)
    # write to log
    with open(log_file, 'a') as f:
        print(f"AUTOMATIC THRESHOLD SELECTION\n", file=f)
        print(f"-t\tmaximum non-covered region", file=f)
        # calculate distance between passing potential primer regions
        while args.threshold < 1:
            distances = calculate_distances(highest_frequencies, deletion_frequency, args.threshold)
            max_distance_between_passing, previous_passing = 0, 0
            # check if the distance between potential primer regions is not larger than:
            #
            for idx, dis in enumerate(distances):
                if sum(distances[idx:idx + 1 + args.n_ambig]) >= config.PRIMER_SIZES[1]:
                    original_index = sum(distances[:idx+1])
                    current_dis = original_index - distances[idx] - previous_passing
                    if max_distance_between_passing < current_dis:
                        max_distance_between_passing = current_dis
                    print(current_dis)
                    previous_passing = original_index
            # check if the distance is acceptable
            if max_distance_between_passing < args.opt_length/2:
                # write each iteration to log
                print(round(args.threshold, 3), max_distance_between_passing, sep="\t", file=f)
                args.threshold += 0.005
            # or reset to the param of the prior iteration
            else:
                args.threshold -= 0.005
                break
        print(f"Automatic parameter selection set -t {round(args.threshold, 3)} at -a {args.n_ambig}.", file=f)

    return round(args.threshold, 2)
