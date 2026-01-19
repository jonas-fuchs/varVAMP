"""
estimate varVAMP threshold if not given
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

    How it works:

    - calculate a np array with frequencies per position for ATCG and gaps (-)
    - increment for different threshold starting at 0.1
    - for each increment calculate:
        - the lengths of conseq. nucleotides that are >= threshold
        - skip all positions with gaps >= 1 - threshold
        - this results in an array of stretch lengths
        - Then check for each stretch if it could be considered as a primer region:
            e.g. 4 ambiguous bases:
                stretches: [1,1,1,1,6,10,5,3,1,1]
                - first check for the current stretch if the sum of the next n stretches (with n being
                the n ambiguous bases + 1) is larger than the minimal length of a primer
                - this would result that the with stars marked positions would pass:
                [1,1,1,1*,6*,10*,5,3,1,1]
        - now calculate the distance between the current start of a stretch and the previous stop
        - if this is larger than the minimal amplicon length - exit the optimization
        - reset the threshold to second prior iteration to make it more robust (allow more primer regions) as
        sometimes the optimization would fail with just one iteration

    Manual optimization can be beneficial in some cases.
    """
    # set coverage to max
    args.threshold = 0.1
    # calc freqs
    highest_frequencies, deletion_frequency = calculate_frequencies(preprocessed_alignment)
    # write to log
    with open(log_file, 'a') as f:
        print(f"AUTOMATIC THRESHOLD SELECTION\n", file=f)
        print(f"-t\tmaximum non-covered region", file=f)
        # calculate distance between passing potential primer regions
        previous_stop = 0
        while args.threshold < 1:
            distances = calculate_distances(highest_frequencies, deletion_frequency, args.threshold)
            max_distance_between_passing, previous_passing = 0, 0
            # check if the distance between potential primer regions is not larger than:
            # args.opt_length - 2 * args.overlap
            for idx, dis in enumerate(distances):
                if sum(distances[idx:idx + 1 + args.n_ambig]) >= config.PRIMER_SIZES[1]:
                    # the stretch start in the gap-excluded alignment is the sum of all prior distances including the current
                    # minus the distance of the current stretch
                    stretch_start = sum(distances[:idx])
                    # then the distance between the prior stop and current start is calculated
                    current_dis = stretch_start - previous_stop
                    # and the max is updated if necessary
                    if max_distance_between_passing < current_dis:
                        max_distance_between_passing = current_dis
                    # update previous stop position
                    previous_stop = stretch_start + distances[idx]
            # write each iteration to log
            print(round(args.threshold, 2), max_distance_between_passing, sep="\t", file=f)
            # check if the distance is acceptable
            distance_threshold = args.opt_length - 2 * args.overlap if args.mode == 'tiled' else args.opt_length
            if max_distance_between_passing < distance_threshold:
                # never exceed 0.99
                if args.threshold != 0.99:
                    args.threshold += 0.01
            # or reset to the param of the two previous iterations
            else:
                args.threshold -= 0.02
                break
        print(f"Automatic parameter selection set -t {round(args.threshold, 2)} at -a {args.n_ambig}.", file=f)

    return round(args.threshold, 2)
