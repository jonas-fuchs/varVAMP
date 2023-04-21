"""
estimate varVAMP threshold and max n of ambiguous bases if none are given
"""

# varVAMP
from varvamp.scripts import config
from varvamp.scripts import alignment


def calculate_frequencys(aln):
    """
    calculate the max nucleotide freq at each pos
    """
    all_freqs = []

    for i in range(0, len(aln[0][1])):
        nuc_dict = {
            "a": 0,
            "c": 0,
            "t": 0,
            "g": 0,
        }
        for seq in aln:
            if seq[1][i].lower() in nuc_dict:
                nuc_dict[seq[1][i].lower()] += 1
        highest_freq = max(nuc_dict.values())
        all_freqs.append(highest_freq/len(aln))

    return all_freqs


def calculate_distances(all_freqs, threshold):
    """
    calc the distance for each freq to the prior
    freq that fell below the cutoff
    """

    current_dis = 0
    previous_dis = 0
    all_dis = []

    for freq in all_freqs:
        if freq < threshold:
            current_dis = 0
        current_dis += 1
        if current_dis <= previous_dis:
            all_dis.append(previous_dis)
        previous_dis = current_dis

    return all_dis


def get_parameters(alignment_path, threshold=None, n_ambig=None):
    """
    give an estimate for number of ambiguous chars and/or threshold
    """
    # set coverage to max
    coverage = 1
    # read in the alignment and calc freqs
    aln = alignment.read_alignment(alignment_path)
    frequencys = calculate_frequencys(aln)
    # if no args for both threshold and n_ambig are given
    # set the n_ambig to 2 and optimize threshold
    if threshold is None:
        threshold = 0.1
        if n_ambig is None:
            n_ambig = 2
        fixed = False
    # if threshold is given, optimize n_ambig (number of ambiguous chars)
    else:
        n_ambig = config.PRIMER_SIZES[0]
        fixed = True
    # optimize until less than 50 % is covered
    while coverage >= 0.5 and threshold < 1:
        distances = calculate_distances(frequencys, threshold)
        # calculate the cummulative sum of the sum of n conseq. streches
        # that are together larger than the min primer length
        covered_pos = sum(
            [distances[x] for x in range(0, len(distances)) if sum(distances[x:x+n_ambig+1]) >= config.PRIMER_SIZES[0]]
        )
        # calculate coverage
        coverage = covered_pos/len(aln[0][1])
        # change the non fixed param if threshold has not been reached
        if coverage >= 0.5:
            if fixed:
                n_ambig -= 1
            else:
                threshold += 0.01
        # or reset to the param of the prior iteration
        else:
            if fixed:
                n_ambig += 1
            else:
                threshold -= 0.01

    return round(threshold, 2), n_ambig
