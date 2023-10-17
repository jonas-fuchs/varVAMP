"""
alignment preprocessing
"""

# BUILT-INS
import re
import multiprocessing

# varVAMP
from varvamp.scripts import config

# LIBS
from Bio import AlignIO
from Bio.Seq import Seq


def read_alignment(alignment_path):
    """
    read alignment with AlignIO and
    convert to list of lists
    """
    alignment_list = []

    for sequence in AlignIO.read(alignment_path, "fasta"):
        alignment_list.append([sequence.id, str(sequence.seq)])

    return alignment_list


def preprocess(alignment_path):
    """
    force nucleotides to lower and
    back transcripe if its RNA
    """

    preprocessed_alignment = []
    # read alignment
    alignment = read_alignment(alignment_path)

    for sequence in alignment:
        seq = Seq(sequence[1])
        seq = seq.lower()
        if "u" in seq:
            seq = seq.back_transcribe()
        preprocessed_alignment.append([sequence[0], str(seq)])

    return preprocessed_alignment


def find_internal_gaps(unique_gaps, gap):
    """
    find all unique gaps that
    lie within the current gap
    """
    overlapping_gaps = []

    if gap[1] - gap[0] == 0:
        # if the gap length = 1 there are
        # no overlapping gaps
        overlapping_gaps = [gap]
    else:
        # for each unique gap check if the intersection with the
        # gap is the same as the unique gap -> internal gap of
        # the current gap
        for unique_gap in unique_gaps:
            unique_set = set(range(unique_gap[0], unique_gap[1]))
            current_range = range(gap[0], gap[1])
            intersection = unique_set.intersection(current_range)
            if not intersection:
                continue
            if min(intersection) == unique_gap[0] and max(intersection)+1 == unique_gap[1]:
                overlapping_gaps.append(unique_gap)

    return overlapping_gaps


def find_overlapping_gaps_worker(gap_list, unique_gaps):
    """
    Worker function to find overlapping gaps and count their occurrences.
    """
    gap_dict_part = {}
    for gap in gap_list:
        overlapping_gaps = find_internal_gaps(unique_gaps, gap)
        for overlapping_gap in overlapping_gaps:
            if overlapping_gap in gap_dict_part:
                gap_dict_part[overlapping_gap] += 1
            else:
                gap_dict_part[overlapping_gap] = 1
    return gap_dict_part


def create_gap_dictionary(unique_gaps, all_gaps, n_threads):
    """
    Creates a dictionary with all gap counts.
    Counts also all overlapping gaps per gap.
    Uses multiprocessing for parallelization.
    """

    with multiprocessing.Pool(processes=n_threads) as pool:
        results = pool.starmap(find_overlapping_gaps_worker, [(gap_list, unique_gaps) for gap_list in all_gaps])

    gap_dict = {}
    for gap_dict_part in results:
        for gap, count in gap_dict_part.items():
            if gap in gap_dict:
                gap_dict[gap] += count
            else:
                gap_dict[gap] = count

    return gap_dict



def find_gaps_to_mask(gap_dict, cutoff):
    """
    filters gaps for their freq cutoff.
    condenses final gaps if there is
    an overlap.
    """
    gaps_to_mask = []
    potential_gaps = []
    opened_region = []

    # check for each region if it is covered
    # by enough sequences
    for gap in gap_dict:
        if gap_dict[gap] > cutoff:
            potential_gaps.append(gap)

    # sort by start and stop
    potential_gaps = sorted(potential_gaps)

    # get the min and max of overlapping gaps
    for i, region in enumerate(potential_gaps):
        region = list(region)
        if opened_region:
            # write the opened region if the start of the current region
            # > opened_region[stop] and the last still opened region
            if region[0] > opened_region[1] or i == len(potential_gaps)-1:
                gaps_to_mask.append(opened_region)
                opened_region = region
            else:
                # 1 case: same start and further stop -> new stop
                if region[0] == opened_region[0]:
                    opened_region[1] = region[1]
                # 2 case: further start and further stop -> new stop
                if region[0] > opened_region[0] and region[1] > opened_region[1]:
                    opened_region[1] = region[1]
        else:
            opened_region = region

    return gaps_to_mask


def clean_gaps(alignment, gaps_to_mask):
    """
    clean an alignment of large common deletions.
    """
    cleaned_alignment = []

    for sequence in alignment:
        start = 0
        masked_seq = str()
        for region in gaps_to_mask:
            # check if it is three bases or more and mask with 2 Ns
            if region[1]-region[0] >= config.QAMPLICON_DEL_CUTOFF:
                mask = "NN"
            # or mask with one N (small deletion)
            else:
                mask = "N"
            stop = region[0]
            masked_seq_temp = sequence[1][start:stop]
            # check if the deletion is at the start
            if start == 0 and len(masked_seq_temp) == 0:
                masked_seq = mask
            # check if deletion is not at start
            elif start == 0 and len(masked_seq_temp) != 0:
                masked_seq = masked_seq_temp
            # else we are in the middle of the alignment
            else:
                masked_seq = masked_seq + mask + masked_seq_temp
            start = region[1]+1
        if max(gaps_to_mask)[1] < len(sequence[1])-1:
            # append the last seq if no gap is at
            # the end of the sequence
            start = max(gaps_to_mask)[1]
            stop = len(sequence[1])-1
            masked_seq_temp = sequence[1][start:stop]
            masked_seq = masked_seq + mask + masked_seq_temp
        else:
            # append the mask to the end of the seq
            masked_seq = masked_seq + mask

        cleaned_alignment.append([sequence[0], masked_seq])

    return cleaned_alignment


def process_alignment(preprocessed_alignment, threshold, n_threads):
    """
    proprocesses alignment and cleans gaps
    """
    all_gaps = []

    gap_cutoff = len(preprocessed_alignment)*(1-threshold)
    for seq in preprocessed_alignment:
        # find all gaps for all sequences with regular expression -{min}
        all_gaps.append(
            [(gap.start(0), gap.end(0) - 1) for gap in re.finditer(
                "-{1,}", seq[1])]
        )
    unique_gaps = list(set(gaps for gap_list in all_gaps for gaps in gap_list))

    if unique_gaps:
        gap_dic = create_gap_dictionary(unique_gaps, all_gaps, n_threads)
        gaps_to_mask = find_gaps_to_mask(gap_dic, gap_cutoff)
        if gaps_to_mask:
            alignment_cleaned = clean_gaps(
                preprocessed_alignment, gaps_to_mask
            )
        else:
            alignment_cleaned = preprocessed_alignment
    else:
        gaps_to_mask = []
        alignment_cleaned = preprocessed_alignment

    return alignment_cleaned, gaps_to_mask


def calculate_total_masked_gaps(gaps_to_mask):
    """
    calculates the cummulative length of gaps
    that were masked.
    """
    if gaps_to_mask:
        sum_gaps = 0
        for region in gaps_to_mask:
            sum_gaps += region[1] - region[0] + 1
        return sum_gaps
    else:
        return 0
