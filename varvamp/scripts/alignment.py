"""
alignment preprocessing
"""

# varVAMP
from varvamp.scripts import config

# LIBS
import numpy as np
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
    back transcribe if its RNA
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


def clean_gaps(alignment, gaps_to_mask):
    """
    Clean an alignment of large common deletions based on gaps_to_mask.
    gaps_to_mask: list of [start, end] (inclusive), sorted by start.
    """
    cleaned_alignment = []
    gaps_to_mask = sorted(gaps_to_mask, key=lambda x: x[0])

    for seq_id, seq in alignment:
        start = 0
        pieces = []
        # for each seq in the alignment, mask the regions
        for region_start, region_end in gaps_to_mask:
            # mask length for this region
            if (region_end - region_start + 1) >= config.QAMPLICON_DEL_CUTOFF:
                mask = "NN"
            else:
                mask = "N"
            # part before region
            pieces.append(seq[start:region_start])
            # mask for region
            pieces.append(mask)
            # next start is after region
            start = region_end + 1

        # tail after last masked region
        if start < len(seq):
            pieces.append(seq[start:])

        cleaned_alignment.append([seq_id, "".join(pieces)])

    return cleaned_alignment


def process_alignment(preprocessed_alignment, threshold, terminal_threshold=config.TERMINAL_MASKING_THRESHOLD):
    """
    - build an array of shape (n_seq, seq_len)
    - for each column, count how many sequences are '-'
    - mark columns to mask if count > cutoff
    - handle terminal gaps differently: count > terminal_cutoff
    - turn those columns into contiguous regions
    """

    # build char array
    seqs = [seq for seq_id, seq in preprocessed_alignment]
    arr = np.array([list(s) for s in seqs], dtype="U1")
    n_seq, len_seq = arr.shape

    # determine n terminal gaps
    for i, row in enumerate(arr):
        # mask forward terminal gaps with ~
        for j in range(0, len_seq):
            if arr[i][j] == "-":
                arr[i][j] = '~'
            else:
                break
        # mask reverse terminal gaps with ~
        for j in range(len_seq - 1, -1, -1):
            if arr[i][j] == "-":
                arr[i][j] = '~'
            else:
                break
    # define the gaps that should be masked, this is defined as the number of sequences with gaps
    # that are above the 1-threshold * number of total sequences that do not have terminal gaps.
    # only consider columns where we have enough non-terminal sequences to make a decision
    # exclude columns that will be handled by terminal masking later on
    n_terminal_per_col = (arr == "~").sum(axis=0)
    n_non_terminal_per_col = n_seq - n_terminal_per_col
    is_terminal = n_terminal_per_col > n_seq * (1 - terminal_threshold)
    gaps_above_threshold = (arr == "-").sum(axis=0) > n_non_terminal_per_col * (1 - threshold)
    cols_to_mask = gaps_above_threshold & (n_non_terminal_per_col > 0) & ~is_terminal

    # convert bool mask into list of (start, end) regions (end inclusive)
    gaps_to_mask = []
    in_gap = False
    start = None
    for i, is_gap in enumerate(cols_to_mask):
        if is_gap and not in_gap:
            in_gap = True
            start = i
        elif not is_gap and in_gap:
            in_gap = False
            gaps_to_mask.append([start, i - 1])

    # define which side of the terminal gaps have to be masked
    if any(is_terminal):
        non_mask = np.where(np.diff(is_terminal))[0]
        # case one - we have terminal gaps on both sides in high enough frequency - mask both sides
        if len(non_mask) == 2:
            gaps_to_mask = [[0, non_mask[0]]] + gaps_to_mask + [[non_mask[1] + 1, len_seq - 1]]
        elif len(non_mask) == 1:
            # determine which end to mask based on the state at position 0
            if is_terminal[0]:
                # starts with terminal gaps -> mask start
                gaps_to_mask = [[0, non_mask[0]]] + gaps_to_mask
            else:
                # starts without terminal gaps -> mask end
                gaps_to_mask = gaps_to_mask + [[non_mask[0] + 1, len_seq - 1]]

    # return alignment if no regions need to be masked
    if not gaps_to_mask:
        return preprocessed_alignment, []

    return clean_gaps(preprocessed_alignment, gaps_to_mask), gaps_to_mask


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
