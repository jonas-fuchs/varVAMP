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


def process_alignment(preprocessed_alignment, threshold):
    """
    - build an array of shape (n_seq, seq_len)
    - for each column, count how many sequences are '-'
    - mark columns to mask if count > cutoff
    - turn those columns into contiguous regions
    """

    # build char array
    seqs = [seq for seq_id, seq in preprocessed_alignment]
    arr = np.array([list(s) for s in seqs], dtype="U1")
    n_seq, len_seq = arr.shape

    # per-column gap counts
    cols_to_mask = (arr == "-").sum(axis=0) > n_seq * (1 - threshold)

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
    if in_gap:
        gaps_to_mask.append([start, len_seq - 1])

    if not gaps_to_mask:
        return preprocessed_alignment, []

    alignment_cleaned = clean_gaps(preprocessed_alignment, gaps_to_mask)

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
