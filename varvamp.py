#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
VARVAMP primer design for viruses with highly variable genomes.
"""

# INFO
__author__ = "Dr. Jonas Fuchs"
__copyright__ = "Copyright 2023"
__license__ = "GPL"
__version__ = "0.1"
__email__ = "jonas.fuchs@uniklinik-freiburg.de"
__status__ = "Development"

# BUILT-INS
import sys
import os
import shutil
import time
import re
import argparse
import collections

# LIBS
import primer3 as p3
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq

# CUSTOM
from scr import config

# DEFINITIONS
# def for logging and raising errors:
def varvamp_progress(progress=0, job="", progress_text="", out=sys.stdout):
    """
    progress bar and logging
    """
    barLength = 40
    block = int(round(barLength*progress))

    if progress == 0:
        print(
            f"\nStarting \033[31m\033[1mvarVAMP ◥(ºwº)◤\033[0m primer design\n",
            file = out,
            flush = True
        )
        if not os.path.exists(results):
            os.makedirs(results)
        else:
            shutil.rmtree(results)
            os.makedirs(results)
        with open(results+"/varvamp_log.txt", 'w') as f:
            f.write('VARVAMP log \n')
    else:
        if progress == 1:
            stop_time = str(round(time.process_time() - start_time, 2))
            progress_text = "all done \n\n\rvarVAMP created an amplicon scheme in " + stop_time + " sec!\n"
            job = "Finalizing output"
        print(
            f"\rJob:\t\t "+job+"\nProgress: \t [{0}] {1}%".format("█"*block + "-"*(barLength-block), progress*100) + "\t" + progress_text,
            file = out,
            flush= True
        )
        with open(results+"/varvamp_log.txt", 'a') as f:
            print(
                "\rJob:\t" + job + "\nResult:\t" + progress_text,
                file = f
            )

def raise_arg_errors(args):
    """
    checks arguments for non-valid input
    """
    # threshold error
    if threshold > 1 or threshold < 0:
        print(
            "\033[31m\033[1mError:\033[0m Threshold can only be between 0-1",
            file=sys.stderr
        )
        sys.exit()

# defs for alignment preprocessing and gap cleaning
def read_alignment(alignment_path):
    """
    read alignment with AlignIO and
    convert to list of lists
    """
    alignment_list = []

    for sequence in AlignIO.read(alignment_path, "fasta"):
        alignment_list.append([sequence.id, str(sequence.seq)])

    return alignment_list

def determine_gap_cutoff(n_seqs):
    """
    determine the cutoff for gaps that
    are covered by at least n seqs
    """
    return int(n_seqs*(1-threshold))

def preprocess_alignment(alignment):
    """
    force nucleotides to lower and
    back transcripe if its RNA
    """
    preprocessed_alignment = []

    for sequence in alignment:
        seq = Seq(sequence[1])
        seq = seq.lower()
        if "u" in seq:
            seq = seq.back_transcribe()
        preprocessed_alignment.append([sequence[0], str(seq)])

    return preprocessed_alignment

def find_gaps_in_alignment(alignment):
    """
    find all gaps for each sequence in alignment
    """
    all_gaps = []

    for seq in alignment:
        # find all gaps for all sequences with regular expression -{min}
        all_gaps.append(
            [(gap.start(0), gap.end(0)-1) for gap in re.finditer(
                "-{"+str(config.DELETION_LENGTH_MIN)+",}", seq[1])]
            )

    return all_gaps

def find_unique_gaps(all_gaps):
    """
    get all unique gaps
    """
    result = list(set(gaps for gap_list in all_gaps for gaps in gap_list))
    return result

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
            if intersection:
                if min(intersection) == unique_gap[0] and max(intersection)+1 == unique_gap[1]:
                    overlapping_gaps.append(unique_gap)

    return overlapping_gaps

def create_gap_dictionary(unique_gaps, all_gaps):
    """
    creates a dictionary with gap counts.
    counts also all overlapping gaps per gap.
    """

    gap_dict = {}

    for gap_list in all_gaps:
        for gap in gap_list:
            overlapping_gaps = find_internal_gaps(unique_gaps, gap)
            for overlapping_gap in overlapping_gaps:
                if overlapping_gap in gap_dict:
                    gap_dict[overlapping_gap] += 1
                else:
                    gap_dict[overlapping_gap] = 1

    return gap_dict

def find_gaps_to_mask(gap_dict, lower_cutoff):
    """
    filters gaps for their freq cutoff.
    condenses final gaps if there is
    an overlap.
    """
    gaps_to_mask = []
    potential_gaps = []

    # check for each region if it is covered
    # by enough sequences
    for gap in gap_dict:
        if gap_dict[gap] > lower_cutoff:
            potential_gaps.append(gap)

    # sort by start and stop
    potential_gaps = sorted(potential_gaps)

    # get the min and max of overlapping gaps
    opened_region = []
    gaps_to_mask = []
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
    mask = config.MASK_LENGTH*"N"

    for sequence in alignment:
        start = 0
        masked_seq = str()
        for region in gaps_to_mask:
            stop = region[0]
            masked_seq_temp = sequence[1][start:stop]
            # check if the deletion is at the start
            if len(masked_seq_temp) != 0:
                masked_seq = (masked_seq + mask + masked_seq_temp)
            start = region[1]+1
        if max(gaps_to_mask)[1] < len(sequence[1])-1:
        # append the last gaps if it is not
        # the end of the sequence
            start = max(gaps_to_mask)[1]
            stop = len(sequence[1])-1
            masked_seq_temp = sequence[1][start:stop]
            masked_seq = (masked_seq + mask + masked_seq_temp)
        else:
        # append the mask to the end of the seq
            masked_seq = masked_seq + mask

        cleaned_alignment.append([sequence[0], masked_seq])

    return cleaned_alignment

def calculate_total_masked_gaps(gaps_to_mask):
    """
    calculates the cummulative length of gaps
    that were masked.
    """
    if gaps_to_mask:
        sum_gaps = 0
        for region in gaps_to_mask:
            sum_gaps += region[1]- region[0] + 1
        return sum_gaps
    else:
        return 0

def preprocess_and_clean_alignment(alignment_path):
    """
    proprocesses alignment and cleans gaps
    """
    alignment = read_alignment(alignment_path)
    alignment_preprocessed = preprocess_alignment(alignment)
    n_seqs = len(alignment)
    gap_cutoff = determine_gap_cutoff(n_seqs)
    all_gaps = find_gaps_in_alignment(alignment_preprocessed)
    unique_gaps = find_unique_gaps(all_gaps)
    if unique_gaps:
        gap_dic = create_gap_dictionary(unique_gaps, all_gaps)
        gaps_to_mask = find_gaps_to_mask(gap_dic, gap_cutoff)
        alignment_cleaned = clean_gaps(alignment_preprocessed, gaps_to_mask)
    else:
        gaps_to_mask = []
        alignment_cleaned = alignment_preprocessed

    return alignment_cleaned, gaps_to_mask, n_seqs

# defs for consensus creation
def determine_consensus_cutoff(n_seqs):
    """
    determine the cutoff to consider a nuc conserved
    """
    return int(n_seqs*threshold)

def determine_nucleotide_counts(alignment, idx):
    """
    count the number of each nucleotides at
    an idx of the alignment. return sorted dic.
    handels ambiguous nucleotides in sequences.
    also handels gaps.
    """
    nucleotide_list = []

    # get all nucleotides
    for sequence in alignment:
        nucleotide_list.append(sequence[1][idx])
    # count occurences of nucleotides
    counter = dict(collections.Counter(nucleotide_list))
    # get permutations of an ambiguous nucleotide
    to_delete = []
    temp_dict = {}
    for nucleotide in counter:
        if nucleotide in config.ambig:
            to_delete.append(nucleotide)
            permutations = config.ambig[nucleotide]
            adjusted_freq = 1/len(permutations)
            for permutation in permutations:
                if permutation in temp_dict:
                    temp_dict[permutation] += adjusted_freq
                else:
                    temp_dict[permutation] = adjusted_freq
        if nucleotide == "-":
            to_delete.append(nucleotide)

    # drop ambiguous entrys and add adjusted freqs to
    if to_delete:
        for i in to_delete:
            counter.pop(i)
        for nucleotide in temp_dict:
            if nucleotide in counter:
                counter[nucleotide] += temp_dict[nucleotide]
            else:
                counter[nucleotide] = temp_dict[nucleotide]

    return dict(sorted(counter.items(), key=lambda x:x[1], reverse = True))

def get_consensus_nucleotides(nucleotide_counts):
    """
    get a list of nucleotides for the consensus seq
    """
    n = 0
    consensus_cutoff = determine_consensus_cutoff(n_seqs)

    consensus_nucleotides = []
    for nuc in nucleotide_counts:
        n += nucleotide_counts[nuc]
        consensus_nucleotides.append(nuc)
        if n >= consensus_cutoff:
            break

    return consensus_nucleotides

def get_ambiguous_char(nucleotides):
    """
    get ambiguous char from a list of nucleotides
    """
    for ambiguous_nuc, permutations in config.ambig.items():
        if set(permutations) == set(nucleotides):
            return ambiguous_nuc

def create_consensus_sequences(alignment):
    """
    build a majority sequence and a sequence that
    has ambiguous chars as determined by the freq
    threshold.
    """

    # ini the consensus seq
    ambiguous_consensus = str()
    majority_consensus = str()

    # define length of the consensus from the first seq in alignment
    length_consensus = len(alignment[0][1])
    # built consensus sequences
    for idx in range(length_consensus):
        nucleotide_counts = determine_nucleotide_counts(alignment, idx)
        consensus_nucleotide = get_consensus_nucleotides(nucleotide_counts)
        if len(consensus_nucleotide) > 1:
            amb_consensus_nucleotide = get_ambiguous_char(consensus_nucleotide)
            ambiguous_consensus = ambiguous_consensus + amb_consensus_nucleotide
        else:
            ambiguous_consensus = ambiguous_consensus + consensus_nucleotide[0]

        majority_consensus = majority_consensus + consensus_nucleotide[0]

    return majority_consensus, ambiguous_consensus

# defs for conserved region search:
def find_conserved_regions(consensus_amb):
    """
    finds conserved regions as specified by a
    certain amount of ambiguous bases in a given
    sequence length
    """
    # init the variables
    current_window = []
    windows_written = 0
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
                if current_window[1] >= config.ALLOWED_N_AMB or nuc == "N":
                # check if there were too many previous amb char in subwindow
                # and make it writable. Always make it writeable if N is
                # reached
                    writable = True
                    if amb_to_amb_len >= config.CONSERVED_MIN_LENGTH and nuc != "N":
                    # check if the last amb is sufficiently far, if yes keep
                    # window open and set amb counter to 0, reset also the
                    # list of amb positions and track only the current pos
                        current_window[1] = 0
                        writable = False
                        amb_pos = [idx]

                current_window[1] += 1

                if writable:
                    writable = False
                    window_length = idx-1-current_window[0]
                    if window_length >= config.CONSERVED_MIN_LENGTH:
                    # check if the writable window has a sufficient length.
                        conserved_regions.append([current_window[0], idx-1,])
                        windows_written += 1
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

def mean_conserved(conserved_regions, consensus):
    """
    calculate the percentage of regions
    that are conserved
    """
    sum = 0
    for region in conserved_regions:
        sum += region[1]-region[0]
    return round(sum/len(consensus)*100,1)

# defs for kmer production
def digest_seq(seq, kmer_size):
    """
    digest the sequence into kmers
    """
    return[[seq[i:i+kmer_size],i, i+len(seq[i:i+kmer_size])] for i in range(len(seq)-kmer_size+1)]

def produce_kmers(conserved_regions):
    """
    produce kmers for all conserved regions
    """
    kmers = []

    for region in conserved_regions:
        sliced_seq = majority_consensus[region[0]:region[1]]
        for kmer_size in range(config.PRIMER_SIZE[0], config.PRIMER_SIZE[1]+1):
            kmers_temp = digest_seq(sliced_seq, kmer_size)
            # adjust the start and stop position of the kmers
            for kmer_temp in kmers_temp:
                kmer_temp[1] = kmer_temp[1]+region[0]
                kmer_temp[2] = kmer_temp[2]+region[0]
            kmers += kmers_temp

    return kmers

# defs for primer calculation and scoring
def calc_gc(primer):
    """
    calculate the gc of a sequence
    """
    return 100*(primer.count("g")+primer.count("G")+primer.count("c")+primer.count("C"))/len(primer)

def calc_temp(primer):
    """
    calculate the melting temperature
    """
    return p3.calcTm(
            primer,
            mv_conc = config.MV_CONC,
            dv_conc = config.DV_CONC,
            dntp_conc = config.DNTP_CONC,
            dna_conc = config.DNA_CONC
        )

def calc_hairpin(primer):
    """
    calculates hairpins
    """
    return p3.calcHairpin(
            primer,
            mv_conc = config.MV_CONC,
            dv_conc = config.DV_CONC,
            dntp_conc = config.DNTP_CONC,
            dna_conc = config.DNA_CONC
        )

def calc_max_polyx(primer):
    """
    calculate maximum polyx in primers
    """
    previous_nuc = primer[0]
    counter = 1
    max_polyx = 1
    for nuc in primer[1:]:
        if nuc == previous_nuc:
            counter += 1
        else:
            counter = 1
            previous_nuc = nuc
        if counter > max_polyx:
            max_polyx = counter
    return(max_polyx)

def calc_max_dinuc_repeats(primer):
    """
    calculate the amount of repeating dinucleotides
    """
    for s in [primer, primer[1:]]:
        previous_dinuc = s[0:2]
        max_dinuc = 0
        counter = 0
        for i in range(2,len(s),2):
            if s[i:i+2] == previous_dinuc:
                counter += 1
            else:
                if counter > max_dinuc:
                    max_dinuc = counter
                counter = 0
                previous_dinuc = s[i:i+2]
    return max_dinuc

def three_prime_ambiguous(amb_primer):
    """
    determine if a sequence contains an ambiguous char at the 3'prime
    """
    len_3_prime = config.MIN_3_WITHOUT_AMB

    if len_3_prime != 0:
        for nuc in amb_primer[len(amb_primer)-len_3_prime:]:
            if nuc not in config.nucs:
                is_amb = True
                break
            else:
                is_amb = False
    else:
        is_amb = False

    return is_amb

def calc_base_penalty(primer):
    """
    Calculate intrinsic primer penalty.
    """
    penalty = 0

    tm = calc_temp(primer)
    gc = calc_gc(primer)
    size = len(primer)

    # TEMP penalty
    if tm > config.PRIMER_TMP[2]:
        penalty += config.PRIMER_TM_PENALTY*(
            tm - config.PRIMER_TMP[2]
            )
    if tm < config.PRIMER_TMP[2]:
        penalty += config.PRIMER_TM_PENALTY*(
            config.PRIMER_TMP[2] - tm
            )
    # GC penalty
    if gc > config.PRIMER_GC_RANGE[2]:
        penalty += config.PRIMER_GC_PENALTY*(
            gc - config.PRIMER_GC_RANGE[2]
            )
    if gc < config.PRIMER_GC_RANGE[2]:
        penalty += config.PRIMER_GC_PENALTY*(
            config.PRIMER_GC_RANGE[2] - gc
        )
    # SIZE penalty
    if size > config.PRIMER_SIZE[2]:
        penalty += config.PRIMER_SIZE_PENALTY*(
            size - config.PRIMER_SIZE[2]
        )
    if size < config.PRIMER_SIZE[2]:
        penalty += config.PRIMER_SIZE_PENALTY * (
            config.PRIMER_SIZE[2] - size
        )

    return penalty

def rev_complement(seq):
    """
    reverse complement a sequence
    """
    return(str(Seq(seq).reverse_complement()))

def hardfilter_primers(primer):
    """
    hard filter primers for temperature, gc content,
    poly x and dinucleotide repeats.
    """
    return(
        (config.PRIMER_TMP[0] <= calc_temp(primer) <= config.PRIMER_TMP[1]) and
        (config.PRIMER_GC_RANGE[0] <= calc_gc(primer) <= config.PRIMER_GC_RANGE[1]) and
        (calc_max_polyx(primer) <= config.MAX_POLYX) and
        (calc_max_dinuc_repeats(primer) <= config.MAX_DINUC_REPEATS)
    )

def filter_primer_direction_specific(direction, primer):
    """
    filter for 3'ambiguous and hairpin - this differs
    depending on the direction of the primer.
    """
    if direction == "LEFT":
        amb_kmer = ambiguous_consensus[primer[1]:primer[2]]
        hairpin_tm = calc_hairpin(primer[0]).tm
    elif direction == "RIGHT":
        amb_kmer = rev_complement(ambiguous_consensus[primer[1]:primer[2]])
        hairpin_tm = calc_hairpin(rev_complement(primer[0])).tm
    if hairpin_tm <= config.PRIMER_HAIRPIN:
        if not three_prime_ambiguous(amb_kmer):
            return primer

def find_lowest_scoring(direction, hardfiltered_kmers):
    """
    sort the primers by start and base penalty.
    for primers that have the same start, retain only
    the lowest scoring and now give a fixed number.
    """
    candidates = []

    # sort for start of the primer and score
    if direction == "LEFT":
        start = 1   # start index of the forward primer
    elif direction == "RIGHT":
        start = 2   # stop index asit is the start of the reverse primer
    hardfiltered_kmers.sort(key = lambda x:(x[start], x[3]))

    # pick the lowest scoring primer for primers that have the same start
    primer_index = 1
    prev_start = -1
    for primer in hardfiltered_kmers:
        if primer[start] != prev_start:
            primer.append(direction+"_"+str(primer_index))
            candidates.append(primer)
            prev_start = primer[start]
            primer_index += 1

    return candidates

def primer_per_base_mismatch(primer, alignment):
    """
    calculate for a given primer with [seq, start, stop]
    percent mismatch per primer pos with the alignment.
    considers if primer or sequences have an amb nuc.
    """
    primer_per_base_mismatch = len(primer[0])*[0]

    for sequence in alignment:
        # slice each sequence of the alignment for the primer
        # start and stop positions
        seq_slice = sequence[1][primer[1]:primer[2]]
        for idx, slice_nuc in enumerate(seq_slice):
            # find the respective nuc to that of the slice
            current_primer_pos = primer[0][idx]
            if slice_nuc != current_primer_pos:
                # check if the slice nucleotide is an amb pos
                if slice_nuc in config.ambig:
                    # check if the primer has an amb pos
                    if current_primer_pos in config.ambig:
                        slice_nuc_set = set(config.ambig[slice_nuc])
                        pri_set = set(config.ambig[current_primer_pos])
                        # check if these sets have no overlap
                        # -> mismatch
                        if len(slice_nuc_set.intersection(pri_set)) == 0:
                            primer_per_base_mismatch[idx] += 1
                    # if no amb pos is in primer then check if primer nuc
                    # is part of the amb slice nuc
                    elif current_primer_pos not in config.ambig[slice_nuc]:
                            primer_per_base_mismatch[idx] += 1
                # check if primer has an amb pos but the current
                # slice_nuc is not part of this amb nucleotide
                elif current_primer_pos in config.ambig:
                        if slice_nuc not in config.ambig[current_primer_pos]:
                            primer_per_base_mismatch[idx] += 1
                # mismatch
                else:
                    primer_per_base_mismatch[idx] += 1

    # gives a percent mismatch over all positions of the primer from 5' to 3'
    primer_per_base_mismatch = [round(x/n_seqs,2) for x in primer_per_base_mismatch]

    return primer_per_base_mismatch

def penalty_3_prime(direction, primer):
    """
    calculate the penalty for mismatches at the 3' end.
    the more mismatches are closer to the 3' end of the primer,
    the higher the penalty.
    """
    penalty = 0

    for i in range(0,len(config.PRIMER_3_PENALTY)):
        if direction == "RIGHT":
            penalty += primer[5][i]*config.PRIMER_3_PENALTY[i]
        elif direction == "LEFT":
            penalty += primer[5][len(primer[0])-1-i]*config.PRIMER_3_PENALTY[i]

    return(penalty)

def find_potenital_primers(kmers):
    """
    hardfilter kmers and further filter
    for potential primers
    """

    hardfiltered_left_kmers = []
    hardfiltered_right_kmers = []

    for kmer in kmers:
        if hardfilter_primers(kmer[0]):
            base_penalty = calc_base_penalty(kmer[0])
            if base_penalty <= config.PRIMER_MAX_BASE_PENALTY:
                if filter_primer_direction_specific("LEFT", kmer):
                    hardfiltered_left_kmers.append(
                        [kmer[0], kmer[1], kmer[2], base_penalty]
                    )
                if filter_primer_direction_specific("RIGHT", kmer):
                    hardfiltered_right_kmers.append(
                        [kmer[0], kmer[1], kmer[2], base_penalty]
                    )

    # filter kmers and complement kmers for possible primers
    left_primer_candidates = find_lowest_scoring("LEFT", hardfiltered_left_kmers)
    right_primer_candidates = find_lowest_scoring("RIGHT", hardfiltered_right_kmers)

    # now calculate the mismatches for each position in the primer.
    # based on this score calculate the 3' penalty and add to base penalty.
    for direction, primer_candidates in [("LEFT", left_primer_candidates), ("RIGHT", right_primer_candidates)]:
        for primer in primer_candidates:
            primer.append(primer_per_base_mismatch(primer, alignment_cleaned))
            primer[3] = primer[3] + penalty_3_prime(direction, primer)

    return left_primer_candidates, right_primer_candidates

if __name__ == "__main__":

    # arg parsing
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "alignment",
        help = "alignment to design primers on"
    )
    parser.add_argument(
        "results",
        help = "path for results dir"
    )
    parser.add_argument(
        "-t",
        "--threshold",
        type = float,
        default = config.FREQUENCY_THRESHOLD,
        help = "threshold for nucleotides in alignment to be considered conserved"
    )

    # define argument variables and verify
    args = parser.parse_args()
    results = args.results
    threshold = args.threshold
    raise_arg_errors(args)

    # ini progress
    start_time = time.process_time()
    varvamp_progress()

    # preprocess and clean alignment of gaps
    alignment_cleaned, gaps_to_mask, n_seqs = preprocess_and_clean_alignment(args.alignment)

    # progress update
    varvamp_progress(
        0.1,
        "Preprocessing alignment and cleaning gaps.",
        str(len(gaps_to_mask)) + " gaps with " +
        str(calculate_total_masked_gaps(gaps_to_mask) ) + " nucleotides"
    )

    # create consensus sequences
    majority_consensus, ambiguous_consensus = create_consensus_sequences(alignment_cleaned)

    # progress update
    varvamp_progress(
        0.2,
        "Creating consensus sequences.",
        "length of the consensus is " + str(len(majority_consensus)) + " nt"
    )

    # generate conserved region list
    conserved_regions = find_conserved_regions(ambiguous_consensus)

    # progress update
    varvamp_progress(
        0.3,
        "Finding conserved regions.",
        str(mean_conserved(conserved_regions, majority_consensus))+"% conserved"
    )

    # produce kmers for all conserved regions
    kmers = produce_kmers(conserved_regions)

    # progress update
    varvamp_progress(
        0.4,
        "Digesting into kmers.",
        str(len(kmers))+" kmers"
    )

    # hard filter kmers for gc, size, temp and base penatly
    # filter direction specific
    left_primer_candidates, right_primer_candidates = find_potenital_primers(kmers)

    # progress update
    varvamp_progress(
        0.5,
        "Filtering for primers.",
        str(len(left_primer_candidates)) +
        " fw and " + str(len(right_primer_candidates)) +
        " rw potential primers"
    )

    # final progress
    varvamp_progress(1)





########## TODO LIST ###########

##### Primer design #####
# TODO: Create graph
# TODO: Score Degeneracy
# TODO: Overthink the final max score
# TODO: Heuristic
# TODO: Find best path

##### Output #####
# TODO: write BED files for conserved regions
# TODO: write consensus sequences
# TODO: write all primers as BED
# TODO: write final primers and amplicons (formatted for bioinformatic workflow)

##### Reporting #####
# TODO: Write used settings
# TODO: Create html
# TODO: Create plot for alignment entropy and primer location
# TODO: Create visualization of primer mismatches

##### Finalize #####
# TODO: New Github repo with current dev branch
# TODO: Error logging
# TODO: Confirm config
# TODO: readme and further documentation
