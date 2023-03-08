#!/usr/bin/env python
"""
varVAMP primer design for viruses with highly variable genomes. varVAMP
first preprocesses the alignment and then creates consensus sequences
that can contain ambiguous characters. Then it searches for conserved
regions as defined by a user defined amount of ambiguous charaters within
the min length of a primer. The conserved regions of a consensus sequence
containing the most prevalent nucleotide (no wobbels) is then digested into
kmers that which are tested for potential primers. Each primer is given a
penalty score.
"""

# INFO
__author__ = "Dr. Jonas Fuchs"
__copyright__ = "Copyright 2023"
__license__ = "GPL"
__version__ = "0.2"
__email__ = "jonas.fuchs@uniklinik-freiburg.de"
__status__ = "Development"

# BUILT-INS
import sys
import os
import shutil
import time
import argparse

# varVAMP
from scr import config
from scr import alignment
from scr import consensus
from scr import conserved
from scr import primers
from scr import scheme


# DEFs
def varvamp_progress(progress=0, job="", progress_text="", out=sys.stdout):
    """
    progress bar, logging and folder creation
    """
    barLength = 40
    block = int(round(barLength*progress))

    if progress == 0:
        if args.console:
            print(
                "\nStarting \033[31m\033[1mvarVAMP ◥(ºwº)◤\033[0m primer design\n",
                file=out,
                flush=True
            )
        if not os.path.exists(results_dir):
            os.makedirs(results_dir)
        else:
            shutil.rmtree(results_dir)
            os.makedirs(results_dir)
        os.makedirs(all_data_dir)
        with open(results_dir+"varvamp_log.txt", 'w') as f:
            f.write('VARVAMP log \n')
    else:
        if progress == 1:
            stop_time = str(round(time.process_time() - start_time, 2))
            progress_text = f"all done \n\n\rvarVAMP created an amplicon scheme in {stop_time} sec!\n"
            job = "Finalizing output."
        if args.console:
            print(
                "\rJob:\t\t " + job + "\nProgress: \t [{0}] {1}%".format("█"*block + "-"*(barLength-block), progress*100) + "\t" + progress_text,
                file=out,
                flush=True
            )
        with open(results_dir+"varvamp_log.txt", 'a') as f:
            print(
                f"\rJob:\t {job} \nResult:\t {progress_text}",
                file=f
            )


def raise_arg_errors(args):
    """
    checks arguments for non-valid input
    """
    # threshold error
    if args.threshold > 1 or args.threshold < 0:
        sys.exit("\n\033[31m\033[1mERROR:\033[0m threshold can only be between 0-1\n")
    if args.allowed_ambiguous < 0:
        sys.exit("\n\033[31m\033[1mERROR:\033[0m set the number of ambiguous nucleotides >= 0.\n")
    if args.allowed_ambiguous > 4:
        print("\n\033[31m\033[1mWARNING:\033[0m high number of ambiguous nucleotides in primer leads to a high degeneracy. Condider reducing.")
    if args.opt_length > args.max_length:
        sys.exit("\n\033[31m\033[1mERROR:\033[0m Optimal length can not be higher than the maximum amplicon length.\n")
    if args.opt_length < 0 or args.max_length < 0:
        sys.exit("\n\033[31m\033[1mERROR:\033[0m amplicon lengths can not be negative.\n")
    if args.opt_length < 200 or args.max_length < 200:
        print("\n\033[31m\033[1mWARNING:\033[0m your amplicon lengths might be to small. Consider increasing")
    if args.overlap < 0:
        sys.exit("\n\033[31m\033[1mERROR:\033[0m overlap size can not be negative.\n")
    if args.overlap < 50:
        print("\n\033[31m\033[1mWARNING:\033[0m small overlaps might hinder downstream analyses. Consider increasing.")
    if args.overlap > args.opt_length:
        sys.exit("\n\033[31m\033[1mERROR:\033[0m overlaps can not be higher than the length of amplicons.\n")


if __name__ == "__main__":

    # arg parsing
    parser = argparse.ArgumentParser()
    if len(sys.argv[1:]) < 1:
        parser.print_help()
        sys.exit("\033[31m\033[1mError:\033[0m No arguments")

    parser.add_argument(
        "alignment",
        help="alignment to design primers on"
    )
    parser.add_argument(
        "results",
        help="path for results dir"
    )
    parser.add_argument(
        "-ol",
        "--opt-length",
        help="optimal length of the amplicons",
        type=int,
        default=config.OPT_AMPLICON_LENGTH
    )
    parser.add_argument(
        "-ml",
        "--max-length",
        help="max length of the amplicons",
        type=int,
        default=config.MAX_AMPLICON_LENGTH
    )
    parser.add_argument(
        "-o",
        "--overlap",
        type=float,
        default=config.MIN_OVERLAP,
        help="min overlap of the amplicons"
    )
    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        default=config.FREQUENCY_THRESHOLD,
        help="threshold for nucleotides in alignment to be considered conserved"
    )
    parser.add_argument(
        "-a",
        "--allowed-ambiguous",
        type=int,
        default=config.ALLOWED_N_AMB,
        help="number of ambiguous characters that are allowed within a primer"
    )

    parser.add_argument(
        "--console",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="show varvamp console output"
    )
    # define argument variables and verify
    args = parser.parse_args()
    raise_arg_errors(args)

    # ini progress
    start_time = time.process_time()
    results_dir = args.results+"/"
    all_data_dir = args.results+"all_data/"
    varvamp_progress()

    # check if config is ok
    config.confirm_config()

    varvamp_progress(
        0.1,
        "Checking config.",
        "config file passed"
    )

    # preprocess and clean alignment of gaps
    alignment_cleaned, gaps_to_mask = alignment.process_alignment(
        args.alignment,
        args.threshold
    )

    # progress update
    varvamp_progress(
        0.2,
        "Preprocessing alignment and cleaning gaps.",
        f"{len(gaps_to_mask)} gaps with {alignment.calculate_total_masked_gaps(gaps_to_mask)} nucleotides"
    )

    # create consensus sequences
    majority_consensus, ambiguous_consensus = consensus.create_consensus(
        alignment_cleaned,
        args.threshold
    )

    # progress update
    varvamp_progress(
        0.3,
        "Creating consensus sequences.",
        f"length of the consensus is {len(majority_consensus)} nt"
    )

    # generate conserved region list
    conserved_regions = conserved.find_regions(
        ambiguous_consensus,
        args.allowed_ambiguous
    )

    # raise error if no conserved regions were found
    if not conserved_regions:
        sys.exit("\n\033[31m\033[1mERROR:\033[0m nothing conserved. Lower the threshod!\n")

    # progress update
    varvamp_progress(
        0.4,
        "Finding conserved regions.",
        f"{conserved.mean(conserved_regions, majority_consensus)} % conserved"
    )

    # produce kmers for all conserved regions
    kmers = conserved.produce_kmers(
        conserved_regions,
        majority_consensus
    )

    # progress update
    varvamp_progress(
        0.5,
        "Digesting into kmers.",
        f"{len(kmers)} kmers"
    )

    # find potential primers
    left_primer_candidates, right_primer_candidates = primers.find_primers(
        kmers,
        ambiguous_consensus,
        alignment_cleaned
    )

    # raise error if no primers were found
    for type, primer_candidates in [("LEFT", left_primer_candidates),("RIGHT", right_primer_candidates)]:
        if not primer_candidates:
            sys.exit(f"\n\033[31m\033[1mERROR:\033[0m no {type} primers found.\n")

    # progress update
    varvamp_progress(
        0.6,
        "Filtering for primers.",
        f"{len(left_primer_candidates)} fw and {len(right_primer_candidates)} rw potential primers"
    )

    # find best primers
    left_primer_candidates, right_primer_candidates = primers.find_best_primers(
        left_primer_candidates,
        right_primer_candidates
    )

    # progress update
    varvamp_progress(
        0.7,
        "Considering only high scoring primers.",
        f"{len(left_primer_candidates)} fw and {len(right_primer_candidates)} rw primers"
    )

    # find all possible amplicons
    amplicons = scheme.find_amplicons(
        left_primer_candidates,
        right_primer_candidates,
        args.opt_length,
        args.max_length
    )

    # raise error if no amplicons were found
    if not amplicons:
        sys.exit("\n\033[31m\033[1mERROR:\033[0m no amplicons found. Increase the max amplicon length or lower threshold!\n")

    # build the amplicon graph
    amplicon_graph = scheme.create_amplicon_graph(amplicons, args.overlap)

    # progress update
    varvamp_progress(
        0.8,
        "Finding potential amplicons.",
        str(len(amplicons)) + " potential amplicons"
    )

    coverage, amplicon_scheme = scheme.find_best_covering_scheme(
        amplicons,
        amplicon_graph
    )

    percent_coverage = round(coverage/len(ambiguous_consensus)*100, 2)

    varvamp_progress(
        0.9,
        "Creating amplicon scheme.",
        f"{percent_coverage} % total coverage with {len(amplicon_scheme)} amplicons"
    )

    # raise low coverage warning
    if percent_coverage < 70:
        print("\n\033[31m\033[1mWARNING:\033[0m coverage < 70 %. Possible solutions:")
        print("\t - lower threshold")
        print("\t - increase amplicons lengths")
        print("\t - increase number of ambiguous nucleotides")
        print("\t - relax primer settings (not recommended) \n")

    # final progress
    varvamp_progress(1)
