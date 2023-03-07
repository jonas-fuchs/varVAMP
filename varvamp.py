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
__version__ = "0.1"
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


# DEFs
def varvamp_progress(progress=0, job="", progress_text="", out=sys.stdout):
    """
    progress bar and logging
    """
    barLength = 40
    block = int(round(barLength*progress))

    if progress == 0:
        print(
            "\nStarting \033[31m\033[1mvarVAMP ◥(ºwº)◤\033[0m primer design\n",
            file=out,
            flush=True
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
            "\rJob:\t\t " + job + "\nProgress: \t [{0}] {1}%".format("█"*block + "-"*(barLength-block), progress*100) + "\t" + progress_text,
            file=out,
            flush=True
        )
        with open(results+"/varvamp_log.txt", 'a') as f:
            print(
                "\rJob:\t" + job + "\nResult:\t" + progress_text,
                file=f
            )


def raise_arg_errors(args):
    """
    checks arguments for non-valid input
    """
    # threshold error
    if args.threshold > 1 or args.threshold < 0:
        sys.exit("\n\033[31m\033[1mERROR:\033[0m Threshold can only be between 0-1\n")
    if args.allowed_ambiguous < 0:
        sys.exit("\n\033[31m\033[1mERROR:\033[0m Set the number of ambiguous nucleotides >= 0.\n")
    if args.allowed_ambiguous > 4:
        print("\n\033[31m\033[1mWARNING:\033[0m High number of ambiguous nucleotides in primer leads to a high degeneracy. Condider reducing.")


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

    # define argument variables and verify
    args = parser.parse_args()
    raise_arg_errors(args)

    # ini progress
    results = args.results
    start_time = time.process_time()
    varvamp_progress()

    # preprocess and clean alignment of gaps
    alignment_cleaned, gaps_to_mask = alignment.process_alignment(
        args.alignment,
        args.threshold
    )

    # progress update
    varvamp_progress(
        0.1,
        "Preprocessing alignment and cleaning gaps.",
        str(len(gaps_to_mask)) + " gaps with "
        + str(alignment.calculate_total_masked_gaps(gaps_to_mask)) + " nucleotides"
    )

    # create consensus sequences
    majority_consensus, ambiguous_consensus = consensus.create_consensus(
        alignment_cleaned,
        args.threshold
    )

    # progress update
    varvamp_progress(
        0.2,
        "Creating consensus sequences.",
        "length of the consensus is " + str(len(majority_consensus)) + " nt"
    )

    # generate conserved region list
    conserved_regions = conserved.find_regions(
        ambiguous_consensus,
        args.allowed_ambiguous
    )

    # progress update
    varvamp_progress(
        0.3,
        "Finding conserved regions.",
        str(conserved.mean(conserved_regions, majority_consensus))+"% conserved"
    )

    # produce kmers for all conserved regions
    kmers = conserved.produce_kmers(
        conserved_regions,
        majority_consensus
    )

    # progress update
    varvamp_progress(
        0.4,
        "Digesting into kmers.",
        str(len(kmers))+" kmers"
    )

    # find potential primers
    left_primer_candidates, right_primer_candidates = primers.find_primers(
        kmers,
        ambiguous_consensus,
        alignment_cleaned
    )

    # progress update
    varvamp_progress(
        0.5,
        "Filtering for primers.",
        str(len(left_primer_candidates))
        + " fw and " + str(len(right_primer_candidates))
        + " rw potential primers"
    )

    # final progress
    varvamp_progress(1)
