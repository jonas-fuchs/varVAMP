"""
main workflow
"""

# BUILT-INS
import sys
import os
import time
import argparse

# varVAMP
from . import _program
from varvamp import __version__
from varvamp.scripts import logging
from varvamp.scripts import alignment
from varvamp.scripts import config
from varvamp.scripts import consensus
from varvamp.scripts import conserved
from varvamp.scripts import primers
from varvamp.scripts import reporting
from varvamp.scripts import scheme


# DEFs
def get_args(sysargs):
    """
    arg parsing for varvamp
    """
    parser = argparse.ArgumentParser(
        prog=_program,
        description='varvamp: variable virus amplicon design',
        usage='''varvamp <alignment> <output dir> [options]''')

    parser.add_argument(
        "input",
        nargs=2,
        help="alignment file and dir to write results"
    )
    parser.add_argument(
        "-ol",
        "--opt-length",
        help="optimal length of the amplicons",
        type=int,
        default=config.AMPLICON_OPT_LENGTH
    )
    parser.add_argument(
        "-ml",
        "--max-length",
        help="max length of the amplicons",
        type=int,
        default=config.AMPLICON_MAX_LENGTH
    )
    parser.add_argument(
        "-o",
        "--overlap",
        type=float,
        default=config.AMPLICON_MIN_OVERLAP,
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
        default=config.PRIMER_ALLOWED_N_AMB,
        help="number of ambiguous characters that are allowed within a primer"
    )
    parser.add_argument(
        "--console",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="show varvamp console output"
    )
    parser.add_argument(
        "-v",
        "--version",
        action='version',
        version=f"varvamp {__version__}"
    )

    if len(sysargs) < 1:
        parser.print_help()
        sys.exit(-1)
    else:
        return parser.parse_args(sysargs)


def main(sysargs=sys.argv[1:]):
    """
    main varvamp workflow
    """
    # start varVAMP
    args = get_args(sysargs)
    if not args.console:
        sys.stdout = open(os.devnull, 'w')
    start_time = time.process_time()
    results_dir, data_dir, log_file = logging.create_dir_structure(args.input[1])
    logging.raise_arg_errors(args, log_file)
    logging.varvamp_progress(log_file)
    # config check
    logging.confirm_config(args, log_file)
    logging.varvamp_progress(
        log_file,
        progress=0.1,
        job="Checking config.",
        progress_text="config file passed"
    )
    # preprocess and clean alignment of gaps
    alignment_cleaned, gaps_to_mask = alignment.process_alignment(
        args.input[0],
        args.threshold
    )
    logging.varvamp_progress(
        log_file,
        progress=0.2,
        job="Preprocessing alignment and cleaning gaps.",
        progress_text=f"{len(gaps_to_mask)} gaps with {alignment.calculate_total_masked_gaps(gaps_to_mask)} nucleotides"
    )
    # create consensus sequences
    majority_consensus, ambiguous_consensus = consensus.create_consensus(
        alignment_cleaned,
        args.threshold
    )
    logging.varvamp_progress(
        log_file,
        progress=0.3,
        job="Creating consensus sequences.",
        progress_text=f"length of the consensus is {len(majority_consensus)} nt"
    )
    # generate conserved region list
    conserved_regions = conserved.find_regions(
        ambiguous_consensus,
        args.allowed_ambiguous
    )
    if not conserved_regions:
        logging.raise_error(
            "nothing conserved. Lower the threshold!",
            log_file,
            exit=True
        )
    logging.varvamp_progress(
        log_file,
        progress=0.4,
        job="Finding conserved regions.",
        progress_text=f"{conserved.mean(conserved_regions, majority_consensus)} % conserved"
    )
    # produce kmers for all conserved regions
    kmers = conserved.produce_kmers(
        conserved_regions,
        majority_consensus
    )
    logging.varvamp_progress(
        log_file,
        progress=0.5,
        job="Digesting into kmers.",
        progress_text=f"{len(kmers)} kmers"
    )
    # find potential primers
    left_primer_candidates, right_primer_candidates = primers.find_primers(
        kmers,
        ambiguous_consensus,
        alignment_cleaned
    )
    for type, primer_candidates in [("+", left_primer_candidates), ("-", right_primer_candidates)]:
        if not primer_candidates:
            logging.raise_error(
                f"no {type} primers found.\n",
                log_file,
                exit=True
            )
    logging.varvamp_progress(
        log_file,
        progress=0.6,
        job="Filtering for primers.",
        progress_text=f"{len(left_primer_candidates)} fw and {len(right_primer_candidates)} rw potential primers"
    )
    # find best primers and create primer dict
    all_primers = primers.find_best_primers(left_primer_candidates, right_primer_candidates)
    logging.varvamp_progress(
        log_file,
        progress=0.7,
        job="Considering only high scoring primers.",
        progress_text=f"{len(all_primers['+'])} fw and {len(all_primers['-'])} rw primers"
    )
    # find all possible amplicons
    amplicons = scheme.find_amplicons(
        all_primers,
        args.opt_length,
        args.max_length
    )
    if not amplicons:
        logging.raise_error(
            "no amplicons found. Increase the max "
            "amplicon length or lower threshold!\n",
            log_file,
            exit=True
        )
    amplicon_graph = scheme.create_amplicon_graph(amplicons, args.overlap)
    logging.varvamp_progress(
        log_file,
        progress=0.8,
        job="Finding potential amplicons.",
        progress_text=str(len(amplicons)) + " potential amplicons"
    )
    # search for amplicon scheme
    coverage, amplicon_scheme = scheme.find_best_covering_scheme(
        amplicons,
        amplicon_graph,
        all_primers
    )
    dimers_not_solved = scheme.check_and_solve_heterodimers(
        amplicon_scheme,
        left_primer_candidates,
        right_primer_candidates,
        all_primers)
    if dimers_not_solved:
        logging.raise_error(
            f"varVAMP found {len(dimers_not_solved)} primer dimers without replacements. Check the dimer file and perform the PCR for incomaptible amplicons in a sperate reaction.",
            log_file
        )
        reporting.write_dimers(dir, dimers_not_solved)
    percent_coverage = round(coverage/len(ambiguous_consensus)*100, 2)
    logging.varvamp_progress(
        log_file,
        progress=0.9,
        job="Creating amplicon scheme.",
        progress_text=f"{percent_coverage} % total coverage with {len(amplicon_scheme[0]) + len(amplicon_scheme[1])} amplicons"
    )
    if percent_coverage < 70:
        logging.raise_error(
            "coverage < 70 %. Possible solutions:\n"
            "\t - lower threshold\n"
            "\t - increase amplicons lengths\n"
            "\t - increase number of ambiguous nucleotides\n"
            "\t - relax primer settings (not recommended)\n",
            log_file
        )
    # write files
    reporting.write_alignment(data_dir, alignment_cleaned)
    reporting.write_fasta(data_dir, "majority_consensus", majority_consensus)
    reporting.write_fasta(results_dir, "ambiguous_consensus", ambiguous_consensus)
    reporting.write_conserved_to_bed(conserved_regions, data_dir)
    reporting.write_all_primers(data_dir, all_primers)
    reporting.write_scheme_to_files(
        results_dir,
        amplicon_scheme,
        ambiguous_consensus
    )
    reporting.varvamp_plot(
        results_dir,
        args.threshold,
        alignment_cleaned,
        conserved_regions,
        all_primers,
        amplicon_scheme,
    )
    logging.varvamp_progress(log_file, progress=1, start_time=start_time)
