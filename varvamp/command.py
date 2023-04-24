"""
main workflow
"""

# BUILT-INS
import sys
import os
import time
import argparse

# varVAMP
from varvamp.scripts import alignment
from varvamp.scripts import config
from varvamp.scripts import consensus
from varvamp.scripts import regions
from varvamp.scripts import logging
from varvamp.scripts import param_estimation
from varvamp.scripts import primers
from varvamp.scripts import qpcr
from varvamp.scripts import reporting
from varvamp.scripts import scheme
from varvamp import __version__
from . import _program


def get_args(sysargs):
    """
    arg parsing for varvamp
    """
    parser = argparse.ArgumentParser(
        prog=_program,
        usage='''\tvarvamp <mode> --help\n\tvarvamp <mode> [mode optional arguments] <alignment> <output dir>''')
    mode_parser = parser.add_subparsers(
        title="varvamp mode",
        dest="mode",
    )
    SANGER_parser = mode_parser.add_parser(
        "sanger",
        help="design primers for sanger sequencing",
        usage="varvamp sanger [optional arguments] <alignment> <output dir>"
    )
    TILED_parser = mode_parser.add_parser(
        "tiled",
        help="design primers for whole genome sequencing",
        usage="varvamp tiled [optional arguments] <alignment> <output dir>"
    )
    QPCR_parser = mode_parser.add_parser(
        "qpcr",
        help="design qPCR primers",
        usage="varvamp qpcr [optional arguments] <alignment> <output dir>"
    )
    parser.add_argument(
        "input",
        nargs=2,
        help="alignment file and dir to write results"
    )
    for par in (SANGER_parser, TILED_parser, QPCR_parser):
        par.add_argument(
            "-t",
            "--threshold",
            metavar="",
            type=float,
            default=None,
            help="threshold for consensus nucleotides"
        )
        par.add_argument(
            "-a",
            "--n-ambig",
            metavar="",
            type=int,
            default=None,
            help="max number of ambiguous characters in a primer"
        )
    for par in (SANGER_parser, TILED_parser):
        par.add_argument(
            "-ol",
            "--opt-length",
            help="optimal length of the amplicons",
            metavar="1000",
            type=int,
            default=1000
        )
        par.add_argument(
            "-ml",
            "--max-length",
            help="max length of the amplicons",
            metavar="1500",
            type=int,
            default=1500
        )
    TILED_parser.add_argument(
        "-o",
        "--overlap",
        type=int,
        metavar="100",
        default=100,
        help="min overlap of the amplicons"
    )
    SANGER_parser.add_argument(
        "-n",
        "--report-n",
        type=int,
        metavar="inf",
        default=float("inf"),
        help="report the top n best hits"
    )
    QPCR_parser.add_argument(
        "-pa",
        "--pn-ambig",
        metavar="1",
        type=int,
        default=None,
        help="max number of ambiguous characters in a probe"
    )
    QPCR_parser.add_argument(
        "-n",
        "--test-n",
        type=int,
        metavar="50",
        default=50,
        help="test the top n qPCR amplicons for secondary structures at the minimal primer temperature"
    )
    QPCR_parser.add_argument(
        "-d",
        "--deltaG",
        type=int,
        metavar="-1",
        default=-1,
        help="minimum free energy (kcal/mol/K) cutoff at the lowest primer melting temp"
    )
    parser.add_argument(
        "--verbose",
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


def shared_workflow(args, log_file):
    """
    part of the workflow that is shared by all modes
    """
    # start varvamp
    logging.varvamp_progress(log_file, mode = args.mode)

    # read in alignment and preprocess
    preprocessed_alignment = alignment.preprocess(args.input[0])

    # estimate threshold or number of ambiguous bases if args were not supplied
    if args.threshold is None or args.n_ambig is None:
        args.threshold, args.n_ambig = param_estimation.get_parameters(preprocessed_alignment, args, log_file)
    if args.mode == "qpcr" and args.n_ambig >= 1 and args.pn_ambig == None:
        args.pn_ambig = args.n_ambig - 1

    # check arguments
    logging.raise_arg_errors(args, log_file)

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
        preprocessed_alignment,
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

    # generate primer region list
    primer_regions = regions.find_regions(
        ambiguous_consensus,
        args.n_ambig
    )
    if not primer_regions:
        logging.raise_error(
            "no primer regions found. Lower the threshold!",
            log_file,
            exit=True
        )
    logging.varvamp_progress(
        log_file,
        progress=0.4,
        job="Finding primer regions.",
        progress_text=f"{regions.mean(primer_regions, majority_consensus)} % of the consensus sequence will be evaluated for primers"
    )

    # produce kmers for all primer regions
    kmers = regions.produce_kmers(
        primer_regions,
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

    return alignment_cleaned, majority_consensus, ambiguous_consensus, primer_regions, left_primer_candidates, right_primer_candidates


def sanger_and_tiled_shared_workflow(args, left_primer_candidates, right_primer_candidates, log_file):
    """
    part of the workflow shared by the sanger and tiled mode
    """

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
            "no amplicons found. Increase the max amplicon length or \
            number of ambiguous bases or lower threshold!\n",
            log_file,
            exit=True
        )
    logging.varvamp_progress(
        log_file,
        progress=0.8,
        job="Finding potential amplicons.",
        progress_text=f"{len(amplicons)} potential amplicons"
    )
    return all_primers, amplicons


def sanger_workflow(args, amplicons, all_primers, log_file):
    """
    workflow part specific for sanger mode
    """

    amplicon_scheme = scheme.find_sanger_amplicons(amplicons, all_primers, args.report_n)
    logging.varvamp_progress(
        log_file,
        progress=0.9,
        job="Finding low scoring amplicons.",
        progress_text=f"{len(amplicon_scheme[0])} low scoring amplicons."
    )

    return amplicon_scheme


def tiled_workflow(args, amplicons, left_primer_candidates, right_primer_candidates, all_primers, ambiguous_consensus, log_file):
    """
    part of the workflow specific for the tiled mode
    """

    # create graph
    amplicon_graph = scheme.create_amplicon_graph(amplicons, args.overlap)

    # search for amplicon scheme
    coverage, amplicon_scheme = scheme.find_best_covering_scheme(
        amplicons,
        amplicon_graph,
        all_primers
    )

    # check for dimers
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

    # evaluate coverage
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
    return amplicon_scheme


def qpcr_workflow(args, alignment_cleaned, ambiguous_consensus, majority_consensus, left_primer_candidates, right_primer_candidates, log_file):
    """
    part of the workflow specific for the tiled mode
    """
    # find regions for qPCR probes
    probe_regions = regions.find_regions(
        ambiguous_consensus,
        args.pn_ambig
    )
    if not probe_regions:
        logging.raise_error(
            "no regions that fullfill probe criterias! lower threshold or increase number of ambiguous chars in probe\n",
            log_file,
            exit=True
        )

    # digest probe regions
    probe_kmers = regions.produce_kmers(
        probe_regions,
        majority_consensus,
        config.QPROBE_SIZES
    )
    # find potential probes
    qpcr_probes = qpcr.get_qpcr_probes(probe_kmers, ambiguous_consensus, alignment_cleaned)
    if not qpcr_probes:
        logging.raise_error(
            "no qpcr probes found\n",
            log_file,
            exit=True
        )
    logging.varvamp_progress(
        log_file,
        progress=0.7,
        job="Finding qPCR probes.",
        progress_text=f"{len(qpcr_probes)} potential qPCR probes"
    )

    # find unique high scoring amplicons with internal probe
    qpcr_scheme_candidates = qpcr.find_qcr_schemes(qpcr_probes, left_primer_candidates, right_primer_candidates, majority_consensus)
    if not qpcr_scheme_candidates:
        logging.raise_error(
            "no qPCR scheme candidates found. lower threshold or increase number of ambiguous chars in primer and/or probe\n",
            log_file,
            exit=True
        )
    logging.varvamp_progress(
        log_file,
        progress=0.8,
        job="Finding unique amplicons with probe.",
        progress_text=f"{len(qpcr_scheme_candidates)} unique amplicons with internal probe"
    )

    # test amplicons for deltaG
    final_schemes = qpcr.test_amplicon_deltaG(qpcr_scheme_candidates, majority_consensus, args.test_n, args.deltaG)
    if not final_schemes:
        logging.raise_error(
            "no qPCR amplicon passed the deltaG threshold\n",
            log_file,
            exit=True
        )
    logging.varvamp_progress(
        log_file,
        progress=0.9,
        job="Filtering amplicons for deltaG.",
        progress_text=f"{len(final_schemes)} non-overlapping qPCR schemes that passed deltaG cutoff"
    )

    return probe_regions, final_schemes


def main(sysargs=sys.argv[1:]):
    """
    main varvamp workflow
    """

    # start varVAMP
    args = get_args(sysargs)
    if not args.verbose:
        sys.stdout = open(os.devnull, 'w')
    start_time = time.process_time()
    results_dir, data_dir, log_file = logging.create_dir_structure(args.input[1])

    # mode unspecific part of the workflow
    alignment_cleaned, majority_consensus, ambiguous_consensus, primer_regions, left_primer_candidates, right_primer_candidates = shared_workflow(args, log_file)

    # write files that are shared in all modes
    reporting.write_regions_to_bed(primer_regions, data_dir)
    reporting.write_alignment(data_dir, alignment_cleaned)
    reporting.write_fasta(data_dir, "majority_consensus", majority_consensus)
    reporting.write_fasta(results_dir, "ambiguous_consensus", ambiguous_consensus)

    # SANGER/TILED mode
    if args.mode == "tiled" or args.mode == "sanger":
        all_primers, amplicons = sanger_and_tiled_shared_workflow(
            args,
            left_primer_candidates,
            right_primer_candidates,
            log_file
        )
        if args.mode == "sanger":
            amplicon_scheme = sanger_workflow(
                args,
                amplicons,
                all_primers,
                log_file
            )
        elif args.mode == "tiled":
            amplicon_scheme = tiled_workflow(
                args,
                amplicons,
                left_primer_candidates,
                right_primer_candidates,
                all_primers,
                ambiguous_consensus,
                log_file
            )
        # write files
        reporting.write_all_primers(data_dir, all_primers)
        reporting.write_scheme_to_files(
            results_dir,
            amplicon_scheme,
            ambiguous_consensus,
            args.mode
        )
        reporting.varvamp_plot(
            results_dir,
            alignment_cleaned,
            primer_regions,
            all_primers=all_primers,
            amplicon_scheme=amplicon_scheme,
        )
        reporting.per_base_mismatch_plot(results_dir, amplicon_scheme, args.threshold)

    # QPCR mode
    if args.mode == "qpcr":
        probe_regions, final_schemes = qpcr_workflow(
            args,
            alignment_cleaned,
            ambiguous_consensus,
            majority_consensus,
            left_primer_candidates,
            right_primer_candidates,
            log_file
        )
        reporting.write_regions_to_bed(probe_regions, data_dir, "probe")
        reporting.write_qpcr_to_files(results_dir, final_schemes, ambiguous_consensus)
        reporting.varvamp_plot(
            results_dir,
            alignment_cleaned,
            primer_regions,
            probe_regions=probe_regions,
            amplicon_scheme=final_schemes
        )
        reporting.per_base_mismatch_plot(results_dir, final_schemes, args.threshold, mode="QPCR")

    # varVAMP finished
    logging.varvamp_progress(log_file, progress=1, start_time=start_time)
