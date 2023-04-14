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
from varvamp.scripts import consensus
from varvamp.scripts import conserved
from varvamp.scripts import logging
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
            default=0.89,
            help="threshold for conserved nucleotides"
        )
        par.add_argument(
            "-a",
            "--n-ambig",
            metavar="",
            type=int,
            default=4,
            help="max number of ambiguous characters in a primer"
        )
    for par in (SANGER_parser, TILED_parser):
        par.add_argument(
            "-ol",
            "--opt-length",
            help="optimal length of the amplicons",
            metavar="",
            type=int,
            default=1000
        )
        par.add_argument(
            "-ml",
            "--max-length",
            help="max length of the amplicons",
            metavar="",
            type=int,
            default=1500
        )
    TILED_parser.add_argument(
        "-o",
        "--overlap",
        type=int,
        metavar="",
        default=100,
        help="min overlap of the amplicons"
    )
    SANGER_parser.add_argument(
        "-n",
        "--report-n",
        type=int,
        metavar="",
        default=float("inf"),
        help="report the top n best hits"
    )
    QPCR_parser.add_argument(
        "-pa",
        "--pn-ambig",
        metavar="",
        type=int,
        default=1,
        help="max number of ambiguous characters in a probe"
    )
    QPCR_parser.add_argument(
        "-n",
        "--test-n",
        type=int,
        metavar="",
        default=50,
        help="test the top n qPCR amplicons for secondary structures at the minimal primer temperature"
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
        args.n_ambig
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
        job="Finding conserved primer regions.",
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
    # the next part is specific for either SANGER/TILED mode
    if args.mode == "tiled" or args.mode == "sanger":
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

        if args.mode == "tiled":
            amplicon_graph = scheme.create_amplicon_graph(amplicons, args.overlap)
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
        elif args.mode == "sanger":
            amplicon_scheme = scheme.find_sanger_amplicons(amplicons, all_primers, args.report_n)
            logging.varvamp_progress(
                log_file,
                progress=0.9,
                job="Finding low scoring amplicons.",
                progress_text="The lowest scoring amplicon is amplicon_0."
            )

        # write  mode specific files
        reporting.write_all_primers(data_dir, all_primers)
        reporting.write_scheme_to_files(
            results_dir,
            amplicon_scheme,
            ambiguous_consensus,
            args.mode
        )
        reporting.varvamp_plot(
            results_dir,
            args.threshold,
            alignment_cleaned,
            conserved_regions,
            all_primers,
            amplicon_scheme,
        )
    # specific for QPCR mode
    if args.mode == "qpcr":
        # find regions for qPCR probes
        probe_conserved_regions = conserved.find_regions(
            ambiguous_consensus,
            args.pn_ambig
        )
        if not probe_conserved_regions:
            logging.raise_error(
                "no conserved regions that fullfill probe criterias! lower threshold or increase number of ambiguous chars in probe\n",
                log_file,
                exit=True
            )
        # digest probe regions
        probe_kmers = conserved.produce_kmers(
            probe_conserved_regions,
            majority_consensus
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
            progress_text=f"{len(qpcr_probes)} potential qPCR probes with max {args.pn_ambig} ambiguous chars"
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
        final_schemes = qpcr.test_amplicon_deltaG(qpcr_scheme_candidates, majority_consensus, args.test_n)
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
    # write files that are shared in all modes
    reporting.write_conserved_to_bed(conserved_regions, data_dir)
    reporting.write_alignment(data_dir, alignment_cleaned)
    reporting.write_fasta(data_dir, "majority_consensus", majority_consensus)
    reporting.write_fasta(results_dir, "ambiguous_consensus", ambiguous_consensus)
    logging.varvamp_progress(log_file, progress=1, start_time=start_time)
