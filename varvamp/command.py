"""
main workflow
"""

# BUILT-INS
import sys
import os
import shutil
import time
import datetime
import argparse

# varVAMP
from . import _program
from varvamp import __version__
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
    # arg parsing
    parser = argparse.ArgumentParser(
        prog=_program,
        description='varvamp: variable virus amplicon design',
        usage='''varvamp <alignment> <output> [options]''')

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


def create_dir_structure(dir):
    """
    create output folders and log file
    """
    cwd = os.getcwd()
    results_dir = os.path.join(cwd, dir)
    data_dir = os.path.join(results_dir, "data/")
    # create folders
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    else:
        shutil.rmtree(results_dir)
        os.makedirs(results_dir)
    os.makedirs(data_dir)

    log_file = os.path.join(results_dir, "varvamp_log.txt")

    return results_dir, data_dir, log_file


def varvamp_progress(log_file, start_time=None, progress=0, job="", progress_text=""):
    """
    progress bar, main progress logging and folder creation
    """
    barLength = 40
    block = int(round(barLength*progress))

    if progress == 0:
        print(
            "\nStarting \033[31m\033[1mvarVAMP ◥(ºwº)◤\033[0m primer design\n",
            flush=True
        )
        with open(log_file, 'w') as f:
            f.write('VARVAMP log \n\n')
    else:
        if progress == 1:
            stop_time = str(round(time.process_time() - start_time, 2))
            progress_text = f"all done \n\n\rvarVAMP created an amplicon scheme in {stop_time} sec!\n{datetime.datetime.now()}"
            job = "Finalizing output."
        print(
            "\rJob:\t\t " + job + "\nProgress: \t [{0}] {1}%".format("█"*block + "-"*(barLength-block), progress*100) + "\t" + progress_text,
            flush=True
        )
        with open(log_file, 'a') as f:
            print(
                f"\rJob:\t {job} \nResult:\t {progress_text}",
                file=f
            )


def raise_error(message, log_file, exit=False):
    """
    raises warnings or errors, writes to log
    """
    # print to log
    with open(log_file, 'a') as f:
        if exit:
            print(f"ERROR: {message}", file=f)
        else:
            print(f"WARNING: {message}", file=f)
    # print to console
    if exit:
        sys.exit(f"\n\033[31m\033[1mERROR:\033[0m {message}")
    else:
        print(f"\033[31m\033[1mWARNING:\033[0m {message}")


def raise_arg_errors(args, log_file):
    """
    checks arguments for non-valid input and raises warnings
    """
    # threshold error
    if args.threshold > 1 or args.threshold < 0:
        raise_error(
            "threshold can only be between 0-1",
            log_file,
            exit=True
        )
    if args.allowed_ambiguous < 0:
        raise_error(
            "number of ambiguous chars can not be negative",
            log_file,
            exit=True
        )
    if args.allowed_ambiguous > 4:
        raise_error(
            "high number of ambiguous nucleotides in primer leads to a high "
            "degeneracy. Consider reducing.",
            log_file
        )
    if args.opt_length > args.max_length:
        raise_error(
            "optimal length can not be higher than the maximum amplicon length.",
            log_file,
            exit=True
        )
    if args.opt_length < 0 or args.max_length < 0:
        raise_error(
            "amplicon lengths can not be negative.",
            log_file,
            exit=True
        )
    if args.opt_length < 200 or args.max_length < 200:
        raise_error(
            "your amplicon lengths might be to small. Consider increasing",
            log_file
        )
    if args.overlap < 0:
        raise_error(
            "overlap size can not be negative.",
            log_file,
            exit=True
        )
    if args.overlap < 50:
        raise_error(
            "small overlaps might hinder downstream analyses. Consider increasing.",
            log_file
        )
    if args.overlap > args.max_length/2 - config.PRIMER_SIZES[1]:
        raise_error(
            "min overlap must be lower than half of your maximum length - maximum primer length. To achieve optimal results reduce it to at least half of your optimal length",
            log_file,
            exit=True
        )
    if args.overlap > args.opt_length:
        raise_error(
            "overlap can not be higher than your optimal length.",
            log_file,
            exit=True
        )
    if args.overlap > args.opt_length/2:
        raise_error(
            "your intended overlap is higher than half of your optimal length. This reduces how well varvamps will find overlapping amplicons. Consider decreasing.",
            log_file
        )


def confirm_config(args, log_file):
    """
    checks the config. raises error and warnings
    if nececarry. writes settings to log
    """
    error = False

    # check if all variables exists
    all_vars = [
        # arg dependent
        "FREQUENCY_THRESHOLD",
        "PRIMER_ALLOWED_N_AMB",
        "AMPLICON_OPT_LENGTH",
        "AMPLICON_MAX_LENGTH",
        "AMPLICON_MIN_OVERLAP",
        # arg independent
        "PRIMER_TMP",
        "PRIMER_GC_RANGE",
        "PRIMER_SIZES",
        "PRIMER_HAIRPIN",
        "PRIMER_MAX_POLYX",
        "PRIMER_MAX_DINUC_REPEATS",
        "PRIMER_MAX_DIMER_TMP",
        "PRIMER_MIN_3_WITHOUT_AMB",
        "PCR_MV_CONC",
        "PCR_DV_CONC",
        "PCR_DNTP_CONC",
        "PCR_DNA_CONC",
        "PRIMER_TM_PENALTY",
        "PRIMER_GC_PENALTY",
        "PRIMER_SIZE_PENALTY",
        "PRIMER_MAX_BASE_PENALTY",
        "PRIMER_3_PENALTY",
        "PRIMER_PERMUTATION_PENALTY",
    ]

    for var in all_vars:
        if var not in vars(config):
            raise_error(
                f"{var} does not exist in config!",
                log_file
            )
            error = True
    # exit if variables are not defined
    if error:
        raise_error(
            "config is missing parameters. Look at the above warnings!",
            log_file,
            exit=True
        )
    # confirm tuples
    for type, tup in [("temp", config.PRIMER_TMP), ("gc", config.PRIMER_GC_RANGE), ("size", config.PRIMER_SIZES)]:
        if len(tup) != 3:
            raise_error(
                f"{type} tuple has to have the form (min, max, opt)!",
                log_file
            )
            error = True
        if tup[0] > tup[1]:
            raise_error(
                f"min {type} should not exeed max {type}!",
                log_file
            )
            error = True
        if tup[0] > tup[2]:
            raise_error(
                f"min {type} should not exeed opt {type}!",
                log_file
            )
            error = True
        if tup[2] > tup[1]:
            raise_error(
                f"opt {type} should not exeed max {type}!",
                log_file
            )
            error = True
        if any(map(lambda var: var < 0, tup)):
            raise_error(
                f"{type} can not contain negative values!",
                log_file
            )
            error = True

    # check values that cannot be zero
    non_negative_var = [
        ("max polyx repeats", config.PRIMER_MAX_POLYX),
        ("max dinucleotide repeats", config.PRIMER_MAX_DINUC_REPEATS),
        ("max GCs at the 3' end", config.PRIMER_MAX_GC_END),
        ("GC clamp", config.PRIMER_GC_CLAMP),
        ("min number of 3 prime nucleotides without ambiguous nucleotides", config.PRIMER_MIN_3_WITHOUT_AMB),
        ("monovalent cation concentration", config.PCR_MV_CONC),
        ("divalent cation concentration", config.PCR_DV_CONC),
        ("dNTP concentration", config.PCR_DNTP_CONC),
        ("primer temperatur penalty", config.PRIMER_TM_PENALTY),
        ("primer gc penalty", config.PRIMER_GC_PENALTY),
        ("primer size penalty", config.PRIMER_SIZE_PENALTY),
        ("max base penalty", config.PRIMER_MAX_BASE_PENALTY),
        ("primer permutation penalty", config.PRIMER_PERMUTATION_PENALTY)
    ]
    for type, var in non_negative_var:
        if var < 0:
            raise_error(
                f"{type} can not be negative!",
                log_file
            )
            error = True
    if any(map(lambda var: var < 0, config.PRIMER_3_PENALTY)):
        raise_error(
            "3' penalties can not be zero!",
            log_file
        )
        error = True
    # exit if variables are not properly defined
    if error:
        raise_error(
            "config has flaws. Look at the above warnings!",
            log_file,
            exit=True
        )
    # specific warnings
    if config.PRIMER_HAIRPIN < 0:
        raise_error(
            "decreasing hairpin melting temp to negative values "
            "will influence successful primer search!",
            log_file
        )
    if config.PRIMER_MAX_DIMER_TMP < 0:
        raise_error(
            "there is no need to set max dimer melting temp below 0.",
            log_file
        )
    if config.PRIMER_MAX_BASE_PENALTY < 8:
        raise_error(
            "decreasing the base penalty will hardfilter more primers.",
            log_file
        )
    if config.PRIMER_GC_CLAMP > 3:
        raise_error(
            "large GC clamps will results in too high 3'end stability",
            log_file
        )
    if config.PRIMER_MAX_GC_END < 5 and config.PRIMER_MAX_GC_END < config.PRIMER_GC_CLAMP:
        raise_error(
            f"GC clamp of {config.PRIMER_GC_CLAMP} length will not be enforced as there are only {config.PRIMER_MAX_GC_END} gc characters allowed at the 3' end",
            log_file
        )
    if config.PRIMER_MAX_GC_END > 5:
        raise_error(
            "only the last 5 nucleotides of the 3' end are considered for GC 3'end calculation.",
            log_file
        )

    # write all settings to file
    var_dic = vars(config)
    with open(log_file, 'a') as f:
        print(
            "settings that can be adjusted via arguments\n",
            f"PRIMER_OPT_LENGTH = {args.opt_length}",
            f"PRIMER_MAX_LENGTH = {args.max_length}",
            f"PRIMER_MIN_OVERLAP = {args.overlap}",
            f"PRIMER_THRESHOLD = {args.threshold}",
            f"PRIMER_ALLOWED_N_AMB = {args.allowed_ambiguous}",
            "\nconfig settings\n",
            sep="\n",
            file=f
        )
        for var in all_vars[5:]:
            print(f"{var} = {var_dic[var]}", file=f)


def main(sysargs=sys.argv[1:]):
    """
    main varvamp workflow
    """

    # start varVAMP
    # - arg parsing
    args = get_args(sysargs)
    # - supress console output
    if not args.console:
        sys.stdout = open(os.devnull, 'w')
    # - ini time
    start_time = time.process_time()
    # - create folder paths
    results_dir, data_dir, log_file = create_dir_structure(args.input[1])
    # raise arg errors
    raise_arg_errors(args, log_file)
    # - update progress
    varvamp_progress(log_file)

    # config check
    confirm_config(args, log_file)
    # - progress update
    varvamp_progress(
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
    # - write alignment
    reporting.write_alignment(data_dir, alignment_cleaned)
    # - progress update
    varvamp_progress(
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
    # - write consensus sequence
    reporting.write_fasta(data_dir, "majority_consensus", majority_consensus)
    reporting.write_fasta(results_dir, "ambiguous_consensus", ambiguous_consensus)
    # - progress update
    varvamp_progress(
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
    # - raise error if no conserved regions were found
    if not conserved_regions:
        raise_error(
            "nothing conserved. Lower the threshold!",
            log_file,
            exit=True
        )
    # - write conserved regions to bed file
    reporting.write_conserved_to_bed(conserved_regions, data_dir)
    # - progress update
    varvamp_progress(
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
    # - progress update
    varvamp_progress(
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
    # - raise error if no primers were found
    for type, primer_candidates in [("LEFT", left_primer_candidates), ("RIGHT", right_primer_candidates)]:
        if not primer_candidates:
            raise_error(
                f"no {type} primers found.\n",
                log_file,
                exit=True
            )
    # - progress update
    varvamp_progress(
        log_file,
        progress=0.6,
        job="Filtering for primers.",
        progress_text=f"{len(left_primer_candidates)} fw and {len(right_primer_candidates)} rw potential primers"
    )

    # - find best primers and create primer dict
    left_primer_candidates = primers.find_best_primers(left_primer_candidates, "LEFT")
    right_primer_candidates = primers.find_best_primers(right_primer_candidates, "RIGHT")
    # - write primers
    reporting.write_all_primers(data_dir, left_primer_candidates, right_primer_candidates)
    # - progress update
    varvamp_progress(
        log_file,
        progress=0.7,
        job="Considering only high scoring primers.",
        progress_text=f"{len(left_primer_candidates)} fw and {len(right_primer_candidates)} rw primers"
    )

    # find all possible amplicons
    amplicons = scheme.find_amplicons(
        left_primer_candidates,
        right_primer_candidates,
        args.opt_length,
        args.max_length
    )
    # - raise error if no amplicons were found
    if not amplicons:
        raise_error(
            "no amplicons found. Increase the max "
            "amplicon length or lower threshold!\n",
            log_file,
            exit=True
        )
    # - build the amplicon graph
    amplicon_graph = scheme.create_amplicon_graph(amplicons, args.overlap)
    # - progress update
    varvamp_progress(
        log_file,
        progress=0.8,
        job="Finding potential amplicons.",
        progress_text=str(len(amplicons)) + " potential amplicons"
    )

    # search for amplicon scheme
    coverage, amplicon_scheme = scheme.find_best_covering_scheme(
        amplicons,
        amplicon_graph
    )
    percent_coverage = round(coverage/len(ambiguous_consensus)*100, 2)
    # - write all relevant files for the scheme
    reporting.write_scheme_to_files(
        results_dir,
        amplicon_scheme,
        amplicons,
        ambiguous_consensus,
        left_primer_candidates,
        right_primer_candidates
    )
    # - progress update
    varvamp_progress(
        log_file,
        progress=0.9,
        job="Creating amplicon scheme.",
        progress_text=f"{percent_coverage} % total coverage with {len(amplicon_scheme)} amplicons"
    )
    # - raise low coverage warning
    if percent_coverage < 70:
        raise_error(
            "coverage < 70 %. Possible solutions:\n"
            "\t - lower threshold\n"
            "\t - increase amplicons lengths\n"
            "\t - increase number of ambiguous nucleotides\n"
            "\t - relax primer settings (not recommended)\n",
            log_file
        )

    # plotting
    reporting.varvamp_plot(
        results_dir,
        args.threshold,
        alignment_cleaned,
        conserved_regions,
        amplicon_scheme,
        amplicons,
        left_primer_candidates,
        right_primer_candidates
    )
    # - final progress
    varvamp_progress(log_file, progress=1, start_time=start_time)
