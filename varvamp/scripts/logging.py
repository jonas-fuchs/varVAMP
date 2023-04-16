"""
varVAMP logging and raising errors
"""

# BUILT-INS
import sys
import os
import shutil
import time
import datetime

# varVAMP
from varvamp.scripts import config


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
            progress_text = f"all done \n\n\rvarVAMP finished in {stop_time} sec!\n{datetime.datetime.now()}"
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
    if args.n_ambig < 0:
        raise_error(
            "number of ambiguous chars can not be negative",
            log_file,
            exit=True
        )
    if args.n_ambig > 4:
        raise_error(
            "high number of ambiguous nucleotides in primer leads to a high "
            "degeneracy. Consider reducing.",
            log_file
        )
    if args.mode in ("tiled", "sanger"):
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
    # SANGER specific warnings
    if args.mode == "sanger":
        if args.report_n < 1:
            raise_error(
                "number of reported amplicons cannot be below 1.",
                log_file,
                exit=True
            )
    # TILED specific warnings
    if args.mode == "tiled":
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
        # QPCR specific warnings
        if args.mode == "qpcr":
            if args.pn_ambig < 0:
                raise_error(
                    "number of ambiguous characters in the qPCR probe cannot be 0.",
                    log_file,
                    exit=True
                )
            if args.pn_ambig > args.n_ambig:
                raise_error(
                    "the degeneracy of qPCR probes should not be higher than that of qPCR primers to retain specificity",
                    log_file,
                    exit=True
                )
            if args.test_n <= 0:
                raise_error(
                    "if the number of tested qPCR amplicons is 0 or lower, no amplicons will be reported.",
                    log_file,
                    exit=True
                )

            if args.test_n > 200:
                raise_error(
                    "checking the deltaG of amplicons is computationally intensive. setting this higher than 200 will likely take a bit to compute.",
                    log_file
                )
            if args.pn_ambig > 2:
                raise_error(
                    "setting the degeneracy of probes higher than 2 can result in variable probe design. Consider reducing!",
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
        # arg independent TILED, SANGER mode
        (
            "PRIMER_TMP",
            "PRIMER_GC_RANGE",
            "PRIMER_SIZES",
            "PRIMER_HAIRPIN",
            "PRIMER_MAX_POLYX",
            "PRIMER_MAX_DINUC_REPEATS",
            "PRIMER_GC_END",
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
            "PRIMER_PERMUTATION_PENALTY"
        ),
        # arg independent QPCR mode
        (
            "QPROBE_TMP",
            "QPROBE_SIZES",
            "QPROBE_GC_RANGE",
            "QPROBE_GC_END",
            "QPRIMER_DIFF",
            "QPROBE_TEMP_DIFF",
            "QPROBE_DISTANCE",
            "QAMPLICON_LENGTH",
            "QAMPLICON_GC",
            "QAMPLICON_DELTAG_CUTOFF"
        )
    ]

    for tup in all_vars:
        for var in tup:
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
    # confirm tuples with 3 elements
    for type, tup in [("primer temp", config.PRIMER_TMP), ("primer gc", config.PRIMER_GC_RANGE), ("primer size", config.PRIMER_SIZES), ("probe temp", config.QPROBE_TMP), ("probe size", config.QPROBE_SIZES), ("probe gc", config.QPROBE_GC_RANGE)]:
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
    # confirm tuples with two elements
    for type, tup in [("primer gc at the 3'", config.PRIMER_GC_END), ("probe and primer temp diff", config.QPROBE_TEMP_DIFF), ("distance probe to primer", config.QPROBE_DISTANCE), ("qpcr amplicon length", config.QAMPLICON_LENGTH), ("qpcr amplicon gc", config.QAMPLICON_GC), ("probe gc at the 3'", config.PRIMER_GC_END)]:
        if len(tup) != 2:
            raise_error(
                f"{type} tuple has to have the form (min, max)!",
                log_file
            )
            error = True
        if tup[0] > tup[1]:
            raise_error(
                f"min {type} should not exeed max {type}!",
                log_file
            )
            error = True
        if any(map(lambda var: var < 0, tup)):
            raise_error(
                f"{type} can not contain negative values!",
                log_file
            )
            error = True
    # check single values that cannot be negative
    non_negative_var = [
        ("max polyx repeats", config.PRIMER_MAX_POLYX),
        ("max dinucleotide repeats", config.PRIMER_MAX_DINUC_REPEATS),
        ("min number of 3 prime nucleotides without ambiguous nucleotides", config.PRIMER_MIN_3_WITHOUT_AMB),
        ("monovalent cation concentration", config.PCR_MV_CONC),
        ("divalent cation concentration", config.PCR_DV_CONC),
        ("dNTP concentration", config.PCR_DNTP_CONC),
        ("primer temperatur penalty", config.PRIMER_TM_PENALTY),
        ("primer gc penalty", config.PRIMER_GC_PENALTY),
        ("primer size penalty", config.PRIMER_SIZE_PENALTY),
        ("max base penalty", config.PRIMER_MAX_BASE_PENALTY),
        ("primer permutation penalty", config.PRIMER_PERMUTATION_PENALTY),
        ("qpcr flanking primer difference", config.QPRIMER_DIFF),
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
            "decreasing the base penalty will filter out more primers.",
            log_file
        )
    if config.PRIMER_GC_END[0] > 3 or config.QPROBE_GC_END[0] > 3:
        raise_error(
            "large GC clamps will results in too high 3'end stability",
            log_file
        )
    if config.PRIMER_GC_END[1] > 5 or config.QPROBE_GC_END[1] > 5:
        raise_error(
            "only the last 5 nucleotides of the 3' end are considered for GC 3'end calculation.",
            log_file
        )
    if config.QAMPLICON_DELTAG_CUTOFF > 0:
        raise_error(
            "there is no need to set the deltaG cutoff higher than 0 (0 will already avoid any secondary structures at the primer melting temp).",
            log_file
        )

    # write all settings to file
    var_dic = vars(config)
    with open(log_file, 'a') as f:
        print(
            f"MODE = {args.mode}",
            sep="\n",
            file=f
        )
        print(
            "\nsettings that can be adjusted via arguments\n",
            f"THRESHOLD = {args.threshold}",
            f"PRIMER_ALLOWED_N_AMB = {args.n_ambig}",
            sep="\n",
            file=f
        )
        if args.mode in ("tiled", "sanger"):
            print(
                f"AMPLICON_OPT_LENGTH = {args.opt_length}",
                f"AMPLICON_MAX_LENGTH = {args.max_length}",
                sep="\n",
                file=f
            )
        if args.mode == "tiled":
            print(
                f"MIN_OVERLAP = {args.overlap}",
                sep="\n",
                file=f
            )
        if args.mode == "sanger":
            print(
                f"REPORT_N_AMPLICONS = {args.report_n}",
                sep="\n",
                file=f
            )
        if args.mode == "qpcr":
            print(
                f"PROBE_ALLOWED_N_AMB = {args.pn_ambig}",
                f"TEST_DELTAG_N_AMPLICONS = {args.test_n}",
                sep="\n",
                file=f
            )
        print(
            "\nconfig settings\n",
            sep="\n",
            file=f
        )
        for var in all_vars[0]:
            print(f"{var} = {var_dic[var]}", file=f)
        if args.mode == "qpcr":
            for var in all_vars[1]:
                print(f"{var} = {var_dic[var]}", file=f)
        print("\nprogress", file=f)
