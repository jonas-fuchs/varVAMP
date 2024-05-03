"""
BLAST feature of varVAMP. Checks primers of potential amplicons for hits
against a db and evals if hits are in proximity to potentially
produce off-target amplicons.
"""

# BUILT-INS
import itertools
import multiprocessing
import os
import subprocess
from shutil import which

#LIBS
import pandas as pd

# varVAMP
from varvamp.scripts import config
from varvamp.scripts import logging


def check_BLAST_installation(log_file):
    """
     simple check if BLAST is installed
    """
    if which("blastn") is not None:
        print("\nINFO: BLASTN is installed.")
    else:
        logging.raise_error("BLASTN is not installed", log_file, exit=True)


def create_BLAST_query(amplicons, data_dir, mode="single_tiled"):
    """
    create a query for the BLAST search
    """
    query_path = os.path.join(data_dir, "BLAST_query.fasta")
    if mode == "single_tiled":
        primer_types = ["LEFT", "RIGHT"]
    elif mode == "qpcr":
        primer_types = ["PROBE", "LEFT", "RIGHT"]
    already_written = set()

    with open(query_path, "w") as query:
        for amp in amplicons:
            for primer_type in primer_types:
                name = f"{primer_type}_{amp[primer_type][1]}_{amp[primer_type][2]}"
                if name not in already_written:
                    print(f">{name}\n{amp[primer_type][0]}", file=query)
                    already_written.add(name)
    return query_path


def run_BLAST(query, blast_db, data_dir, n_threads):
    """
    runs a BLAST search on a search query.
    """
    basename = os.path.basename(os.path.normpath(query))
    blast_out = os.path.join(data_dir, f"{basename.strip('.fasta')}_result.tabular")

    blast_command_options = {
        "query": query,
        "db": blast_db,
        "out": blast_out,
        "task": "blastn-short",
        "num_threads": n_threads,
    }
    blast_command_options.update(config.BLAST_SETTINGS)
    args = ['blastn']
    for k, v in blast_command_options.items():
        args.append('-' + k)
        args.append(str(v))
    # run BLAST
    subprocess.run(args, capture_output=True, check=True)

    return blast_out


def parse_and_filter_BLAST_output(blast_out):
    """
    create a BLAST hit database for each primer. filter for mismatches.
    returns a prefiltered pandas df
    """
    columns = ["query",
               "ref",
               "query_len",
               "aln_len",
               "mismatch",
               "gaps",
               "ref_start",
               "ref_end",
               "strand"
               ]

    blast_df = pd.read_table(blast_out, names=columns)
    blast_df = blast_df[
        blast_df["query_len"] - blast_df["aln_len"] + blast_df["mismatch"] + blast_df["gaps"] <=
        blast_df["query_len"]-(round(blast_df["query_len"] * config.BLAST_MAX_DIFF))
        ]
    # removes tabular output
    os.remove(blast_out)

    return blast_df


def check_off_targets(df_amp_primers_sorted, max_length, primers):
    """
    determines off targets in subset blast df
    """
    combinations = list(itertools.combinations(primers, 2))

    for ref in set(df_amp_primers_sorted["ref"]):
        df_ref_subset = df_amp_primers_sorted[df_amp_primers_sorted["ref"] == ref]
        # reindex to remember the correct index later on
        df_ref_subset.reset_index(inplace=True, drop=True)
        # ini the primer search
        indices = []  # remember the index of the subset df
        start = -float("inf")
        for row in df_ref_subset.itertuples():
            # for the current row check if start of the current batch is close enough
            if row[8] - start <= config.BLAST_SIZE_MULTI * max_length:
                indices.append(row[0])
            else:
                # reset start to the next element in the index list if list has more than
                # one element
                if len(indices) > 1:
                    indices.pop(0)
                    start = df_ref_subset.iloc[indices[0]]["ref_start"]
                # else reset to the current row and empty index list
                else:
                    start = row[7]
                    indices = []
            # check if in the df subset is more than one primer
            if len(set(df_ref_subset.iloc[indices]["query"])) <= 1:
                continue
            # subset of df with potential off-targets
            possible_off_targets = df_ref_subset.iloc[indices]
            # check if any two primers in scheme bind different directions
            for combi in combinations:
                direction_fw = set(possible_off_targets[possible_off_targets["query"] == combi[0]]["strand"])
                direction_rv = set(possible_off_targets[possible_off_targets["query"] == combi[1]]["strand"])
                # do we have one primer in subset that binds both directions or are the
                # direction sets different? --> 2 different primers bind on the same chrom,
                # are close enough and are in opposite direction ->[seq]<-
                # >this classifies as an off-target<
                if len(direction_fw) > 1 or len(direction_rv) > 1 or direction_fw != direction_rv:
                    return True
    return False


def predict_non_specific_amplicons_worker(amp, blast_df, max_length, mode):
    """
    Worker function to predict unspecific targets for a single amplicon.
    """
    # get correct primers
    if mode == "single_tiled":
        primer_types = ["LEFT", "RIGHT"]
    elif mode == "qpcr":
        primer_types = ["PROBE", "LEFT", "RIGHT"]
    primers = []
    for primer_type in primer_types:
        primers.append(f"{primer_type}_{amp[primer_type][1]}_{amp[primer_type][2]}")
    # subset df for primers
    df_amp_primers = blast_df[blast_df["query"].isin(primers)]
    # sort by reference and ref start
    df_amp_primers_sorted = df_amp_primers.sort_values(["ref", "ref_start"])
    # check for off-targets for specific primers
    if check_off_targets(df_amp_primers_sorted, max_length, primers):
        amp["off_targets"] = True
    else:
        amp["off_targets"] = False
    return amp


def predict_non_specific_amplicons(amplicons, blast_df, max_length, mode, n_threads):
    """
    Main function to predict unspecific targets within a size range and give
    these primers a high penalty. Uses multiprocessing for parallelization.
    """
    # process amplicons concurrently
    with multiprocessing.Pool(processes=n_threads) as pool:
        annotated_amps = [
            result for result in pool.starmap(
                predict_non_specific_amplicons_worker,
                [(amp, blast_df, max_length, mode) for amp in amplicons]
            ) if result is not None
        ]
    n_off_targets = sum(amp["off_targets"] for amp in annotated_amps)
    return n_off_targets, annotated_amps


def primer_blast(data_dir, db, query_path, amplicons, max_length, n_threads, log_file, mode):
    """
    performs the blast search for the single or tiled workflow
    """
    print("\n#### Starting varVAMP primerBLAST. ####\n")
    print("Running BLASTN...")
    try:
        blast_out = run_BLAST(
            query_path,
            db,
            data_dir,
            n_threads
        )
    except subprocess.CalledProcessError as e:
        logging.raise_error(
            "Failed to run BLAST! Error message: " + str(
                e.stderr
            ) + " Exit code: " + str(
                e.returncode
            ),
            log_file,
            exit=True
        )
    finally:
        # remove query input
        os.remove(query_path)

    blast_df = parse_and_filter_BLAST_output(blast_out)
    print("Predicting non-specific amplicons...")
    n_off_targets, amplicons = predict_non_specific_amplicons(
        amplicons,
        blast_df,
        max_length,
        mode,
        n_threads
    )
    if n_off_targets > 0:
        success_text = f"varVAMP predicted non-specific amplicons:\n\t> {n_off_targets}/{len(amplicons)} amplicons could produce amplicons with the blast db.\n\t> will attempt to avoid them in the final list of amplicons"
    else:
        success_text = f"NO off-target amplicons found with the blast db and a total of {len(amplicons)} amplicons"
    print(success_text)
    with open(log_file, 'a') as f:
        print(
            f"\nBLAST results: {success_text}",
            file=f
        )
    print("\n#### off-target search finished ####\n")

    return amplicons

