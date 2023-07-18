"""
BLAST feature of varVAMP. Checks primers of potential amplicons for hits
against a db and evals if hits are in close proximity to potentially
produce off-target amplicons.
"""

##########################   Next steps  #####################################

# - check for amplicons unspecific amplification (max size 2x amplicon length)
# - if so add a high penalty for this amplicon (+100 (?) - which means it will
# be used if no other amplicon is available). in this case write a warning to
# log file

###############################################################################
# BUILT-INS
import subprocess
import os
import sys
from sys import platform

#LIBS
import pandas as pd
from Bio.Blast.Applications import NcbiblastnCommandline

# varVAMP
from varvamp.scripts import config

def check_BLAST_installation():
    """
     simple check if BLAST is installed
    """
    if platform == "win32":
        blast_loc = subprocess.getoutput("where blastn")
    else:
        blast_loc = subprocess.getoutput("which blastn")
    if blast_loc:
        print("BLASTN is installed. Continuing...")
    else:
        sys.exit("ERROR: BLASTN is not installed")


def create_BLAST_query(all_primers, data_dir):
    """
    create a query for the BLAST search
    """
    query_path = os.path.join(data_dir, "BLAST_query.fasta")
    with open(query_path, "w") as query:
        for strand in all_primers:
            for primer in all_primers[strand]:
                print(f">{primer}\n{all_primers[strand][primer][0]}", file=query)

    return query_path

def run_BLAST(query, blast_db, data_dir):
    """
    runs a BLAST search on a search query.
    """
    basename = os.path.basename(os.path.normpath(query))
    blast_out = os.path.join(data_dir, f"{basename.strip('.fasta')}_result.tabular")

    blast_command = NcbiblastnCommandline(
        query=query,
        db=blast_db,
        out=blast_out,
        task="blastn-short",
        num_threads=4,
        **config.BLAST_SETTINGS
        )
    # run BLAST
    stdout, stderr = blast_command()
    # remove query input
    os.remove(query)

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
               ]

    blast_df = pd.read_table(blast_out, names=columns)
    blast_df = blast_df[
        blast_df["query_len"] - blast_df["aln_len"] + blast_df["mismatch"] + blast_df["gaps"] <= round(
            blast_df["query_len"] * config.BLAST_MAX_DIFF)
        ]
    # removes tabular output
    os.remove(blast_out)

    return(blast_df)


def predict_non_specific_amplicons(amplicons, blast_df, max_length):
    """
    for a given primer pair, predict unspecific targets within a size
    range and give these primers a high penalty.
    """
    off_target_amp = []

    for amp in amplicons:
        fw_primer = amplicons[amp][2]
        rw_primer = amplicons[amp][3]
        # subset df for primers
        df_amp_primers = blast_df[blast_df["query"].isin([fw_primer, rw_primer])]
        # sort by reference and ref start (so one can easily iterate over rows)
        df_amp_primers_sorted = df_amp_primers.sort_values(["ref", "ref_start"])
        # iterate over ref for each primer pair
        for ref in set(df_amp_primers_sorted["ref"]):
            df_ref_subset = df_amp_primers_sorted[df_amp_primers_sorted["ref"] == ref]
            # reindex to remember the correct index
            df_ref_subset.reset_index(inplace=True, drop=True)
            # ini the primer search
            indices = []  # remember the index of the subset df
            start = -float("inf")
            off_target = False
            for i, row in df_ref_subset.iterrows():
                # for the current row check if start of the current batch is close enough
                if row[7] - start <= config.BLAST_SIZE_MULTI * max_length:
                    indices.append(i)
                else:
                    # reset start to the next element in the index list if list has more than
                    # one element
                    if len(indices) > 1:
                        indices.pop(0)
                        start = df_ref_subset.iloc[indices[0]]["ref_start"]
                    # else reset to the current row and empty index list
                    else:
                        start = row[6]
                        indices = []
                # do we need to check the current index list?
                if len(indices) <= 1:
                    continue
                # check if in the df subset is more than one query
                if len(set(df_ref_subset.iloc[indices]["query"])) > 1:
                    off_target = True
                    break
            if off_target:
                amplicons[amp][5] = amplicons[amp][5] + config.BLAST_PENALTY
                off_target_amp.append(amp)
                break

    return(off_target_amp, amplicons)


def write_non_specific_warnings():
    """
    for each primer pair that has potential unspecific amplicons
    write warnings to file.
    """
    print('TODO')