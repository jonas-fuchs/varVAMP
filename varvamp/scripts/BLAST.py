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
    for strand in all_primers:
        query_path = os.path.join(data_dir, "BLAST_query.fasta")
        with open(query_path, "w") as query:
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
    os.remove(blast_out)

    return(blast_df)


def predict_non_specific_amplicons():
    """
    for a given primer pair, predict unspecific targets within a size
    range and give these primers a high penalty.
    """
    print('TODO')


def write_non_specific_warnings():
    """
    for each primer pair that has potential unspecific amplicons
    write warnings to file.
    """
    print('TODO')