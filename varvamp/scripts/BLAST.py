"""
BLAST feature of varVAMP. Checks primers of potential amplicons for hits
against a db and evals if hits are close enough together to potentially
produce off-target amplicons.
"""

##########################   DEV Notes  #######################################

# CHALLENGES/QUESTIONS TO CONSIDER

# code must address the following challenges:
# - well-defined BLAST criteria for off-target effects based on:
#   - BLAST coverage and mismatches
#   - E-value
#   - closeness of two hits
#   --> kind of similar to the primer BLAST feature
# - hard or soft-masking? If no other amplicons than the one producing
#   the off-target amplificate is possible should it be considered or should
#   varvamp fail
# - off-target effects within the consensus genome? so can the primer also
#   bind a different spot? Probably more relevant for larger alignments with
#   repetitive regions
# - error catching if BLAST is installed or not (windoof as well)
# - output cleanup: Clean the output of useless files (BLAST always needs
#   fasta sequences as input). No way (?) to just stay in ram
# - automatic db download if none is provided?
# - write primers to test into files
# - multithreading?

###############################################################################
# BUILT-INS
import subprocess
from sys import platform

#LIBS
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
    if not blast_loc:
        sys.exit("ERROR: BLASTN is not installed")


def create_BLAST_query(all_primers, data_dir):
    """
    create a query for the BLAST search
    """
    for strand in all_primers:
        if strand == "+":
            BLAST_query = os.path.join(data_dir, "BLAST_query_fw.fasta")
        else:
            BLAST_query = os.path.join(data_dir, "BLAST_query_rw.fasta")

        with open(BLAST_query, "w") as query:
            for primer in all_primers[strand]:
                print(f">{primer}\n{all_primers[strand][primer][0]}", file=query)


def run_BLAST(query, blast_db, data_dir):
    """
    runs a BLAST search on a search query.
    """
    outfile = os.path.join(data_dir, f"{query.strip('.fasta')}_result.tabular")
    blast_command = NcbiblastnCommandline(query=query,
                                          db=blast_db,
                                          out=outfile,
                                          task="blastn-short",
                                          num_threads=4,
                                          **config.blast_settings
                                          )

    stdout, stderr = blast_command()

def parse_BLAST_output_to_dictionary():
    """
    create a BLAST hit database for each primer. this can then be used to
    look up if two primers of an amplicons have potential amplificates
    """
    print('TODO')


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