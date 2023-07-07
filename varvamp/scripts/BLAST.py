"""
BLAST feature of varVAMP. Checks primers of potential amplicons for hits
against a db and evals if hits are close enough together to potentially
produce off-target amplicons.
"""

##########################   DEV Notes  #######################################

# CHALLENGES/QUESTIONS TO CONSIDER

# code must adress the following challenges:
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


# NEEDED FUNCTIONS

# 'def write_primers_to_fasta'
# 'def clean_BLAST_output'
# 'def BLAST_primers'
# 'def parse_BLAST_output'
# 'def find_off_targets'
# 'def mask_primers'
# 'def write_off_target_warnings'
# (...)

###############################################################################