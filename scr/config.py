"""
config file for all settings in varVAMP
"""

# PARAMETERS
# params for cleaning deletions in alignment
DELETION_LENGTH_MIN = 1
MASK_LENGTH = 1
# params for consensus creation
FREQUENCY_THRESHOLD = 0.91
# params for conserved region search
CONSERVED_MIN_LENGTH = 19
ALLOWED_N_AMB = 4
# basic primer parameter (min, max, opt)
PRIMER_TMP = (59.5, 62.5, 61)
PRIMER_GC_RANGE = (30,55,50)
PRIMER_SIZE = (19,34,22)
PRIMER_HAIRPIN = 47
MAX_POLYX = 5
MAX_DINUC_REPEATS = 2
MIN_3_WITHOUT_AMB = 2
# PCR parameters
MV_CONC = 100
DV_CONC = 2
DNTP_CONC = 0.8
DNA_CONC = 15
# parameters for stability
MAX_LOOP = 30
TEMP_C = 37
# multipliers for base penalty scoring
PRIMER_TM_PENALTY = 2
PRIMER_GC_PENALTY = 0.2
PRIMER_SIZE_PENALTY = 0.5
PRIMER_3_PENALTY = (8,7,6) # 3' penalty from the 3'
PRIMER_MAX_BASE_PENALTY = 8

# NUCLEOTIDES
nucs = set("atcg")
ambig = {"r": ["a", "g"],
        "y":["c", "t"],
        "s":["g", "c"],
        "w":["a", "t"],
        "k":["g", "t"],
        "m":["a", "c"],
        "b":["c", "g", "t"],
        "d":["a", "g", "t"],
        "h":["a", "c", "t"],
        "v":["a", "c", "g"],
        "n":["a", "c", "g", "t"]}
