"""
varVAMP config file
"""

# PARAMETERS - can be changed
# alignment and consensus creation threshold
FREQUENCY_THRESHOLD = 0.91
# allowed number of ambiguous chars in min primer length
ALLOWED_N_AMB = 4
# basic primer parameter (min, max, opt)
PRIMER_TMP = (59.5, 62.5, 61)
PRIMER_GC_RANGE = (30,55,50)
PRIMER_SIZES = (19,34,22)
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
# multipliers for primer penalties
PRIMER_TM_PENALTY = 2
PRIMER_GC_PENALTY = 0.2
PRIMER_SIZE_PENALTY = 0.5
PRIMER_3_PENALTY = (8,7,6) # 3' penalty from the 3'
PRIMER_MAX_BASE_PENALTY = 8


# nucleotide definitions - do not change
nucs = set("atcg")
ambig_nucs = {
    "r": ["a", "g"],
    "y":["c", "t"],
    "s":["g", "c"],
    "w":["a", "t"],
    "k":["g", "t"],
    "m":["a", "c"],
    "b":["c", "g", "t"],
    "d":["a", "g", "t"],
    "h":["a", "c", "t"],
    "v":["a", "c", "g"],
    "n":["a", "c", "g", "t"]
}
