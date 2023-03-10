"""
This contains all varVAMP settings.
"""

# CAN BE CHANGED
# alignment and consensus creation threshold
FREQUENCY_THRESHOLD = 0.9  # freq at which a nucleotide is considered conserved
ALLOWED_N_AMB = 4  # allowed number of ambiguous chars in primer

# basic primer parameters
PRIMER_TMP = (57, 63, 60)  # temperatur (min, max, opt)
PRIMER_GC_RANGE = (40, 60, 50)  # gc (min, max, opt)
PRIMER_SIZES = (18, 27, 20)  # size (min, max, opt)
PRIMER_HAIRPIN = 47  # max melting temp for secondary structures
MAX_POLYX = 5  # max number of polyx
MAX_DINUC_REPEATS = 2  # max number of dinucleotide repeats
MAX_DIMER_TMP = 21  # max melting temp for dimers (homo- or heterodimers)
MIN_3_WITHOUT_AMB = 2  # min len of 3' without ambiguous charaters

# PCR parameters - adjust to your PCR
MV_CONC = 50  # monovalent cations mM
DV_CONC = 2  # divalent cations mM
DNTP_CONC = 0.8  # dntp concentration mM
DNA_CONC = 50  # primer concentration nM

# multipliers for primer base penalties
PRIMER_TM_PENALTY = 2  # temperature penalty
PRIMER_GC_PENALTY = 0.2  # gc penalty
PRIMER_SIZE_PENALTY = 0.5  # size penalty
PRIMER_MAX_BASE_PENALTY = 8  # penalty for primer hardfiltering
PRIMER_3_PENALTY = (10, 10, 10)  # penalties for 3' mismatches
PRIMER_PERMUTATION_PENALTY = 0.1  # penalty for the number of permutations

# amplicon settings
MIN_OVERLAP = 100
OPT_AMPLICON_LENGTH = 1000
MAX_AMPLICON_LENGTH = 2000

# DO NOT CHANGE
# nucleotide definitions
nucs = set("atcg")
ambig_nucs = {
    "r": ["a", "g"],
    "y": ["c", "t"],
    "s": ["g", "c"],
    "w": ["a", "t"],
    "k": ["g", "t"],
    "m": ["a", "c"],
    "b": ["c", "g", "t"],
    "d": ["a", "g", "t"],
    "h": ["a", "c", "t"],
    "v": ["a", "c", "g"],
    "n": ["a", "c", "g", "t"]
}
