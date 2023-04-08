"""
This contains all varVAMP parameters. Options that can be adjusted by arguments
are FREQUENCY_THRESHOLD, PRIMER_ALLOWED_N_AMB, AMPLICON_MIN_OVERLAP, AMPLICON_OPT_LENGTH,
AMPLICON_MAX_LENGTH.
"""

# CAN BE CHANGED
# basic primer parameters
PRIMER_TMP = (57, 63, 60)  # melting temperatur (min, max, opt)
PRIMER_GC_RANGE = (35, 65, 50)  # gc (min, max, opt)
PRIMER_SIZES = (18, 24, 21)  # size (min, max, opt)
PRIMER_MAX_POLYX = 3  # max number of polyx repeats
PRIMER_MAX_DINUC_REPEATS = 3  # max number of dinucleotide repeats
PRIMER_HAIRPIN = 47  # max melting temp for secondary structures
PRIMER_MAX_GC_END = 3  # max GCs in the last 5 bases of the primer
PRIMER_GC_CLAMP = 1  # min number of GCs in the last 5 bases of the primer
PRIMER_MIN_3_WITHOUT_AMB = 3  # min len of 3' without ambiguous charaters
PRIMER_MAX_DIMER_TMP = 47  # max melting temp for dimers (homo- or heterodimers)

# QPCR probe parameters
QPROBE_TMP = (63, 69, 66) # melting temperatures
QPROBE_SIZES = (20, 30, 25) # size range
QPROBE_GC_RANGE = (30, 70, 50) # GC content
QPROBE_AMB = 1  # n of ambig nucs

# PCR parameters
PCR_MV_CONC = 100  # monovalent cations mM
PCR_DV_CONC = 2  # divalent cations mM
PCR_DNTP_CONC = 0.8  # dntp concentration mM
PCR_DNA_CONC = 15  # primer concentration nM

# multipliers for primer base penalties
PRIMER_TM_PENALTY = 2  # temperature penalty
PRIMER_GC_PENALTY = 0.2  # gc penalty
PRIMER_SIZE_PENALTY = 0.5  # size penalty
PRIMER_MAX_BASE_PENALTY = 8  # max base penalty for a primer
PRIMER_3_PENALTY = (32, 16, 8, 4, 2)  # penalties for 3' mismatches
PRIMER_PERMUTATION_PENALTY = 0.1  # penalty for the number of permutations


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
