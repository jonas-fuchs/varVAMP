"""
This contains all varVAMP parameters.
"""

# CAN BE CHANGED, DO NOT DELETE
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

# QPCR parameters
# base probe parameters to consider
QPROBE_TMP = (64, 70, 67)  # mean 7°C higher than the primer temp
QPROBE_SIZES = (20, 30, 25)
QPROBE_GC_RANGE = (40, 80, 60)
QPROBE_MAX_GC_END = 4
QPROBE_GC_CLAMP = 0
# constraints for amplicon design
QPRIMER_DIFF = 2  # maximal temperature diff of qPCR primers
QPROBE_TEMP_DIFF = (5, 10)  # min/max temp diff between probe and primers
QPROBE_DISTANCE = (4, 15)  # min/max distance to the primer on the same strand
QAMPLICON_LENGTH = (70, 200)  # min/max length of the qPCR amplicon
QAMPLICON_GC = (40, 60)  # GC min/max of the qPCR amplicon
QAMPLICON_DELTAG_CUTOFF = -1  # minimum free energy (kcal/mol/K) at the lowest primer temp

# PCR parameters
PCR_MV_CONC = 100  # monovalent cations mM
PCR_DV_CONC = 2  # divalent cations mM
PCR_DNTP_CONC = 0.8  # dntp concentration mM
PCR_DNA_CONC = 15  # primer concentration nM

# multipliers for primer and qpcr probe penalties
PRIMER_TM_PENALTY = 2  # temperature penalty
PRIMER_GC_PENALTY = 0.2  # gc penalty
PRIMER_SIZE_PENALTY = 0.5  # size penalty
PRIMER_MAX_BASE_PENALTY = 8  # max base penalty for a primer
PRIMER_3_PENALTY = (32, 16, 8, 4, 2)  # penalties for 3' mismatches
PRIMER_PERMUTATION_PENALTY = 0.1  # penalty for the number of permutations


# DO NOT CHANGE
# nucleotide definitions
NUCS = set("atcg")
AMBIG_NUCS = {
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
