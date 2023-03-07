"""
varVAMP config file
"""

# CAN BE CHANGED
# alignment and consensus creation threshold
FREQUENCY_THRESHOLD = 0.91  # freq at which a nucleotide is considered conserved
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
PRIMER_TM_PENALTY = 2
PRIMER_GC_PENALTY = 0.2
PRIMER_SIZE_PENALTY = 0.5
PRIMER_MAX_BASE_PENALTY = 8
PRIMER_3_PENALTY = (10, 10, 10)  # 3' penalty from the 3' - each must be > 1
PRIMER_PERMUTATION_PENALTY = 0.1  # penalty for the number of permutations of a primer must be >0

# amplicon settings
MAX_AMPLICON_LENGTH = 2000
OPT_AMPLICON_LENGTH = 1000
MIN_OVERLAP = 100

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

"""
-------------------------------
EXPLANATION for penalty scoring
-------------------------------

- base penalty:
    each primer is scored for its deviation from the optimal temperature, gc,
    and size. primer base penalties are higher than the max base penalty are
    hardfiltered.
- primer 3' penalty:
    each position in the primer is scored for mismatches in all sequences.
    if a 3' penalty is given the first in the tuple is multiplied with the
    freq mismatch at the very 3' end. the next is multiplied with the -1 freq
    and so on. increase penalty if you want to shift amplicons torwards best
    3' matching. empty the tuple () if you do not care about 3' mismatches.
- permutation penalty:
    the number permutations of a primer is multiplied by the penalty. for
    example 24 permutations and a penalty of 0.1 will yield a penalty of
    2.4. set to 0 if you do not care about permutations.

In the end all scores are summed up. the score for each amplicon is
then the score of its LEFT + RIGHT primers.

"""
