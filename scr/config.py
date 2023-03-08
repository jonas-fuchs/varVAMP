"""
varVAMP config file
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
PRIMER_TM_PENALTY = 2
PRIMER_GC_PENALTY = 0.2
PRIMER_SIZE_PENALTY = 0.5
PRIMER_MAX_BASE_PENALTY = 8
PRIMER_3_PENALTY = (10, 10, 10)  # 3' penalty from the 3' - each must be > 0
PRIMER_PERMUTATION_PENALTY = 0.1  # penalty for the number of permutations of a primer must be >0

# amplicon settings
MIN_OVERLAP = 100
OPT_AMPLICON_LENGTH = 1000
MAX_AMPLICON_LENGTH = 2000

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
    3' matching. set to 0 if you do not care about 3' mismatches.
- permutation penalty:
    the number permutations of a primer is multiplied by the penalty. for
    example 24 permutations and a penalty of 0.1 will yield a penalty of
    2.4. set to 0 if you do not care about the number of permutations.

In the end all scores are summed up. The score for each amplicon is
then the score of its LEFT + RIGHT primers multiplied by the fold
increase of the amplicon length comapred to the optional length.
This insures that in the final scheme not only large amplicons are
used.
"""

# DO NOT CHANGE
import sys


def confirm_config():
    """
    checks the config. raises error and warnings
    if nececarry.
    """
    # check if all variables exists
    all_vars = [
        "FREQUENCY_THRESHOLD",
        "ALLOWED_N_AMB",
        "PRIMER_TMP",
        "PRIMER_GC_RANGE",
        "PRIMER_SIZES",
        "PRIMER_HAIRPIN",
        "MAX_POLYX",
        "MAX_DINUC_REPEATS",
        "MAX_DIMER_TMP",
        "MIN_3_WITHOUT_AMB",
        "MV_CONC",
        "DV_CONC",
        "DNTP_CONC",
        "DNA_CONC",
        "PRIMER_TM_PENALTY",
        "PRIMER_GC_PENALTY",
        "PRIMER_SIZE_PENALTY",
        "PRIMER_MAX_BASE_PENALTY",
        "PRIMER_3_PENALTY",
        "PRIMER_PERMUTATION_PENALTY",
        "MAX_AMPLICON_LENGTH",
        "OPT_AMPLICON_LENGTH",
        "MIN_OVERLAP"
    ]

    for var in all_vars:
        if var not in globals():
            sys.exit(f"\n\033[31m\033[1mERROR:\033[0m {var} does not exist in config\n")

    # confirm tuples
    for type, tup in [("temp", PRIMER_TMP), ("gc",PRIMER_GC_RANGE), ("size", PRIMER_SIZES)]:
        if len(tup) != 3:
            sys.exit(f"\n\033[31m\033[1mERROR:\033[0m {type} tuple has to have the form (min, max, opt)\n")
        if PRIMER_TMP[0] > PRIMER_TMP[1]:
            sys.exit(f"\n\033[31m\033[1mERROR:\033[0m min {type} should not exeed max {type}\n")
        if PRIMER_TMP[0] > PRIMER_TMP[2]:
            sys.exit(f"\n\033[31m\033[1mERROR:\033[0m min {type} should not exeed opt {type}\n")
        if PRIMER_TMP[2] > PRIMER_TMP[1]:
            sys.exit(f"\n\033[31m\033[1mERROR:\033[0m opt {type} should not exeed max {type}\n")
        if any(map(lambda var: var < 0, tup)):
            sys.exit(f"\n\033[31m\033[1mERROR:\033[0m {type} can not contain negative values\n")

    # check values that cannot be zero
    non_negative_var = [
        ("max polyx nucleotides", MAX_POLYX),
        ("max polyx nucleotides", MAX_DINUC_REPEATS),
        ("min number of 3 prime nucleotides without ambiguous nucleotides", MIN_3_WITHOUT_AMB),
        ("monovalent cation concentration", MV_CONC),
        ("divalent cation concentration", DV_CONC),
        ("dNTP concentration", DNTP_CONC),
        ("primer temperatur penalty", PRIMER_TM_PENALTY),
        ("primer gc penalty", PRIMER_GC_PENALTY),
        ("primer size penalty", PRIMER_SIZE_PENALTY),
        ("max base penalty", PRIMER_MAX_BASE_PENALTY),
        ("primer permutation penalty", PRIMER_PERMUTATION_PENALTY)
    ]
    for type, var in non_negative_var:
        if var < 0:
            sys.exit(f"\n\033[31m\033[1mERROR:\033[0m {type} can not be negative\n")

    if any(map(lambda var: var < 0, PRIMER_3_PENALTY)):
        sys.exit("\n\033[31m\033[1mERROR:\033[0m 3' penalties can not be zero\n")

    # specific warnings
    if PRIMER_HAIRPIN < 0:
        print("\n\033[31m\033[1mWARNING:\033[0m decreasing hairpin melting temp to negative values will influence successful primer search.")
    if MAX_DIMER_TMP < 0:
        print("\n\033[31m\033[1mWARNING:\033[0m there is no need to set max dimer melting temp below 0.")
    if PRIMER_MAX_BASE_PENALTY < 8:
        print("\n\033[31m\033[1mWARNING:\033[0m decreasing the base penalty will hardfilter more primers.")


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
