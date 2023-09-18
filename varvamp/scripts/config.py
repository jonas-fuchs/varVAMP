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
PRIMER_HAIRPIN = 47  # max melting temp for secondary structure
PRIMER_GC_END = (1, 3)  # min/max GCs in the last 5 bases of the 3' end
PRIMER_MIN_3_WITHOUT_AMB = 3  # min len of 3' without ambiguous charaters
PRIMER_MAX_DIMER_TMP = 47  # max melting temp for dimers (homo- or heterodimers)

# QPCR parameters
# basic probe parameters
QPROBE_TMP = (64, 70, 67)  # mean 7°C higher than the primer temp
QPROBE_SIZES = (20, 30, 25)
QPROBE_GC_RANGE = (40, 80, 60)
QPROBE_GC_END = (0, 4)
# constraints for amplicon design
QPRIMER_DIFF = 2  # maximal temperature diff of qPCR primers
QPROBE_TEMP_DIFF = (5, 10)  # min/max temp diff between probe and primers
QPROBE_DISTANCE = (4, 15)  # min/max distance to the primer on the same strand
END_OVERLAP = 5  # maximum allowed nt overlap between the ends of two primers
QAMPLICON_LENGTH = (70, 200)  # min/max length of the qPCR amplicon
QAMPLICON_GC = (40, 60)  # GC min/max of the qPCR amplicon
QAMPLICON_DEL_CUTOFF = 4  # consider regions of the alignment for deltaG calculation if they have smaller deletions than cutoff

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

# BLAST parameters (ref: PrimerBLAST (YE, Jian, et al. Primer-BLAST: a tool to design
# target-specific primers for polymerase chain reaction. BMC bioinformatics, 2012, 13.
# Jg., S. 1-11.)
BLAST_SETTINGS = {  # blast settings for query search
    "outfmt": "6 qseqid sseqid qlen length mismatch gapopen sstart send sstrand",  # do NOT change
    "evalue": 5000,
    "reward": 1,
    "penalty": -1,
    "gapopen": 2,
    "gapextend": 1
}
BLAST_MAX_DIFF = 0.5  # min percent match between primer and BLAST hit (coverage and/or mismatches)
BLAST_SIZE_MULTI = 2  # multiplier for the max_amp size of off targets (in relation to max amp size)
BLAST_PENALTY = 50  # amplicon score increase -> considered only if no other possibilities

# nucleotide definitions, do NOT change
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
