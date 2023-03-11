## Usage

In the most simple case varvamp will take an alignment and a output directory to store the results in. It uses as standard settings:

* optimal_length = 1000
(optimal amplicon length)
* max_length = 2000
(maximum amplicon length)
* overlap = 100
(minimum overlap length)
* threshold = 0.9
(nucleotide consensus threshold)
* allowed_ambiguous = 4
(number of allowed ambiguous characters in a primer)

These settings are quite relaxed and can produce decent results for diverse viruses (80-90 % sequence identity). However, you can likely optimize the result. To to so, type:
```bash
varvamp --help
```
Full usage:
```
usage: varvamp <alignment> <output> [options]

varvamp: variable virus amplicon design

positional arguments:
  alignment output_dir  alignment to design primers on

optional arguments:
  -h, --help            show this help message and exit
  -ol OPT_LENGTH, --opt-length OPT_LENGTH
                        optimal length of the amplicons
  -ml MAX_LENGTH, --max-length MAX_LENGTH
                        max length of the amplicons
  -o OVERLAP, --overlap OVERLAP
                        min overlap of the amplicons
  -t THRESHOLD, --threshold THRESHOLD
                        threshold for nucleotides in alignment to be considered conserved
  -a ALLOWED_AMBIGUOUS, --allowed-ambiguous ALLOWED_AMBIGUOUS
                        number of ambiguous characters that are allowed within a primer
  --console, --no-console
                        show varvamp console output (default: True)
  -v, --version         show program's version number and exit
```

## Further customization (advanced)

Although, we believe that this will be in the most cases not necessary, you can customize all settings for varVAMP that are not specified via commands in the config.py. Here are also all default parameters stored if no optional arguments are given.

Go to the configs location:
```bash
cd varVAMP/varvamp/scripts/
```
And open the config.py with an text editor, e.g.:
```bash
gedit config.py
```
Here you can adjust various settings including primer parameters and penalties.
```python

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
```
To apply these new settings just repeat the installation procedure in the varVAMP dir:
```bash
python3 install .
```
If you did everything right, varVAMPs config check passes otherwise it will produce an error. If that happens you can simply perform a git pull or adjust the settings that produced a warning.

#### [previous: installation](./installation.md) [next: output](./output.md)
