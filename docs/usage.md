## Usage


In the most simple case varvamp will take an alignment and a output directory to store the results in.

**Minimal usage:**

```shell
varvamp <alignment> <output dir>
```

In this case varVAMP uses as standard settings:

* ```OPT_LENGTH``` = 1000 (optimal amplicon length)
* ```MAX_LENGTH``` = 2000 (maximum amplicon length)
* ```OVERLAP``` = 100 (minimum overlap length)
* ```THRESHOLD``` = 0.9 (nucleotide consensus threshold)
* ```N_AMBIG``` = 4 (number of allowed ambiguous characters in a primer)
* ```MODE``` = TILED

These settings are quite relaxed and can produce decent results for diverse viruses (80-90 % sequence identity). However, you can likely optimize the result.

**Full usage:**
```shell
usage: varvamp <alignment> <output dir> <mode> [options]
```

```
varvamp: variable virus amplicon design

positional arguments:
  input                alignment file and dir to write results
  mode                 varvamp modes: TILED, SANGER

optional arguments:
  -h, --help           show this help message and exit
  -ol , --opt-length   optimal length of the amplicons
  -ml , --max-length   max length of the amplicons
  -t , --threshold     threshold for conserved nucleotides
  -a , --n-ambig       max number of ambiguous characters in a primer
  -o , --overlap       TILED: min overlap of the amplicons
  -n , --report-n      SANGER: report the top n best hits
  --verb, --no-verb    show varvamp console output (default: True)
  -v, --version        show program's version number and exit
```

## Further customization (advanced)

Although, we believe that this will be in the most cases not necessary, you can customize all settings for varVAMP that are not specified via commands in the config.py. Here are also all default parameters stored if no optional arguments are given. [To fully customize varvamp install it directly from this github repository](./installation.md)

Go to the configs location:
```shell
cd varVAMP/varvamp/scripts/
```
And open the config.py with an text editor, e.g.:
```shell
gedit config.py
```
Here you can adjust various settings including primer parameters and penalties.

```python
# basic primer parameters
PRIMER_TMP = (57, 63, 60)  # temperatur (min, max, opt)
PRIMER_GC_RANGE = (40, 60, 50)  # gc (min, max, opt)
PRIMER_SIZES = (17, 27, 20)  # size (min, max, opt)
PRIMER_MAX_POLYX = 4  # max number of polyx repeats
PRIMER_MAX_DINUC_REPEATS = 4  # max number of dinucleotide repeats
PRIMER_HAIRPIN = 47  # max melting temp for secondary structures
PRIMER_MAX_GC_END = 3  # max GCs in the last 5 bases of the primer
PRIMER_GC_CLAMP = 1  # min number of GC nucleotides at the very 3' end
PRIMER_MIN_3_WITHOUT_AMB = 2  # min len of 3' without ambiguous charaters
PRIMER_MAX_DIMER_TMP = 47  # max melting temp for dimers (homo- or heterodimers)

# PCR parameters
PCR_MV_CONC = 50  # monovalent cations mM
PCR_DV_CONC = 2  # divalent cations mM
PCR_DNTP_CONC = 0.8  # dntp concentration mM
PCR_DNA_CONC = 50  # primer concentration nM

# multipliers for primer base penalties
PRIMER_TM_PENALTY = 2  # temperature penalty
PRIMER_GC_PENALTY = 0.2  # gc penalty
PRIMER_SIZE_PENALTY = 0.5  # size penalty
PRIMER_MAX_BASE_PENALTY = 8  # max base penalty for a primer
PRIMER_3_PENALTY = (10, 10, 10)  # penalties for 3' mismatches
PRIMER_PERMUTATION_PENALTY = 0.1  # penalty for the number of permutations
```
To apply these new settings just repeat the installation procedure in the varVAMP dir:
```shell
pip install .
```
If you did everything right, varVAMP's config check passes. Otherwise it will produce an error. If that happens you can simply perform a git pull or adjust the settings that produced a warning.

#### [Previous: Data preparation](./preparing_the_data.md)&emsp;&emsp;[Next: Output](./output.md)
