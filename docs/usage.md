## Usage


In the most simple case, varVAMP will take an alignment and an output directory to store the results.

**Minimal usage:**

```shell
varvamp <mode> <alignment> <output dir>
```

**Full usage:**
```shell
usage: 	varvamp <mode> --help
	varvamp <mode> [mode optional arguments] <alignment> <output dir>
```

```
positional arguments:
  input                 alignment file and dir to write results

optional arguments:
  -h, --help            show this help message and exit
  --verbose, --no-verbose
                        show varvamp console output (default: True)
  -v, --version         show program's version number and exit

varvamp mode:
  {sanger,tiled,qpcr}
    sanger              design primers for sanger sequencing
    tiled               design primers for whole genome sequencing
    qpcr                design qPCR primers

```
**sanger** mode:
```shell
usage: varvamp sanger [optional arguments] <alignment> <output dir>
```
```
optional arguments:
  -h, --help            	show this help message and exit
	-t, --threshold 	threshold for consensus nucleotides
  -a , --n-ambig        	max number of ambiguous characters in a primer
  -ol 1000, --opt-length 1000   optimal length of the amplicons
  -ml 1500, --max-length 1500   max length of the amplicons
  -n inf, --report-n inf	report the top n best hits
```
**tiled** mode:
```shell
usage: varvamp tiled [optional arguments] <alignment> <output dir>
```
```
optional arguments:
  -h, --help            	show this help message and exit
  -t, --threshold 	threshold for consensus nucleotides
  -a , --n-ambig        	max number of ambiguous characters in a primer
  -ol 1000, --opt-length 1000	optimal length of the amplicons
  -ml 1500, --max-length 1500	max length of the amplicons
  -o 100, --overlap 100		min overlap of the amplicons
```
**qpcr** mode:
```shell
usage: varvamp qpcr [optional arguments] <alignment> <output dir>
```
```
optional arguments:
  -h, --help            	show this help message and exit
	-t, --threshold 	threshold for consensus nucleotides
  -a , --n-ambig        	max number of ambiguous characters in a primer
  -pa , --pn-ambig   		 max number of ambiguous characters in a probe
  -n 50, --test-n 50    	test the top n qPCR amplicons for secondary structures at the minimal primer temperature
  -d -1, --deltaG -1    minimum free energy (kcal/mol/K) cutoff at the lowest primer melting temp


```

## Further customization (advanced)

Although we believe that this will be in the most cases not necessary, you can customize all settings for varVAMP that are not specified via commands in the `config.py`. Here are also all default parameters stored if no optional arguments are given for your data. [To fully customize varVAMP, install it directly from this GitHub repository](./installation.md)

Go to the configs location:
```shell
cd varVAMP/varvamp/scripts/
```
And open the `config.py` with a text editor, e.g.:
```shell
gedit config.py
```
Here you can adjust various settings including primer parameters and penalties.

```python
# CAN BE CHANGED, DO NOT DELETE
# basic primer parameters
PRIMER_TMP = (57, 63, 60)  # melting temperatur (min, max, opt)
PRIMER_GC_RANGE = (35, 65, 50)  # gc (min, max, opt)
PRIMER_SIZES = (18, 24, 21)  # size (min, max, opt)
PRIMER_MAX_POLYX = 3  # max number of polyx repeats
PRIMER_MAX_DINUC_REPEATS = 3  # max number of dinucleotide repeats
PRIMER_HAIRPIN = 47  # max melting temp for secondary structure
PRIMER_GC_END = (0, 4)  # min/max GCs in the last 5 bases of the 3' end
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

```
To apply these new settings just repeat the installation procedure in the varVAMP dir:
```shell
pip install .
```
If you did everything right, varVAMP's config check passes. Otherwise it will produce an error. If that happens you can simply perform a git pull or adjust the settings that produced a warning. Please use the [GitHub issues](https://github.com/jonas-fuchs/varVAMP/issues) to report any problems and bugs.

#### [Previous: Data preparation](./preparing_the_data.md)&emsp;&emsp;[Next: Output](./output.md)
