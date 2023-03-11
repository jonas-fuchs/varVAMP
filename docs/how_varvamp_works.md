## How varVAMP works

### Overview
varVAMP primer design for viruses with highly variable genomes. varVAMP first preprocesses the alignment and then creates consensus sequences that can contain ambiguous characters. Then it searches for conserved regions as defined by a user defined amount of ambiguous charaters within
the min length of a primer. The conserved regions of a consensus sequence containing the most prevalent nucleotide is then digested into kmers which are considered potential primers if they pass all primer requirements. These primers are further filtered for high scoring primers for each region. Then it constructs all possible amplicons and determines which amplicons are overlapping. A graph based approach is used to find the best amplicon scheme.

### Alignment preprocessing

### Consensus generation

### Conserved region search

### Primer search

### Amplicon search

### Penalty calculation

```python
PRIMER_MAX_BASE_PENALTY
```
Each primer is scored for its deviation from the optimums. Base penealty is the cummulative penalty of penalties for temperature, gc and size. Primer base penalties higher than the max base penalty are hardfiltered.

```python
PRIMER_3_PENALTY
```
Each position in the primer is scored for mismatches in all sequences. If a 3' penalty is given the first in the tuple is multiplied with the freq mismatch at the very 3' end. The next is multiplied with the -1 freq and so on. Increase penalty if you want to shift amplicons torwards best 3' matching. set to 0 if you do not care about 3' mismatches.

```python3
PRIMER_PERMUTATION_PENALTY
```
The number permutations of a primer is multiplied by the penalty. For example 24 permutations and a penalty of 0.1 will yield a penalty of 2.4. Set to 0 if you do not care about the number of permutations.

In the end all scores of a primer are summed up and yield a final score. The score for each amplicon is then the score of its LEFT + RIGHT primers multiplied by the fold increase of the amplicon length comapred to the optional length. This insures that in the final scheme not only large amplicons are used.

#### [Previous: Output](./output.md)
