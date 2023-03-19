## How varVAMP works

### Overview
varVAMP primer design for viruses with highly variable genomes. varVAMP first preprocesses the alignment and then creates consensus sequences that can contain ambiguous characters. Then it searches for conserved regions as defined by a user defined amount of ambiguous charaters within
the min length of a primer. The conserved regions of a consensus sequence containing the most prevalent nucleotide is then digested into kmers which are considered potential primers if they pass all primer requirements. These primers are further filtered for high scoring primers for each region. Then it constructs all possible amplicons and determines which amplicons are overlapping. A graph based approach is used to find the best amplicon scheme.

---

### Alignment preprocessing
The alignment preprocessing contains three steps.
1. Convert RNA to DNA
2. Force nucleotides to lower charaters
3. Clean gaps. This is the most important part. Larger insertions in single viral sequences enlarge the length of the alignment and therefore the length of amplicons will be overestimated. varVAMP masks regions it deletes with an "N" and knows now that it is not allowed to design primers spanning the potential insertion site.

### Consensus generation
varVAMP creates consensus sequence on the basis of the threshold. If a nucleotide is present in this or a higher frequency in all aligned sequences it is considered a consensus nucleotide. If not, the frequencies of the highest scoring nucleotides are added until the threshold is reached and the appropriate ambiguous nucleotide character is used. Importantly, varVAMP is aware of ambiguous nucleotides in the alignment and handels them by de-multiplexing the character into real nucleotides and adding its portion to the nucleotide counts at the alignment position.

### Conserved region search
varVAMP searches for conserved regions as defined by a user-defined number of ambiguous bases that is allowed within the minimal length of a primer. The algorithm opens windows over the ambiguous consensus sequence and determines if a window satisfies these constraints.

### Primer search
varVAMP uses primer3-py to search for potential primers. Some of the evaluation if primers match certain criteria was adapted from [primalscheme](www.github.com/aresti/primalscheme). The primer search contains multiple steps:
1. Digest the conserved regions into kmers with the min and max length of primers. This is performed on a consensus sequence that does not contain ambiguous characters but is just the majority consensus of the alignment. Therefore, primer parameters (gc, temp, thermodynamic model) will be calculated for the best fitting primer.
2. Evaluate if these kmers are potential primers direction independent (temp, gc, size etc) and direction dependent (secondary structure, gc clamp etc). Hardfilter for kmers that satisfy all constraints and calculate their penalties.
3. Find lowest scoring primer. varVAMP sorts the primers by their score and always takes the best scoring if the primer positions have not been covered by a better primer. This greatly reduces the complexity while only retaining the best scoring primers.

### Amplicon search

#### Amplicon-tiling
To search for the best scoring amplicon, varVAMP uses a graph based approach.
1. Create all possible amplicons with the given length constraints.
2. Create a graph for amplicons that satisfy the overlap constraints. The cost to go to the next node is determined by the amplicon score.
3. Use the dijkstra algorithm to find the shortest paths from a given start node.
4. Determine potential stop nodes as those amplicons that have the furthest stop in the alignment.
5. Determine shortest paths between the start and all potential stop nodes. Get the lowest scoring.
6. Repeat this process for each start node until the best coverage is reached.

This then is the best scoring amplicon scheme!

#### Sanger sequencing
coming soon

#### qPCR
coming soon

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

#### [Previous: Output](./output.md)&emsp;&emsp;[Next: FAQ](./FAQ.md)
