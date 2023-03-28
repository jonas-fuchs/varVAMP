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
1. Digest the conserved regions into kmers with the min and max length of primers. This is performed on a consensus sequence that does not contain ambiguous characters but is just the majority consensus of the alignment. Therefore, primer parameters will be later calculated for the best fitting primer.
2. Evaluate if these kmers are potential primers direction independent (temp, gc, size, poly-x repeats and poly dinucleotide repeats) and direction dependent (secondary structure, gc clamp, number of GCs in the last 5 bases of the 3'end and min 3' nucleotides without a ambiguous base). Filter for kmers that satisfy all constraints and calculate their penalties (explained in the last section).
3. Find lowest scoring primer. varVAMP sorts the primers by their score and always takes the best scoring if the primer positions have not been covered by a better primer. This greatly reduces the complexity of the later amplicon search while only retaining the best scoring primer of a set of overlapping primers.

### Amplicon search

#### Amplicon-tiling
To search for the best scoring amplicon, varVAMP uses a graph based approach.
1. Create all possible amplicons with the given length constraints and ensure that primer pairs are not forming dimers.
2. Create a graph containing all amplicons and their potential neighboring amplicons. To design a good scheme, the next primer has to lie within the second half of the current primer and satisfy the overlap constraint. The cost to go to a neighboring amplicon is determined by the amplicon score.
3. Use the dijkstra algorithm to find the path with the lowest costs from a given start node.
4. Determine potential stop nodes as those amplicons that have the furthest stop in the alignment.
5. Determine shortest paths between the start and all potential stop nodes. Get the lowest scoring.
6. Repeat steps 3-5 for each start node until the best coverage is reached. This then is the best scoring amplicon scheme!
7. Lastly, the best scoring scheme is evaluated for primer dimers in their respective pools. If a primer dimer pair is found, varVAMP evaluates for each primer their overlapping previously not considered primers (primer search step 2) and again minimizes the score. The scheme and all primers are updated. If no alternative primers can be found, varVAMP issues a warning and reports the unsolvable primer dimers.

#### Sanger sequencing
varVAMP sorts all amplicons by their score and takes the non-overlapping amplicon with the lowest score! As varVAMP gives a size penalty to amplicons, varVAMP automatically finds amplicons with low primer scores close to your optimal length (if possible).

#### qPCR
coming soon

### Penalty calculation

```python
PRIMER_MAX_BASE_PENALTY
```
Each primer is scored for its deviation from the optimal size, gc content and temperature. Base penalty is the sum of these deviations. Primer base penalties higher than the max base penalty are excluded.

```python
PRIMER_3_PENALTY
```
Each position in the primer is scored for mismatches in all sequences. If a 3' penalty is given the first in the tuple is multiplied with the freq mismatch at the very 3' end. The next is multiplied with the -1 freq and so on. Increase penalty if you want to shift amplicons towards best 3' matching. set to 0 if you do not care about 3' mismatches.

```python3
PRIMER_PERMUTATION_PENALTY
```
The number permutations of a primer is multiplied by the penalty. For example 24 permutations and a penalty of 0.1 will yield a penalty of 2.4. Set to 0 if you do not care about the number of permutations.

All scores of a primer are summed up and yield a final score. The score for each amplicon is then the score of its LEFT + RIGHT primers multiplied by the fold increase of the amplicon length compared to the optional length. This insures that in the final scheme not only large amplicons are used.

#### [Previous: Output](./output.md)&emsp;&emsp;[Next: FAQ](./FAQ.md)
