## varVAMP output

varVAMP produces multiple main output files:


| Mode | Output | Description |
| --- | --- | --- |
| ALL | ambiguous_consensus.fasta | The consensus sequence containing ambiguous nucleotides. |
| ALL | amplicon_plot.pdf | A nice overview for your final amplicon design. |
| ALL| amplicons.bed | A bed file showing the amplicon location compared to the consensus sequence. |
| ALL| per_base_mismatches.pdf | Barplot of the percent mismatches at each nucleotide position of the primer. |
| TILED/SANGER | primer_to_amplicon_assignments.tabular | Simple tab separated file, which primers belong together. Useful for bioinformatic workflows that include primer trimming |
| ALL | primers.bed | A bed file with the primer locations. Includes the primer score. The lower, the better. |
| ALL | varvamp_log.txt | Log file. |
| TILED | unsolvable_primer_dimers.tsv | Only produced if there are primer dimers without replacements. Tells which primers form dimers and at which temperature.
| TILED/SANGER | primer.tsv | A tab separated file with important parameters for the primers including the sequence with ambiguous nucleotides (already in the right strand) and the gc and temperature of the best fitting primer as well as for the mean for all permutations of the primer. |
| QPCR | qpcr_design.tsv | A tab separated file with important parameters for the qPCR amplicon including the deltaG. |
| QPCR | qpcr_primers.tsv | A tab separated file with important parameters for the primers  and probes including the sequence with ambiguous nucleotides (already in the right strand) and the gc and temperature of the best fitting primer and probe as well as for the mean for all permutations. |


It also produces some secondary output files [*data/*]:

| Mode | Output | Description |
| --- | --- | --- |
| ALL | alignment_cleaned | The preprocessed alignment. |
| ALL | majority_consensus.fasta | Consensus sequence that does not have ambiguous characters but instead has the most prevalent nucleotide at each position. |
| ALL | conserved_primer_regions.bed | A bed file showing where the conserved regions lie on the ambiguous consensus. |
| TILED/SANGER | all_primers.bed | A bed file with all high scoring primers that varVAMP found. |
| qPCR | conserved_probe_regions.bed | A bed file showing where the conserved regions lie on the ambiguous consensus. |

#### [Previous: Usage](./usage.md)&emsp;&emsp;[Next: How varVAMP works](./how_varvamp_works.md)
