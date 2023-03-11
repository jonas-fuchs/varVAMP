## varVAMP output

varVAMP produces multiple main output files:


| Output | Description |
| --- | --- |
| ambiguous_consensus.fasta | The consensus sequence containing ambiguous nucleotides. |
| amplicons.bed | A bed file showing the amplicon location compared to the consensus sequence. |
| primer_to_amplicon_assignments.tabular | Simple tab seperated file, which primers belong together. Useful for bioinformatic workflows that include primer trimming |
| primers.bed | A bed file with the primer locations. Included the primer score. The lower, the better |
| primer.tsv | A tab seperated file with important paramenters for the primers including the sequence. |
| varvamp_log.txt | Log file. |

It also produces some secondary output files [*data/*]:

| Output | Description |
| --- | --- |
| alignment_cleaned | The preprocessed alignment. |
| all_primers.bed | A bed file with all high scoring primers that varVAMP found. |
| conserved_regions.bed | A bed file showing where the conserved regions lie on the ambiguous consensus. |
| majority_consensus.fasta | Consensus sequence that does not have ambiguous characters but instead has the most prevalent nucleotide at each position. |

#### [Previous: Usage](./usage.md)&emsp;&emsp;[Next: How varVAMP works](./how_varvamp_works.md)
