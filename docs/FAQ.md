## FAQ

1. **How do I set settings for varVAMP?**

Settings differ for each mode.

For SANGER mode: Start with setting your optimal length of your amplicon, the max length that you can tolerate and if you want to allow ambiguous characters, set them to 4 as a start.
For TILED mode:  Also set the overlap that you want to achieve.
For QPCR mode: Also set the number of ambiguous character for the probe lower than the value for the primers. You can also increase the number of amplicons tested for their deltaG, but be aware that this might increase running time.

If you are not satisfied with the results try to optimize and run varVAMP again. 

2. **How do I optimize the output?**

It all depends on how many conserved regions varVAMP is able to find! There are two main parameters that influence this. The number of ambiguous bases allowed within a primer and the threshold for consensus nucleotides. Setting the threshold higher or the number of ambiguous bases lower will result in less conserved regions. If you have set the parameters below and get a decent output, increase the threshold until the output gets worse. This will increase the specificity of your primers. Likewise, if you do not have a good output, consider increasing the number of ambiguous bases before you lower the threshold. The console output varVAMP will also give you some suggestions.

3. **varVAMP TILED reported primer dimers. What now?**

In your case varVAMP could not find suitable replacement primers in the TILED mode. You can either rerun varVAMP and try different settings or you can perform a third pool that contains a amplicon that has one of the conflicting dimers. Notably, varVAMP also reports the dimer melting temperature. If it is still reasonable low, using a hot start polymerase might still lead to successful PCR amplification.

4. **I have multiple hits after SANGER/QPCR mode. Which should I use?**

varVAMP sorts all amplicons and qpcr designs by score and always assigns the lowest number to the best one of non-overlapping amplicons/qpcr schemes. If you are not interested in a specific gene region, amplicon_0 or qpcr_scheme_0  are your best candidates!

5. **What is deltaG reported for the QPCR mode?**

DeltaG is the [minimal free energy](https://en.wikipedia.org/wiki/Gibbs_free_energy) of your amplicon. The more positive the value is, the more unlikely it is that the qPCR amplicon forms secondary structures during PCR. Check the structure with [mfold](http://www.unafold.org/Dinamelt/applications/two-state-melting-folding.php). DeltaGs calculated here slightly differ to [seqfold](https://github.com/Lattice-Automation/seqfold) used by varVAMP but are in a similar range. It should give you a well enough estimate of the amplicon folding at the lowest PCR temperature.

6. **Does varVAMP evaluate the off-target specificity of the designed oligos?**

This you have to do [manually](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) at the moment. A BLAST feature is planned for the future.

7. **How fast is varVAMP?**

varVAMP is pretty fast given the complexity of the problem. Running time is depended on the alignment length, number of sequences and the running mode. While the TILED is rather slow, qPCR and SANGER can be faster. An alignment with a few hundred sequences and with a genome size of 10 kb will likely run in under a minute for the TILED mode. For large e.g. DNA viruses (200 kb) it takes considerably longer, but should still finish in minutes. Running time optimizations are planned.


8. **Can I contribute?**

Yes, please. Give feedback to code or settings or also if you have successfully used this to design primers that make your life easier in the lab!



#### [Previous: How it works](./how_varvamp_works.md)
