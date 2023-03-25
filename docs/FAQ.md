## FAQ

1. **How do I set settings for varVAMP?**

Start with setting your optimal length of your amplicon, the max length that you can tolerate and the overlap that you want to achieve. If you want to allow ambiguous characters, set them to 2 as a start. Set the threshold to the mean sequence identity of your alignment. Run varvamp and then optimize the output.

2. **How do I optimize the output?**

It all depends on how many conserved regions varVAMP is able to find! There are two main parameters that influence this. The number of ambiguous bases allowed within a primer and the threshold for consensus nucleotides. Setting the threshold higher or the number of ambiguous bases lower will result in less conserved regions. If you have set the parameters below and get a decent output, increase the threshold until the output gets worse. This will increase the specificity of your primers. Likewise, if you do not have a good output, consider increasing the number of ambiguous bases before you lower the threshold. The console output varVAMP will also give you some suggestions.

3. **varVAMP reported primer dimers. What now?**

In your case varVAMP could not find suitable replacement primers. You can either rerun varVAMP and try different settings or you can perform a third pool that contains a amplicon that has one of the conflicting dimers. Notably, varVAMP also reports the dimer melting temperature. If it is still reasonable low, using a hot start polymerase might still lead to successful PCR amplification.


4. **How fast is varVAMP?**

varVAMP is pretty fast given the complexity of the problem. Running time is depended on the alignment length, number of sequences and the running mode. While the TILED is rather slow, qPCR and SANGER are faster. An alignment with a few hundred sequences and with a genome size of 10 kb will likely run in under a minute for the TILED mode. For large e.g. DNA viruses (200 kb) it takes considerably longer, but should still finish in minutes. Running time optimizations are planned.




#### [Previous: How it works](./how_varvamp_works.md)
