## FAQ

1. **How do I set settings for varVAMP?**

Settings differ for each mode. However, varVAMP comes with some pre-defined values and can optimize the threshold or number of ambiguous bases by itself if no arguments are supplied.

If you are not satisfied with the results try to optimize and run varVAMP again.

2. **How do I optimize the output?**

It all depends on how many primer regions varVAMP is able to find! There are two main parameters that influence this. The number of ambiguous bases allowed within a primer and the threshold for consensus nucleotides. Setting the threshold higher or the number of ambiguous bases lower will result in less primer regions. If you have set the parameters below and get a decent output, increase the threshold until the output gets worse. This will increase the specificity of your primers. Likewise, if you do not have a good output, consider increasing the number of ambiguous bases before you lower the threshold.

The most easy way is to set the number of ambiguous characters you can tolerate in your primers and let varVAMP optimize the threshold. Have a look in the log file, to get detailed information on this. Afterwards you can try to optimize by incrementing the threshold 0.1.

3. **varVAMP TILED reported primer dimers. What now?**

In your case varVAMP could not find suitable replacement primers in the TILED mode. You can either rerun varVAMP and try different settings or you can perform a third pool that contains a amplicon that has one of the conflicting dimers. Notably, varVAMP also reports the dimer melting temperature. If it is still reasonable low, using a hot start polymerase might still lead to successful PCR amplification.

4. **I have multiple hits after SINGLE/QPCR mode. Which should I use?**

varVAMP sorts all amplicons and qpcr designs by their penalty and always assigns the lowest number to the one with the lowest penalty of the non-overlapping amplicons/qpcr schemes. If you are not interested in a specific gene region, amplicon_0 or qpcr_scheme_0  are your best candidates!

5. **What is deltaG reported for the QPCR mode?**

DeltaG is the [minimal free energy](https://en.wikipedia.org/wiki/Gibbs_free_energy) of your amplicon. The more positive the value is, the more unlikely it is that the qPCR amplicon forms secondary structures during PCR. Check the structure with [`mfold`](http://www.unafold.org/Dinamelt/applications/two-state-melting-folding.php). DeltaGs calculated here slightly differ to [`seqfold`](https://github.com/Lattice-Automation/seqfold) used by varVAMP but are in a similar range. It should give you a well enough estimate of the amplicon folding at the lowest PCR temperature.

6. **Does varVAMP evaluate the off-target specificity of the designed oligos?**

Yes, you can run primers against a local [BLAST](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#downloadblastdata) db via the `-db` argument.

7. **Why does the % primer regions differ from the coverage estimated with the automatic parameter selection?**

The coverage is estimated on an alignment that still has gaps. If there are a lot of large gaps the consensus sequence will be shorter and the % of the alignment that is considered for the alignment will be larger. The parameter selection is therefore only a rough estimate to provide settings at which varVAMP will perform reasonably well. However, it is likely that the output can be tweaked to find even better primers.

8. **Can I customize all parameters (e.g. primer temperature)**

[Yes you can!](./usage.md) Just edit the config file as explained. If you want to use your pip installation, you have to make it [editable](https://setuptools.pypa.io/en/latest/userguide/development_mode.html).

9. **How fast is varVAMP?**

varVAMP is pretty fast given the complexity of the problem. Running time is depended on the alignment length, number of sequences and the running mode. While the TILED is rather slow, qPCR and SANGER can be faster. An alignment with a few hundred sequences and with a genome size of 10 kb will likely run in under a minute for the TILED mode. For large e.g. DNA viruses (200 kb) it takes considerably longer, but should still finish in minutes. Running time optimizations are planned.

10. **Can I contribute?**

Yes, please. Give feedback to code or settings or also if you have successfully used this to design primers that make your life easier in the lab!



#### [Previous: How it works](./how_varvamp_works.md)
