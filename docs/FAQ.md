## FAQ

1. **How do I set settings for varVAMP?**
Start with setting your optimal length of your amplicon, the max length that you can tolerate and the overlap that you want to achieve. If you want to allow ambiguous characters, set them to 2 as a start. Set the threshold to the mean sequence identity of your alignment. Run varvamp and then optimize the output.

2. **How do I optimize the output?**
It all depends on how many conserved regions varVAMP is able to find! There are two main parameters that influence this. The number of ambiguous bases allowed within a primer and the threshold for consensus nucleotides. Settings the threshold higher or the number of ambiguous bases lower will result in less conserved regions. If you have set the parameters below and get a decent output, increase the threshold until the output gets worse. This will increase the specificity of your primers. Likewise, if you do not have a good output, consider increasing the number of ambiguous bases before you lower the threshold. The console output varVAMP will also give you some suggestions.

#### [Previous: How it works](./how_varvamp_works.md)
