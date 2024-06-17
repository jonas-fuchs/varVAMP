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
  {single,tiled,qpcr}
    single              design primers for single amplicons
    tiled               design primers for whole genome sequencing
    qpcr                design qPCR primers

```
**single** mode:
```shell
usage: varvamp single [optional arguments] <alignment> <output dir>
```
```
optional arguments:
  -h, --help            	show this help message and exit
  -t, --threshold 	        threshold for consensus nucleotides
  -a , --n-ambig        	max number of ambiguous characters in a primer
  -db None, --database None     location of the BLAST db
  -th 1, --threads 1            number of threads
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
  -t, --threshold 	        threshold for consensus nucleotides
  -a , --n-ambig        	max number of ambiguous characters in a primer
  -db None, --database None     location of the BLAST db
  -th 1, --threads 1	        number of threads
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
  -t, --threshold 		threshold for consensus nucleotides
  -a , --n-ambig        	max number of ambiguous characters in a primer
  -db None, --database None     location of the BLAST db
  -th 1, --threads 1   	        number of threads
  -pa , --pn-ambig   		max number of ambiguous characters in a probe
  -n 50, --test-n 50    	test the top n qPCR amplicons for secondary structures at the minimal primer temperature
  -d -3, --deltaG -3            minimum free energy (kcal/mol/K) cutoff at the lowest primer melting temp


```

## Further customization (advanced)

Although we believe that this will in most cases not be necessary, you can customize all settings for varVAMP that are not specified via command line options in the form of a **custom config file**.

### Format of a custom config file

Custom config files need to follow the format of varVAMP's [default config file](https://github.com/jonas-fuchs/varVAMP/blob/master/varvamp/scripts/default_config.py), which you may want to use as a starting point for modifications, but they do not require an `__all__` list of all existing parameters and may provide new settings for only some of the existing parameters.

If all you want to do, for example, is to change the preferred primer size to 22 nts (from 21), a look into the default config file reveals the setting:

```python
PRIMER_SIZES = (18, 24, 21)  # size (min, max, opt)
```

and your single-line custom config file could look like this (everything following a `#` on a line serves as a comment and will be ignored by varVAMP:

```python
PRIMER_SIZES = (18, 24, 22)  # size (min, max, opt); changed opt from 21 to prefer somewhat longer primers
```

### Passing a custom config file via a shell variable

Now to pass this custom config file to varVAMP and have the single new parameter definition overwrite the one in the default config file, you can prepend `VARVAMP_CONFIG=<path/to/custom_config>` to any regular varvamp command line.

Let's assume you have
- saved your custom config file under the name `custom_config.py` in the same folder as you are running your varvamp command from, and
- the basic command you want to run is `varvamp qpcr input_alignment.fasta my_results`.

Then with

Linux:
```shell
VARVAMP_CONFIG=custom_config.py varvamp qpcr input_alignment.fasta my_results
```

Windows:
```shell
set "VARVAMP_CONFIG=custom_config.py"
varvamp qpcr input_alignment.fasta my_results
```

you run the command with the parameter change applied.

If, as another example, you have several custom config files, each optimized for a specific use case, stored in a folder `/home/me/my varvamp configs`, you might want to run:

```shell
VARVAMP_CONFIG="/home/me/my varvamp configs/sars-cov-2_config.py" varvamp tiled ncov_alignment.fasta my_results
```

Note that, in this example, the quotes around the config file path are necessary to treat it as a single path despite the spaces in the folder name.

### Custom config file priority rules

Any parameter setting(s) defined in a custom config file will always overwrite the corresponding default config settings.

Any parameters not defined in a custom config file, will be taken from the default config instead.

If you screwed up a custom config file by, e.g., breaking the format or mistyping a parameter name, varVAMP will likely let you know in the form of an error message, but it's good practice to examine the run logs carefully when you use a custom config file for the first time to see if all your parameter changes took effect.

Please use the [GitHub issues](https://github.com/jonas-fuchs/varVAMP/issues) to report any problems and bugs.

#### [Previous: Data preparation](./preparing_the_data.md)&emsp;&emsp;[Next: Output](./output.md)
