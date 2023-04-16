**var**iable **V**irus**AMP**licons (varVAMP) is a tool to design primers for highly diverse viruses. The input is an alignment of your viral (full-genome) sequences.

# varVAMP
[![language](https://img.shields.io/badge/python-%3E3.9-green)](https://www.python.org/)
[![License: GPL v3](https://img.shields.io/github/license/jonas-fuchs/varvamp)](https://www.gnu.org/licenses/gpl-3.0)
[![PiPy](https://img.shields.io/pypi/v/varvamp?label=pypi%20version)](https://pypi.org/project/varvamp/)
[![PiPy](https://static.pepy.tech/personalized-badge/varvamp?period=total&units=international_system&left_color=grey&right_color=green&left_text=pypi%20downloads)](https://pypi.org/project/varvamp/)
[![CONDA](https://img.shields.io/conda/v/bioconda/varvamp?label=conda%20version)](https://anaconda.org/bioconda/varvamp)
[![CONDA](https://img.shields.io/conda/dn/bioconda/varvamp?label=conda%20downloads)](https://anaconda.org/bioconda/varvamp)
<img src= https://anaconda.org/conda-forge/r-clv/badges/platforms.svg /> 
[![DOI](https://zenodo.org/badge/606756158.svg)](https://zenodo.org/badge/latestdoi/606756158)

For a lot of virus genera it is difficult to design pan-specific primers. varVAMP solves this by introducing ambiguous characters into primers and minimizes mismatches at the 3' end. Primers might not work for some sequences of your input alignment but should recognize the large majority.

**varVAMP comes in three different flavors:**

<img src="https://github.com/jonas-fuchs/varVAMP/blob/master/docs/varvamp.png" alt="varVAMP logo" />

**SANGER**: varVAMP searches for the very best primers and reports back non-overlapping amplicons which can be used for PCR-based screening approaches.

**TILED**: varVAMP uses a graph based approach to design overlapping amplicons that tile the entire viral genome. This designs amplicons that are suitable for Oxford Nanopore or Illumina based full-genome sequencing.

**QPCR**: varVAMP searches for small amplicons with an optimized internal probe (TaqMan). It minimizes temperature differences between the primers and checks for amplicon secondary structures.

This program is currently being developed and in an alpha state. You are welcome to use this software. If you successfully design primers, drop me a mail. It might be possible to collaborate! Ideas and suggestions are highly welcome.

# Documentation

* [Installation](https://github.com/jonas-fuchs/varVAMP/blob/master/docs/installation.md)
* [Preparing the data](https://github.com/jonas-fuchs/varVAMP/blob/master/docs/preparing_the_data.md)
* [Usage](https://github.com/jonas-fuchs/varVAMP/blob/master/docs/usage.md)
* [Output](https://github.com/jonas-fuchs/varVAMP/blob/master/docs/output.md)
* [How it works](https://github.com/jonas-fuchs/varVAMP/blob/master/docs/how_varvamp_works.md)
* [FAQ](https://github.com/jonas-fuchs/varVAMP/blob/master/docs/FAQ.md)

---

**Important disclaimer:**
*For the primer design, varVAMP uses [primer3](https://pypi.org/project/primer3-py/) to check if digested kmers of a sequence are potential primers. Some of the functions for this were adapted from [primalscheme](www.github.com/aresti/primalscheme) and I do not claim credit.*

*The remaing code is under the GPLv3 licence. The code is WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.*
