**var**iable **V**irus**AMP**licons (varVAMP) is a tool to design primers for highly diverse viruses. The input is an alignment of your viral (full-genome) sequences.

# varVAMP

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

For a lot of virus genera it is difficult to design pan-specific primers. varVAMP solves this, by introducing ambiguous characters into primers and minimizes mismatches at the 3' end. Primers might not work for some sequences of your input alignment but should recognize the large majority.

**varVAMP comes in three different flavors:**

<img src="./docs/varvamp.png" alt="varVAMP logo" />

**SANGER** *(coming soon)*: varVAMP searches for the very best primers and reports back an amplicon which can be used for PCR-based screening approaches.

**TILED**: varVAMP uses a graph based approach to design overlapping amplicons that tile the entire viral genome. This designs amplicons that are suitable for Oxford Nanopore or Illumina based full-genome sequencing.

**QPCR** *(coming soon)*: varVAMP searches for small amplicons with an internal primer for the probe. It minimizes temperature differences between the primers.

This program is currently being developed and in an alpha state. You are welcome to use this software. If you successfully design primers, drop me a mail. It might be possible to collaborate!

# Documentation

* [Installation](docs/installation.md)
* [Preparing the data](docs/preparing_the_data.md)
* [Usage](docs/usage.md)
* [Output](docs/output.md)
* [How it works](docs/how_varvamp_works.md)
* [FAQ](docs/FAQ.md)

---

**Important disclaimer:**
*For the primer design, varVAMP uses [primer3](https://pypi.org/project/primer3-py/) to check if digested kmers of a sequence are potential primers. Some of the functions for this were adapted from [primalscheme](www.github.com/aresti/primalscheme) and I do not claim credit.*

*The remaing code is under the GPLv3 licence. The code is WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.*
