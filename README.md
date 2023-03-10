# varVAMP

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**var**iable **V**irus**AMP**licons is a tool to design primers for highly diverse viruses. The input is an alignment of your viral (full-genome) sequences.

**varVAMP comes in three different flavors:**

<img src="./docs/varvamp.png" alt="varVAMP logo" />

**TILED**: varVAMP uses a graph based approach to design overlapping amplicons that tile the entire viral genome. This designs amplicons that are suitable for Oxford Nanopore or Illumina based full-genome sequencing.

**SANGER** *(coming soon)*: varVAMP searches for the very best primers and reports back an amplicon which can be used for PCR-based screening approaches.

**QPCR** *(coming soon)*: varVAMP searches for small amplicons with an internal primer for the probe. It minimizes temperature differences between the primers.

Features that are planned, but have not been implemented yet:
* Dimer search for amplicon schemes (currently only primer homo-dimers and hetero-dimers within the amplicon are evaluated).
* Primer BLAST.

This program is currently being developed and in an alpha state. You are welcome to use this software. If you successfully design primers, drop me a mail. It might be possible to collaborate!

# Documentation

* [Installation](docs/installation.md)
* [Preparing the data](docs/preparing_the_data.md)
* [Usage](docs/usage.md)
* [Output](docs/output.md)
* [How it works](docs/how_varvamp_works.md)

---

**Important disclaimer:**
*For the primer design, varVAMP uses [primer3](https://pypi.org/project/primer3-py/) to check if digested kmers of a sequence are potential primers. The functions for this were adapted from [primalscheme](www.github.com/aresti/primalscheme) and I do not claim credit.*

*The remaing code is under the GPLv3 licence. The code is WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.*
