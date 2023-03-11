# varVAMP

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**var**iable **V**irus**AMP**licons is a tool to design primers for highly diverse viruses.

**varVAMP comes in three different flavors:**

* Amplicon sequencing
* Sanger sequencing (coming soon)
* qPCR (coming soon)


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
