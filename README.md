# varVAMP

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**var**iable **V**irus**AMP**licons is a tool to design amplicons for highly diverse viruses (like Hepatits E or Poliovirus).

**varVAMP has three different modes:**

* Amplicon sequencing
* Sanger sequencing (not implemented yet)
* qPCR (not implemented yet)


# Documentation

* [Installation](docs/installation.md)
* [Usage](docs/usage.md)
* [Output description](docs/output.md)
* [How it works](docs/how_varvamp_works.md)

---

**Important disclaimer:**
*For the primer design, varVAMP uses [primer3](https://pypi.org/project/primer3-py/) to check if digested kmers of a sequence are potential primers. The functions for this were adapted from [primalscheme](www.github.com/aresti/primalscheme) and I do not claim credit. Primalscheme is a super awsome tool. Go check it out!*
