## UseGalaxy implementation (coming soon)

If you are unfamiliar with command-line programs, do not fear. varVAMP is also available at the community-driven web-based analysis platform [UseGalaxy](https://usegalaxy.eu).

## Requirements
varVAMP runs on UNIX Systems, MacOSX and Windows with Python3 >=3.9 installed.

## Installation

### PyPI:

```shell
pip install varvamp
```

### CONDA:

```shell
conda install -c bioconda varvamp
```

### DOCKER:

```shell
docker pull quay.io/biocontainers/varvamp:<current_tag>
```

That was already it. To check if it worked:

```shell
varvamp -v
```
You should see the current varVAMP version.

### BLAST module

If you want to use varVAMPs blast module to predict off-targets, make sure that you have installed [BLASTN](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#downloadblastdata) and created a local blast database containing sequences of potential off-targets.


## Installation for development

### - via pip (recommended)

```shell
git clone https://github.com/jonas-fuchs/varVAMP
cd varVAMP
pip install .
varvamp -v
```

### - via requirements.txt

```shell
git clone https://github.com/jonas-fuchs/varVAMP
cd varVAMP
pip install -r requirements.txt
python3 varvamp -v
```


#### [Next: Preparing your data](./preparing_the_data.md)
