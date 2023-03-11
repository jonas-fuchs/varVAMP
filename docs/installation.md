## Requirements
varVAMP runs on UNIX Systems and MacOSX. Python3 >=3.9 has to be installed.


## Installation with pip
Install varVAMP by cloning this repository:
```shell
git clone https://github.com/jonas-fuchs/varVAMP
cd varVAMP
python3 install .
```
That was already it. To check if it worked:
```shell
varvamp -v
```
You should see the current varVAMP version.

## Manual installation from source
Install the following dependencies:
* biopython>=1.79
* matplotlib>=3.5.1
* primer3-py>=1.1.0 

## Update
```shell
cd varVAMP
git pull
python3 install .
```

#### [Next: Preparing your data](./preparing_the_data.md)
