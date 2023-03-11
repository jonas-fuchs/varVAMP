## Requirements
varVAMP runs on UNIX Systems and MacOSX. Python3 >=3.9 has to be installed.


## Installation of varvamp with pip

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

## Installation via requirements.txt for developement porposes

```shell
git clone https://github.com/jonas-fuchs/varVAMP
cd varVAMP
pip install -r requirements.txt
python3 varvamp -v
```

## Update
```shell
cd varVAMP
git pull
python3 install .
```

#### [Next: Preparing your data](./preparing_the_data.md)
