## Requirements
varVAMP runs on UNIX Systems, MacOSX and Windows with Python3 >=3.9 installed.


## Installation of varvamp with pip

```shell
git clone https://github.com/jonas-fuchs/varVAMP
cd varVAMP
pip install .
```
That was already it. To check if it worked:
```shell
varvamp -v
```
You should see the current varVAMP version.

## Installation via requirements.txt for development

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
pip install .
```

#### [Next: Preparing your data](./preparing_the_data.md)
