from setuptools import setup, find_packages
from varvamp import __version__, _program

# read the contents of your README file
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name='varvamp',
    long_description=long_description,
    long_description_content_type='text/markdown',
    version=__version__,
    python_requires=">=3.9",
    license_files=('licence.txt'),
    packages=find_packages(),
    install_requires=[
        "biopython>=1.79",
        "matplotlib>=3.5.1",
        "primer3-py>=1.1.0",
        "pandas>=1.4.4",
        "numpy>=1.23.3"
    ],
    description='varvamp',
    url='https://github.com/jonas-fuchs/varVAMP',
    author='Dr. Jonas Fuchs',
    author_email='jonas.fuchs@uniklinik-freiburg.de',
    classifiers=[
        "Programming Language :: Python :: 3.9",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)"
    ],
    entry_points="""
    [console_scripts]
    {program} = varvamp.command:main
    """.format(program=_program),
    include_package_data=True,
    keywords=[],
    zip_safe=False
)
