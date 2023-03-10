from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
import glob
import os
import pkg_resources

from varvamp import __version__, _program

setup(name='varvamp',
      version=__version__,
      packages=find_packages(),
      scripts=["varvamp/scripts/alignment.py",
                "varvamp/scripts/config.py",
                "varvamp/scripts/consensus.py",
                "varvamp/scripts/conserved.py",
                "varvamp/scripts/primers.py",
                "varvamp/scripts/reporting.py",
                "varvamp/scripts/scheme.py"
      ],
      install_requires=[
            "biopython>=1.79",
            "matplotlib>=3.5.1",
            "primer3-py>=1.1.0"
        ],
      description='varvamp',
      url='https://github.com/jonas-fuchs/varVAMP',
      author='Dr. Jonas Fuchs',
      author_email='jonas.fuchs@uniklinik-freiburg.de',
      entry_points="""
      [console_scripts]
      {program} = varvamp.command:main
      """.format(program = _program),
      include_package_data=True,
      keywords=[],
      zip_safe=False)
