from setuptools import setup, find_packages
from varvamp import __version__, _program

setup(
    name='varvamp',
    version=__version__,
    python_requires=">=3.9",
    packages = find_packages(),
    install_requires=[
        "biopython>=1.79",
        "matplotlib>=3.5.1",
        "primer3-py>=1.1.0"
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
    """.format(program = _program),
    include_package_data=True,
    keywords=[],
    zip_safe=False
)
