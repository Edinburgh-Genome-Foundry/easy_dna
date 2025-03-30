from setuptools import setup, find_packages

version = {}
with open("easy_dna/version.py") as fp:
    exec(fp.read(), version)

setup(
    name="easy_dna",
    version=version["__version__"],
    author="Zulko",
    description="Methods for DNA sequence reading, writing and editing.",
    long_description=open("pypi-readme.rst").read(),
    license="MIT",
    url="https://github.com/Edinburgh-Genome-Foundry/easy_dna",
    keywords="DNA sequence Genbank record editing",
    packages=find_packages(exclude="docs"),
    include_package_data=True,
    install_requires=[
        "numpy",
        "Biopython",
        "snapgene_reader",
        "flametree",
        "pandas",
        "crazydoc",
    ],
)
