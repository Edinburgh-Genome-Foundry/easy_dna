# This will try to import setuptools. If not here, it will reach for the embedded
# ez_setup (or the ez_setup package). If none, it fails with a message
try:
    from setuptools import setup, find_packages
except ImportError:
    try:
        import ez_setup

        ez_setup.use_setuptools()
    except ImportError:
        raise ImportError(
            "EasyDNA could not be installed, probably because"
            " neither setuptools nor ez_setup are installed on"
            "this computer. \nInstall ez_setup "
            "([sudo] pip install ez_setup) and try again."
        )

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
