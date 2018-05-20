[![DOI](https://zenodo.org/badge/134107589.svg)](https://zenodo.org/badge/latestdoi/134107589)
# _Odontosyllis undecimdonta_ luciferase

This repository contains analyses required to recreate the analyses found in Schultz, Kotlobay, et al. 2018

Specifically, it contains the scripts to recreate figure 2 (the
transcript alignment figure), as well as the orthology search with other
polychaetes.

## Figure 2

[add figure 2 here]

## Requirements:

This is a reproducible pipeline assuming that all of the required
software and sequencing reads are on your file system.

Required software:
- [Python 3.x - I recommend Anaconda](https://docs.anaconda.com/anaconda/install/)
- [Snakemake](http://snakemake.readthedocs.io/en/stable/tutorial/setup.html)
- [Trimmomatic v0.35](http://www.usadellab.org/cms/?page=trimmomatic)
- [a local copy of the NCBI nr database](https://www.ncbi.nlm.nih.gov/books/NBK52640/)
- [a local operating version of blast](https://www.ncbi.nlm.nih.gov/books/NBK52640/)
- [pauvre plotting software](https://github.com/conchoecia/pauvre)
- [the Trinity RNAseq assembler docker](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-in-Docker)
- [bioawk](https://github.com/lh3/bioawk) and normal awk
- [hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml)
- [bwa](http://bio-bwa.sourceforge.net/)
- [samtools](http://samtools.sourceforge.net/)
- [minimap2](https://github.com/lh3/minimap2)

## How to run

For now you will need to edit `config.yaml` and specify the location
of the reads files from _O. undecimdonta_ sequencing project on ENA/SRA.
An automatic download feature will be added when the reads become available.

After install all of the required software above, download this
repository and execute the following commands:

```
git clone https://github.com/conchoecia/odontosyllis_undecimdonta_luciferase.git
cd odontosyllis_undecimdonta_luciferase/
snakemake --cores <desired number of cores>
```

This analysis took approximately 10 days on a linux server with 90
threads. Most of the computation time is spent assembling transcriptomes.
