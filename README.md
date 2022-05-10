# Single Instruction Multiple Data with DNA (SIMD||DNA)

This software was used for next-generation sequencing dataa analyses as described in the manuscript:

**Parallel molecular computation on digital data stored in DNA**

**Boya Wang, Siyuan Wang, Cameron Chalk, Andrew Ellington, and David Soloveichik**

## Installation

simddna requires Biopython for sequence alignment and parsing FASTQ files:

```
pip install Biopython

```

Installing simddna:

```
python setup.py install

```

## Typical Workflow

Place the raw NGS read files (.fastq or .fastq.gz) into a working directory with an analysis script like example/ngsrunanalysis.py - see the example folder. Simply run the analysis script

```
python ngsrunanalysis.py

```
This produces the following .pkl files:
- Viable reads: Reads that pass the criteria for analysis
- NGS barcode calls: Viable reads organized by their decoded NGS sample barcodes
- SIMD barcode calls: Viable reads organized by both their decoded NGS barcodes and SIMD barcodes
- Interpreted bits: Viable reads organized by their barcodes and interpreted 4-bit values - only these files are needed for downstream analysis and visualization

