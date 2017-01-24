# sfs_utils
Henry Juho Barton

## Introduction

A collection of scripts to extract the site frequency spectrum from various file formats.
Included scripts:
  * vcf2raw_sfs.py - extracts SFS from VCF for specified regions, variant types and chromosomes
  
## vcf2raw_sfs.py

This script uses the pysam python module (<https://github.com/pysam-developers/pysam>), to write site frequencies to stdout for all sites that meet the specified criteria (see usage).

### Usage 

todo

### Examples - commandline

Below are some examples for different extracted site frequencies

#### SFS for all deletions genome wide

```bash
./vcf2raw_sfs.py -vcf data/test_data_sfs.vcf.gz -mode del
```

Which yields:

```bash
0.05
0.35
0.3
0.05
0.1
0.75
0.2
```

#### SFS for intergenic insertions on chromosome 10

```bash
./vcf2raw_sfs.py -vcf data/test_data_sfs.vcf.gz -mode ins -region intergenic -chr chr10
```

Which yields:

```
0.1
0.25
0.45
```

### Examples - calling within python

todo
