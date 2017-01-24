# sfs_utils
Henry Juho Barton

## Introduction

A collection of scripts to extract the site frequency spectrum from various file formats.
Included scripts:
  * vcf2raw_sfs.py - extracts SFS from VCF for specified regions, variant types and chromosomes
  
## vcf2raw_sfs.py

This script uses the pysam python module (<https://github.com/pysam-developers/pysam>), to write site frequencies to stdout for all sites that meet the specified criteria (see usage).

### Usage 

Usage information can be viewed as follows:

```
./vcf2raw_sfs.py -h

usage: vcf2raw_sfs.py [-h] -vcf VCF [-chr CHR] [-region REGION] -mode
                      {snp,ins,del,indel} [-folded]

optional arguments:
  -h, --help            show this help message and exit
  -vcf VCF              VCF file to extract sfs from
  -chr CHR              Chromosome to extract
  -region REGION        Genomic regions to extract, default = ALL
  -mode {snp,ins,del,indel}
                        Variant mode to run in
  -folded               If specified will output minor allele spectrum
```

#### Options

 * ```-vcf``` specifies a VCF to extract site frequencies from
 * ```-chr``` used to specify a specific chromosome to run on, note that for this to work, the VCF file must be compressed with bgzip and indexed with tabix
 * ```-region``` used to specify a specific genomic region to extract the site frequencies from, eg) intron or intergenic, this relies on annotation information to be present in the VCF info field in the form of ```ANNO=region```
 * ```-mode``` is used to determine the variant type to extract, eith 'snp', 'indel', 'del', or 'ins'
 * ```-folded``` if present will output the folded (minor allele) site frequencies (cannot be used in conjunction with -mode del or mode ins)

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
