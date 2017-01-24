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
  -degen {0,2,3,4}      Degeneracy of coding SNPs to extract (must run with
                        -mode snp
  -folded               If specified will output minor allele spectrum
```

#### Options

 * ```-vcf``` specifies a VCF to extract site frequencies from
 * ```-chr``` used to specify a specific chromosome to run on, note that for this to work, the VCF file must be compressed with bgzip and indexed with tabix
 * ```-region``` used to specify a specific genomic region to extract the site frequencies from, eg) intron or intergenic, this relies on annotation information to be present in the VCF info field in the form of ```ANNO=region```
 * ```-mode``` is used to determine the variant type to extract, eith 'snp', 'indel', 'del', or 'ins'
 * ```-degen``` used to specify degeneracy of snps desired for site frequencies, relies on degeneracy information being present in the VCF info field in the form ```DEGEN=int``` eg) ```DEGEN=0``` for zerofold snps. Must be run in conjunction with ```-mode snp``` 
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

#### Sorted, counted SFS for quick plotting

```
./vcf2raw_sfs.py -vcf data/test_data_sfs.vcf.gz -mode indel -folded | sort | uniq -c
```

Which yields:

```
      2 0.05
      2 0.1
      1 0.2
      2 0.25
      1 0.3
      1 0.35
      1 0.45
```

#### folded SFS for zerofold SNPs

```
/vcf2raw_sfs.py -vcf data/test_data_sfs_snp.vcf.gz -mode snp -degen 0 -folded
```

Which returns:

```
0.05
0.1
```

### Examples - calling within python

todo
