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
  -mute_type {WW,SS,SW,WS}
                        Mutation type, use only with mode -snp
  -folded               If specified will output minor allele spectrum
  -multi_allelic        If specified will not restrict output to biallelic
                        sites
  -bed                  If specified will output allele frequncies in bed
				        format,each row specifying chromosome start end
				        allele_frequency

```

#### Options

 * ```-vcf``` specifies a VCF to extract site frequencies from
 * ```-chr``` used to specify a specific chromosome to run on, note that for this to work, the VCF file must be compressed with bgzip and indexed with tabix
 * ```-region``` used to specify a specific genomic region to extract the site frequencies from, eg) intron or intergenic, this relies on annotation information to be present in the VCF info field in the form of ```ANNO=region```
 * ```-mode``` is used to determine the variant type to extract, eith 'snp', 'indel', 'del', or 'ins'
 * ```-degen``` used to specify degeneracy of snps desired for site frequencies, relies on degeneracy information being present in the VCF info field in the form ```DEGEN=int``` eg) ```DEGEN=0``` for zerofold snps. Must be run in conjunction with ```-mode snp``` 
 * ```-mute_type``` used to specify what mutation type is desired S<->S = SS, W<->W = WW, W->S = WS and S->W = SW. Can only be used with ```-mode snp``` and SW and WS will only output if ```-folded``` is not specified 
 * ```-folded``` if present will output the folded (minor allele) site frequencies (cannot be used in conjunction with -mode del or mode ins)
 * ```-multi_allelic``` if specified then will not restrict sites to biallelic (default)
 * ```-bed``` if specified will output allele frequncies in bed format, each row from left to right specifies the chromosome start end allele_frequency


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

#### folded SFS for S<->S SNPs

```
./vcf2raw_sfs.py -vcf data/test_data_sfs_snp.vcf.gz -mode snp -folded -mute_type SS
```

Output:

```
0.1
0.1
```

#### Unfolded SFS for indels output to bed format

```
./vcf2raw_sfs.py -vcf data/test_data_sfs.vcf.gz -mode indel -folded -bed
```

Output:

```
chr10   3638302 3638304 0.05
chr10   3638545 3638546 0.35
chr10   3639132 3639135 0.3
chr10   3641277 3641280 0.05
chr10   3641557 3641558 0.1
chr10   3641563 3641565 0.1
chr10   3641928 3641930 0.25
chr10   3642608 3642609 0.25
chr10   3642947 3642951 0.2
chr10   3643568 3643569 0.45
```

### Examples - calling within python

An example of calling the script within python making use of the subproccess module, allowing the SFS to be read in and stored in a list.

```python
>>> import subprocess
>>> del_sfs_cmd = './vcf2raw_sfs.py -vcf data/test_data_sfs.vcf.gz -region intergenic -mode del' 
>>> del_sfs = subprocess.Popen(del_sfs_cmd, shell=True, stdout=subprocess.PIPE).communicate()[0].split('\n')[:-1]
>>> del_sfs
['0.3', '0.05', '0.1', '0.75', '0.2']
```
