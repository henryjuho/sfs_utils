#!/usr/bin/env python

from __future__ import print_function
import argparse
import pysam
import sys


# functions
def custom_info_expression(expression, vcf_line):

    """
    looks at at custom info fields and if they meet criteria
    :param expression: str
    :param vcf_line: pysam variant
    :return: bool
    """

    if '=' in expression:

        field, filter_val = expression.split('=')
        anno_val = int(vcf_line.info[field])
        if anno_val == int(filter_val):
            return True
        else:
            return False

    elif '>' in expression:

        field, filter_val = expression.split('>')
        anno_val = int(vcf_line.info[field])
        if anno_val > int(filter_val):
            return True
        else:
            return False

    else:

        field, filter_val = expression.split('<')
        anno_val = int(vcf_line.info[field])
        if anno_val < int(filter_val):
            return True
        else:
            return False


def indel_length(vcf_line):

    """
    calculates INDEL length
    :param vcf_line: pysam variant
    :return: int
    """

    ref_seq = vcf_line.ref
    alt_seq = vcf_line.alts[0]

    length = abs(len(ref_seq) - len(alt_seq))

    return length


def get_derived_freq(vcf_line, run_mode, no_samples):

    """
    function takes a pysam vcf variant and a variant mode and returns the derived allele frequency for that variant,
    if it matches the given variant mode (eg ins, if insertion spectrum required)
    :param vcf_line: pysam variant
    :param run_mode: string
    :param no_samples: int
    :return: None or float
    """

    ref_seq = vcf_line.ref
    alt_seq = vcf_line.alts[0]
    alt_freq = round(vcf_line.info['AC'][0]/float(no_samples*2), 3)
    try:
        anc_seq = vcf_line.info['AA']
    except KeyError:  # ie. not polarised
        return None

    # set derived sequence and freq
    if alt_seq == anc_seq:
        derv_seq = ref_seq
        derv_freq = 1 - alt_freq
    elif ref_seq == anc_seq:
        derv_seq = alt_seq
        derv_freq = alt_freq
    else:
        return None

    # determine type
    if len(anc_seq) > len(derv_seq):  # deletion
        if run_mode == 'del':
            return derv_freq
        else:
            return None
    elif len(anc_seq) < len(derv_seq):  # insertion
        if run_mode == 'ins':
            return derv_freq
        else:
            return None
    else:  # snp
        return derv_freq


def get_minor_freq(vcf_line, run_mode, no_samples):

    """
    takes a pysam variant and returns the minor allele frequency
    :param vcf_line: pysam variant
    :param run_mode: string
    :param no_samples: int
    :return: float
    """

    if is_indel(vcf_line):
        variant_type = 'indel'
    else:
        variant_type = 'snp'

    if run_mode != variant_type:
        return None

    alt_allele_freq = round(vcf_line.info['AC'][0]/float(no_samples*2), 3)
    if alt_allele_freq <= 0.5:
        return alt_allele_freq
    else:
        return 1 - alt_allele_freq


def is_indel(variant):

    """
    takes a pysam variant and return whether or not it is an indel
    :param variant: pysam variant
    :return: bool
    """

    if variant.rlen > 1:
        return True

    allele_lengths = [len(allele) for allele in variant.alts]

    if len(set(allele_lengths)) > 1:
        return True

    if allele_lengths[0] > 1:
        return True
    else:
        return False


def get_out_freq(vcf_line, pol, run_mode, no_samples):

    """
    takes pysam variant, polarisation argument and variant run type
    :param vcf_line: pysam variant
    :param pol: bool
    :param run_mode: str
    :param no_samples: int
    :return: None or float
    """

    if pol is True:
        return get_derived_freq(vcf_line, run_mode, no_samples)
    else:
        return get_minor_freq(vcf_line, run_mode, no_samples)


def in_regions(vcf_line, target_regions):

    """
    takes a pysam variant and sees if it falls within a specified genomic region
    :param vcf_line: pysam variant
    :param target_regions: list or None
    :return: bool
    """

    if target_regions is None:
        return True
    else:
        try:
            var_region = vcf_line.info['ANNO']
            if var_region in target_regions:
                return True
            else:
                return False
        except KeyError:
            return False


def is_degen(vcf_line, target_degen):

    """
    takes a pysam variant and desired degeneracy and returns true or false for that variant
    :param vcf_line: pysam variant
    :param target_degen: int
    :return: bool
    """

    if target_degen is None:
        return True
    else:
        try:
            degeneracy = int(vcf_line.info['DEGEN'])
            if target_degen == degeneracy:
                return True
            else:
                return False
        except KeyError:
            return False


def is_mute_type(vcf_line, mute_list, pol):

    """
    takes pysam variant and determins if variant is of type listed in mutation list
    :param vcf_line: pysam variant
    :param mute_list: list
    :param pol: bool
    :return:
    """
    # strong = CG weak = AT
    base_types = {'A': 'W', 'T': 'W', 'C': 'S', 'G': 'S'}
    if mute_list is None:
        return True
    else:
        ref_base = base_types[vcf_line.ref]
        alt_base = base_types[vcf_line.alts[0]]
        if ref_base == alt_base == 'W':
            mutation_type = 'WW'
        elif ref_base == alt_base == 'S':
            mutation_type = 'SS'
        else:
            if pol is False:
                return False
            else:
                try:
                    anc_base = base_types[vcf_line.info['AA']]
                    if anc_base == ref_base:
                        mutation_type = anc_base + alt_base
                    else:
                        mutation_type = alt_base + ref_base
                except KeyError:
                    return False
        if mutation_type in mute_list:
            return True
        else:
            return False


def allele_num_ok(vcf_line, no_samples, multi):

    """
    checks to see if variant is biallelic
    :param vcf_line: pysam variant
    :param no_samples: int
    :param multi: bool
    :return: bool
    """

    if multi is False:
        pos_biallelic_freqs = [round(i/float(2*no_samples), 3) for i in range(1, 2*no_samples)]
        alt_allele_freq = round(vcf_line.info['AC'][0] / float(no_samples * 2), 3)
        if alt_allele_freq in pos_biallelic_freqs:
            return True
        else:
            return False
    else:
        return True


def is_auto(variant_line):

    """
    returns true if autosome
    :param variant_line: pysam variant
    :return: bool
    """
    sex_chromos = {'chrZ', 'Z', 'chrW', 'W', 'X', 'XHet', 'Y', 'YHet'}
    if variant_line.contig not in sex_chromos:
        return True
    else:
        return False


def is_heterozygous(variant_line):

    """
    returns False if all individuals at site are homozygous
    :param variant_line: pysam variant
    :return: bool
    """

    # makes a list of lengths of sets of each genotype and then makes that list of lengths a set and gets the length
    # ie if all individuals are homozygous it will return 1
    homozygous_value = len(set([len(set(x['GT'])) for x in variant_line.samples.values()]))

    if homozygous_value == 1:
        return False
    else:
        return True


def get_homozygosity(variant_line):

    """
    returns number of homozygous individuals
    :param variant_line: pysam variant
    :return: int
    """

    # makes a list of lengths of sets of each genotype ie: [1, 2, 1, 2] == [homo, hetero, homo, hetero]
    homozygous = [len(set(x['GT'])) for x in variant_line.samples.values()].count(1)

    return homozygous


def vcf2sfs(vcf_name, mode, chromo='ALL',
            start=None, stop=None, degen=None, mute_type=None, regions=None,
            fold=False, auto_only=False, multi_allelic=False, skip_hetero=False, bed=False, homozygosity=False,
            lengths=set([]), custom_info=None):

    """
    function that outputs site frequencies from vcf and is called in main()
    :param vcf_name: str
    :param mode: str
    :param chromo: str
    :param start: int
    :param stop: int
    :param degen: int
    :param mute_type: str
    :param regions: list
    :param fold: bool
    :param auto_only: bool
    :param multi_allelic: bool
    :param skip_hetero: bool
    :param bed: bool
    :param homozygosity: bool
    :param lengths: set
    :param custom_info: str
    :return: yields tuples or floats
    """

    # check commandline options
    if mode == 'indel' and fold is False:
        sys.exit('mode indel must be in conjunction with fold True')
    if mode == 'ins' and fold is True or mode == 'del' and fold is True:
        sys.exit('mode ins and mode del cannot be run in conjunction with fold True')
    if degen is not None and mode != 'snp':
        sys.exit('degen can only be specified in conjunction with mode snp')
    if mute_type is not None and mode != 'snp':
        sys.exit('mute_type can only be used with mode snp')
    if homozygosity and not bed:
        sys.exit('homozygosity can only be output in bed mode')
    if len(lengths) > 1 and mode == 'snp':
        sys.exit('len can only be specified with indel ins or del mode')

    # initiate pysam vcf
    if vcf_name != 'stdin':
        vcf_file = pysam.VariantFile(vcf_name)
    else:
        vcf_file = pysam.VariantFile('-')

    # loop through vcf
    if chromo == 'ALL' and vcf_name != 'stdin':
        vcf = vcf_file.fetch()
    elif vcf_name != 'stdin':
        if start is None and stop is None:
            vcf = vcf_file.fetch(chromo)
        else:
            vcf = vcf_file.fetch(chromo, start, stop)
    else:
        vcf = vcf_file

    number_samples = len(vcf_file.header.samples)

    for variant in vcf:
        # catch indels starting before first coord
        if start is not None and stop is not None:
            if variant.pos < start + 1:
                continue

        # gets relevant freq, minor or derived, see functions
        frequency = get_out_freq(variant, not fold, mode, number_samples)
        if frequency is None:  # skips when no freq returned, ie unpolarised or wrong var type
            continue

        # if custom info specified check if match
        if custom_info is not None:
            if not custom_info_expression(custom_info, variant):
                continue

        # gets variant region if regional sfs required
        falls_in_regions = in_regions(variant, regions)

        # gets degeneracy if required
        degen_ok = is_degen(variant, degen)

        # gets mutation type if required
        mutetype_ok = is_mute_type(variant, mute_type, not fold)

        # checks if is biallelic
        alleles_ok = allele_num_ok(variant, number_samples, multi_allelic)

        # checks if auto
        auto = is_auto(variant)
        if auto_only is True and auto is False:
            continue

        # checks if any individuals are heterozygous
        if skip_hetero and is_heterozygous(variant):
            continue

        # checks length specification is met
        if len(lengths) > 0:
            if indel_length(variant) not in lengths:
                continue

        # outputs if all criteria ok
        if falls_in_regions is True and degen_ok is True and mutetype_ok is True and alleles_ok is True:
            if bed:
                if homozygosity:
                    yield variant.contig, variant.start, variant.stop, frequency, get_homozygosity(variant)
                else:
                    yield variant.contig, variant.start, variant.stop, frequency
            else:
                yield frequency


# main call
def main():

    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-vcf', help='VCF file to extract sfs from, if not specified will read from standard in,'
                                     'but must contain the header', default='stdin')
    parser.add_argument('-chr', help='Chromosome to extract', default='ALL')
    parser.add_argument('-start', help='Start coord 0 based', type=int)
    parser.add_argument('-stop', help='End coord, 0 based', type=int)
    parser.add_argument('-region', help='Genomic regions to extract, default = ALL', action='append')
    parser.add_argument('-mode', help='Variant mode to run in', choices=['snp', 'ins', 'del', 'indel'], required=True)
    parser.add_argument('-degen', help='Degeneracy of coding SNPs to extract (must run with -mode snp',
                        choices=[0, 2, 3, 4], type=int)
    parser.add_argument('-mute_type', help='Mutation type, use only with mode -snp',
                        choices=['WW', 'SS', 'SW', 'WS'], action='append')
    parser.add_argument('-folded', help='If specified will output minor allele spectrum',
                        default=False, action='store_true')
    parser.add_argument('-auto_only', help='If specified will exclude sex chromosomes',
                        default=False, action='store_true')
    parser.add_argument('-multi_allelic', help='If specified will not restrict output to biallelic sites',
                        default=False, action='store_true')
    parser.add_argument('-skip_hetero', help='If specified will skip sites with heterozygous individuals',
                        default=False, action='store_true')
    parser.add_argument('-bed', help='If specified will output allele frequencies in bed format,'
                                     'each row specifying chromosome\tstart\tend\tallele_frequency',
                        default=False, action='store_true')
    parser.add_argument('-homozygosity', help='If specified will output homozygosity in bed format, each row specifying'
                                              ' chromosome\tstart\tend\tallele_frequency\thomozygosity',
                        default=False, action='store_true')
    parser.add_argument('-len', help='If specified then outputs INDELs of specified length', action='append',
                        default=[], type=int)
    args = parser.parse_args()

    # call to sfs function
    sfs = vcf2sfs(args.vcf, args.mode, chromo=args.chr, start=args.start, stop=args.stop, degen=args.degen,
                  mute_type=args.mute_type, regions=args.region, fold=args.folded, auto_only=args.auto_only,
                  multi_allelic=args.multi_allelic, skip_hetero=args.skip_hetero, bed=args.bed,
                  homozygosity=args.homozygosity, lengths=set(args.len))

    for site in sfs:
        if args.bed:
            print(*list(site), sep='\t')
        else:
            print(site)

if __name__ == '__main__':
    main()
