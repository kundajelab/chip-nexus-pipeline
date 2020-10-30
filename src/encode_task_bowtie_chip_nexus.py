#!/usr/bin/env python

# ENCODE DCC bowtie wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import re
import argparse
from encode_lib_common import (
    log, ls_l, mkdir_p, rm_f, run_shell_cmd, strip_ext_fastq, strip_ext_tar,
    untar)
from encode_lib_genomic import samtools_sort, bam_is_empty

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE DCC bowtie aligner for ChIP-nexus.',
                                        description='')
    parser.add_argument('bowtie_index_prefix_or_tar', type=str,
                        help='Path for prefix (or a tarball .tar) \
                            for reference bowtie index. \
                            Prefix must be like [PREFIX].sa. \
                            Tar ball must be packed without compression \
                            and directory by using command line \
                            "tar cvf [TAR] [TAR_PREFIX].*".')
    parser.add_argument('fastqs', nargs='+', type=str,
                        help='List of FASTQs (R1 and R2). \
                            FASTQs must be compressed with gzip (with .gz).')
    parser.add_argument('--paired-end', action="store_true",
                        help='Paired-end FASTQs.')
    parser.add_argument('--param', type=str, default='--chunkmbs 512 -k 1 -m 1 -v 2 --best --strata',
                        help='Bowtie command line parameter.')
    parser.add_argument('--nth', type=int, default=1,
                        help='Number of threads to parallelize.')
    parser.add_argument('--mem-gb', type=float,
                        help='Max. memory for samtools sort in GB. '
                        'It should be total memory for this task (not memory per thread).')
    parser.add_argument('--out-dir', default='', type=str,
                            help='Output directory.')
    parser.add_argument('--log-level', default='INFO', 
                        choices=['NOTSET','DEBUG','INFO',
                            'WARNING','CRITICAL','ERROR','CRITICAL'],
                        help='Log level')
    args = parser.parse_args()

    # check if fastqs have correct dimension
    if args.paired_end and len(args.fastqs)!=2:
        raise argparse.ArgumentTypeError('Need 2 fastqs for paired end.')
    if not args.paired_end and len(args.fastqs)!=1:
        raise argparse.ArgumentTypeError('Need 1 fastq for single end.')

    log.setLevel(args.log_level)
    log.info(sys.argv)    
    return args


def bowtie_se(fastq, ref_index_prefix, param, nth, mem_gb, out_dir):
    basename = os.path.basename(strip_ext_fastq(fastq))
    prefix = os.path.join(out_dir,
        strip_merge_fastqs_prefix(basename))        
    unzipped_fq = '{}.tmp'.format(prefix)
    tmp_bam = '{}.bam'.format(prefix)

    run_shell_cmd(
        'bowtie -S -p {nth} {param} {index} {fastq} | '
        'samtools view -1 -F 4 -S /dev/stdin > {tmp_bam}'.format(
            nth=nth,
            param=param,
            index=ref_index_prefix,
            fastq=unzipped_fq,
            tmp_bam=tmp_bam,                        
        )
    )

    bam = samtools_sort(tmp_bam, nth, mem_gb, out_dir)
    rm_f([tmp_bam, unzipped_fq])

    return bam


def bowtie_pe(fastq1, fastq2, ref_index_prefix, param, nth, mem_gb, out_dir):
    raise NotImplementedError


def chk_bowtie_index(prefix):    
    index_sa = '{}.1.ebwt'.format(prefix)
    if not os.path.exists(index_sa):
        raise Exception("bowtie index does not exists. "+
            "Prefix = {}".format(prefix))


def find_bowtie_index_prefix(d):
    """
    Returns:
        prefix of BWA index. e.g. returns PREFIX if PREFIX.sa exists
    Args:
        d: directory to search for .1.bt2 or .1.bt2l file
    """
    if d == '':
        d = '.'
    for f in os.listdir(d):
        if f.endswith('.1.ebwt'):
            return re.sub('\.1\.ebwt$', '', f)


def main():
    # read params
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # declare temp arrays
    temp_files = [] # files to deleted later at the end

    # generate read length file
    log.info('Generating read length file...')
    R1_read_length_file = make_read_length_file(
                            args.fastqs[0], args.out_dir)
    if args.paired_end:
        R2_read_length_file = make_read_length_file(
                            args.fastqs[1], args.out_dir)
    
    # if bowtie index is tarball then unpack it
    if args.bowtie_index_prefix_or_tar.endswith('.tar') or \
            args.bowtie_index_prefix_or_tar.endswith('.tar.gz'):    
        log.info('Unpacking bowtie index tar...')
        tar = args.bowtie_index_prefix_or_tar
        # untar
        untar(tar, args.out_dir)
        bowtie_index_prefix = find_bowtie_index_prefix(args.out_dir)
        temp_files.append('{}.*'.format(
            bowtie_index_prefix))
    else:
        bowtie_index_prefix = args.bowtie_index_prefix_or_tar

    # check if bowties indices are unpacked on out_dir
    chk_bowtie_index(bowtie_index_prefix)

    # bowtie
    log.info('Running bowtie...')
    if args.paired_end:
        bam = bowtie_pe(
            args.fastqs[0], args.fastqs[1],
            bowtie_index_prefix,
            args.param,
            args.nth, args.mem_gb,
            args.out_dir)
    else:
        bam = bowtie_se(
            args.fastqs[0],
            bowtie_index_prefix,
            args.param,
            args.nth, args.mem_gb,
            args.out_dir)

    log.info('Removing temporary files...')
    print(temp_files)
    rm_f(temp_files)

    log.info('Checking if BAM file is empty...')
    if bam_is_empty(bam, args.nth):
        raise ValueError('BAM file is empty, no reads found.')

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')

if __name__=='__main__':
    main()
