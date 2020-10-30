#!/usr/bin/env python

# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
from encode_lib_common import (
    mkdir_p, log, ls_l, rm_f, strip_ext_fastq)
from encode_lib_genomic import (
    get_read_length, remove_chrs_from_bam, samstat, samtools_index)


def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE post align',
                                     description='')
    parser.add_argument('fastq', type=str,
                        help='Path for FASTQ R1')
    parser.add_argument('bam', type=str,
                        help='Path for BAM')
    parser.add_argument(
        '--chrsz', type=str,
        help='2-col chromosome sizes file. If not given then '
             'SAMstats on mito-free BAM will not be calcaulted.')
    parser.add_argument('--mito-chr-name', default='chrM',
                        help='Mito chromosome name.')
    parser.add_argument('--nth', type=int, default=1,
                        help='Number of threads to parallelize.')
    parser.add_argument('--mem-gb', type=float,
                        help='Max. memory for samtools sort in GB. '
                        'It should be total memory for this task (not memory per thread).')
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO',
                                 'WARNING', 'CRITICAL', 'ERROR',
                                 'CRITICAL'],
                        help='Log level')
    args = parser.parse_args()

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args


def make_read_length_file(fastq, out_dir):
    basename = os.path.basename(strip_ext_fastq(fastq))
    prefix = os.path.join(out_dir, basename)
    txt = '{}.read_length.txt'.format(prefix)
    read_length = get_read_length(fastq)
    with open(txt, 'w') as fp:
        fp.write(str(read_length))
    return txt


def main():
    # read params
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # generate read length file
    log.info('Generating read length file...')
    make_read_length_file(
        args.fastq, args.out_dir)

    log.info('Running samtools index...')
    samtools_index(args.bam, args.nth, args.out_dir)

    log.info('SAMstats on raw BAM...')
    samstat(args.bam, args.nth, args.mem_gb, args.out_dir)

    if args.chrsz:
        log.info('SAMstats on non-mito BAM...')
        non_mito_out_dir = os.path.join(args.out_dir, 'non_mito')
        mkdir_p(non_mito_out_dir)
        non_mito_bam = remove_chrs_from_bam(args.bam, [args.mito_chr_name],
                                            args.chrsz,
                                            args.nth,
                                            non_mito_out_dir)
        samstat(non_mito_bam, args.nth, args.mem_gb, non_mito_out_dir)
        rm_f(non_mito_bam)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')


if __name__ == '__main__':
    main()
