#!/usr/bin/env python
# ENCODE DCC adapter trimmer wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
import multiprocessing
import copy
from encode_lib_common import (
    copy_f_to_dir, log, ls_l, mkdir_p, read_tsv, rm_f,
    run_shell_cmd, strip_ext_fastq)
from encode_task_merge_fastq import merge_fastqs


def parse_arguments(debug=False):
    parser = argparse.ArgumentParser(prog='ENCODE DCC adapter trimmer.',
                                        description='')
    parser.add_argument('fastqs', nargs='+', type=str,
                        help='TSV file path or list of FASTQs. \
                            FASTQs must be compressed with gzip (with .gz). \
                            Use TSV for multiple fastqs to be merged later. \
                            row=merge_id, col=end_id).')
    parser.add_argument('--adapter', type=str, default='AGATCGGAAGAGCACACGTCTGGATCCACGACGCTCTTCC',
                        help='Adapter sequence.')
    parser.add_argument('--barcodes', type=str, default='',
                        help='Barcodes (comma-delimited if multiple).')
    parser.add_argument(
        '--cutadapt-param', type=str, default='-e 0.2 -m 22 -O 4',
        help='Parameters for cutadapt. CLI: double-quote it with a leading whitespace. '
             'default: -e 0.2 -m 22 -O 4. '
             'Important cutadapt paramters: '
             '-e: Maximum allowed adapter error rate (no. errors divided by the length of the matching adapter region), '
             '-m: Minimum trim length (throwing away processed reads shorter than this), '
             '-O: Minimum overlap.'
    )
    parser.add_argument(
        '--nimnexus-param', type=str, default='-t 0 -k 18 -r 5',
        help='Parameters for nimnexus. CLI: double-quote it with a leading whitespace. '
             'default: -t 0 -k 18 -r 5. '
             'Important nimnexus paramters: '
             '-t: Pre-trim all reads by this length before processing, '
             '-k: Minimum number of bases required after barcode to keep read. '
             '-r: Number of bases at the start of each read used for random barcode.'
    )

    parser.add_argument('--paired-end', action="store_true",
                        help='Paired-end FASTQs.')
    parser.add_argument('--nth', type=int, default=1,
                        help='Number of threads to parallelize.')
    parser.add_argument('--out-dir', default='', type=str,
                            help='Output directory.')
    parser.add_argument('--log-level', default='INFO', 
                        choices=['NOTSET','DEBUG','INFO',
                            'WARNING','CRITICAL','ERROR','CRITICAL'],
                        help='Log level')
    args = parser.parse_args()

    # parse fastqs command line
    if args.fastqs[0].endswith('.gz') or args.fastqs[0].endswith('.fastq') or \
            args.fastqs[0].endswith('.fq'):  # it's fastq
        args.fastqs = [[f] for f in args.fastqs]  # make it a matrix
    else:  # it's TSV
        args.fastqs = read_tsv(args.fastqs[0])

    adapters = copy.deepcopy(args.fastqs)
    for i, adapters_ in enumerate(adapters):
        for j, adapter in enumerate(adapters_):
            adapters[i][j] = args.adapter

    # check if fastqs, adapers have same/correct dimension
    if len(adapters) != len(args.fastqs):
        raise argparse.ArgumentTypeError(
            'fastqs and adapters dimension mismatch.')
    for i, fastqs in enumerate(args.fastqs):
        if args.paired_end and len(fastqs) != 2:
            raise argparse.ArgumentTypeError(
                'Need 2 fastqs per replicate for paired end.')
        if not args.paired_end and len(fastqs) != 1:
            raise argparse.ArgumentTypeError(
                'Need 1 fastq per replicate for single end.')
        if len(fastqs) != len(adapters[i]):
            raise argparse.ArgumentTypeError(
                'fastqs and adapters dimension mismatch.')
            
    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args, adapters


def trim_adapter_se(fastq, adapter, barcodes, cutadapt_param, nimnexus_param, nth, out_dir):
    if adapter:
        prefix = os.path.join(out_dir,
            os.path.basename(strip_ext_fastq(fastq)))
        trimmed = '{}.trim.fastq.gz'.format(prefix)

        adapter1 = adapter[1:]
        adapter2 = adapter[2:]
        adapter3 = adapter[3:]
        adapter4 = adapter[4:]

        cmd = 'pigz -cd {fastq} | '.format(fastq=fastq)
        if barcodes:
            cmd += 'nimnexus trim {barcodes} {nimnexus_param} | '.format(
                barcodes=barcodes,
                nimnexus_param=nimnexus_param
            )
        cmd += (
            'parallel -j {nth} --pipe -L 4 --block 10M '
            'cutadapt {cutadapt_param} -a {adapter} -a {adapter1} '
            '-a {adapter2} -a {adapter3} -a {adapter4} - | pigz -nc > {trimmed}'
        ).format(
            nth=nth,
            cutadapt_param=cutadapt_param,
            adapter=adapter,
            adapter1=adapter1,
            adapter2=adapter2,
            adapter3=adapter3,
            adapter4=adapter4,
            trimmed=trimmed,
        )

        run_shell_cmd(cmd)
        return trimmed
    else:
        return copy_f_to_dir(fastq, out_dir)


def trim_adapter_pe(fastq1, fastq2, adapter1, adapter2, barcodes,
        cutadapt_param, nimnexus_param, nth, out_dir):
    raise NotImplementedError


def main():
    # read params
    args, adapters = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # declare temp arrays
    temp_files = [] # files to deleted later at the end

    trimmed_fastqs_R1 = []
    trimmed_fastqs_R2 = []

    log.info('Trimming adapters and merging fastqs...')
    for i in range(len(args.fastqs)):
        # for each fastq to be merged later
        fastqs = args.fastqs[i] # R1 and R2
        adapters_ = adapters[i]
        if args.paired_end:
            fastqs = trim_adapter_pe(
                fastqs[0], fastqs[1], 
                adapters_[0], adapters_[1],
                args.barcodes,
                args.cutadapt_param,
                args.nimnexus_param,
                args.nth,
                args.out_dir)
            trimmed_fastqs_R1.append(fastqs[0])
            trimmed_fastqs_R2.append(fastqs[1])
        else:
            fastq = trim_adapter_se(
                fastqs[0],
                adapters_[0],
                args.barcodes,
                args.cutadapt_param,
                args.nimnexus_param,
                args.nth,
                args.out_dir)
            trimmed_fastqs_R1.append(fastq)

    log.info('Merging fastqs...')
    log.info('R1 to be merged: {}'.format(trimmed_fastqs_R1))
    merge_fastqs(trimmed_fastqs_R1, 'R1', args.out_dir)
    if args.paired_end:
        log.info('R2 to be merged: {}'.format(trimmed_fastqs_R2))
        merge_fastqs(trimmed_fastqs_R2, 'R2', args.out_dir)

    temp_files.extend(trimmed_fastqs_R1)
    temp_files.extend(trimmed_fastqs_R2)

    log.info('Removing temporary files...')
    rm_f(temp_files)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)
    ls_l(os.path.join(args.out_dir, 'R1'))
    if args.paired_end:
        ls_l(os.path.join(args.out_dir, 'R2'))

    log.info('All done.')


if __name__=='__main__':
    main()