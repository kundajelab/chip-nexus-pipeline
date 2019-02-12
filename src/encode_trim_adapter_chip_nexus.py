#!/usr/bin/env python

# ENCODE DCC adapter trimmer wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
import multiprocessing
import copy
from encode_common import *

def parse_arguments(debug=False):
    parser = argparse.ArgumentParser(prog='ENCODE DCC adapter trimmer.',
                                        description='')
    parser.add_argument('fastqs', nargs='+', type=str,
                        help='TSV file path or list of FASTQs. \
                            FASTQs must be compressed with gzip (with .gz). \
                            Use TSV for multiple fastqs to be merged later. \
                            row=merge_id, col=end_id).')
    parser.add_argument('--min-overlap', type=int, default=4,
                        help='Minimum overlap for cutadapt -O.')
    parser.add_argument('--min-trim-len', type=int, default=22,
                        help='Minimum trim length for cutadapt -m \
                            (throwing away processed reads shorter than this).')
    parser.add_argument('--err-rate', type=float, default=0.2,
                        help='Maximum allowed adapter error rate for cutadapt -e \
                            (no. errors divided by the length \
                            of the matching adapter region).')
    parser.add_argument('--trim', type=int, default=0,
                        help='Pre-trim all reads by this length before processing '
                        '(nimnexus -t).')
    parser.add_argument('--keep', type=int, default=18,
                        help='Minimum number of bases required after barcode to keep read '
                        '(nimnexus -k).')
    parser.add_argument('--randombarcode', type=int, default=5,
                        help='Number of bases at the start of each read used for random barcode '
                        '(nimnexus -r).')
    parser.add_argument('--adapter', type=str, default='AGATCGGAAGAGCACACGTCTGGATCCACGACGCTCTTCC',
                        help='Adapter sequence.')
    parser.add_argument('--barcodes', type=str, default='',
                        help='Barcodes (comma-delimited if multiple).')
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
    if args.fastqs[0].endswith('.gz'): # it's fastq
         args.fastqs = [[f] for f in args.fastqs] # make it a matrix
    else: # it's TSV
        args.fastqs = read_tsv(args.fastqs[0])
    
    adapters = copy.deepcopy(args.fastqs)
    for i, adapters_ in enumerate(adapters):
        for j, adapter in enumerate(adapters_):
            adapters[i][j] = args.adapter

    # check if fastqs, adapers have same/correct dimension
    if len(adapters)!=len(args.fastqs):
        raise argparse.ArgumentTypeError(
            'fastqs and adapters dimension mismatch.')
    for i, fastqs in enumerate(args.fastqs):
        if args.paired_end and len(fastqs)!=2:
            raise argparse.ArgumentTypeError(
                'Need 2 fastqs per replicate for paired end.')
        if not args.paired_end and len(fastqs)!=1:
            raise argparse.ArgumentTypeError(
                'Need 1 fastq per replicate for single end.')
        if len(fastqs)!=len(adapters[i]):
            raise argparse.ArgumentTypeError(
                'fastqs and adapters dimension mismatch.')
            
    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args, adapters

def trim_adapter_se(fastq, adapter, barcodes, min_overlap, min_trim_len, err_rate,
    trim, keep, randombarcode, nth, out_dir):
    if adapter:
        prefix = os.path.join(out_dir,
            os.path.basename(strip_ext_fastq(fastq)))
        trimmed = '{}.trim.fastq.gz'.format(prefix)

        adapter1 = adapter[1:]
        adapter2 = adapter[2:]
        adapter3 = adapter[3:]
        adapter4 = adapter[4:]

        # cmd = 'zcat -f {} | '
        # if barcodes:
        #     cmd += 'nimnexus trim {} -t {} -k {} -r {} | '.format(barcodes,
        #         trim,
        #         keep,
        #         randombarcode)
        # cmd += 'cutadapt {} {} -e {} -a {} -a {} -a {} -a {} -a {} - | gzip -nc > {}'
        # cmd = cmd.format(
        #     fastq,            
        #     '-O {}'.format(min_overlap) if min_overlap > 0 else '',
        #     '-m {}'.format(min_trim_len) if min_trim_len > 0 else '',
        #     err_rate,
        #     adapter, adapter1, adapter2, adapter3, adapter4,
        #     trimmed)     

        cmd = 'pigz -cd {} | '
        if barcodes:
            cmd += 'nimnexus trim {} -t {} -k {} -r {} | '.format(barcodes,
                trim,
                keep,
                randombarcode)
        cmd += 'parallel -j {} --pipe -L 4 --block 10M '
        cmd += 'cutadapt {} {} -e {} -a {} -a {} -a {} -a {} -a {} - | pigz -nc > {}'
        cmd = cmd.format(
            fastq,
            nth,
            '-O {}'.format(min_overlap) if min_overlap > 0 else '',
            '-m {}'.format(min_trim_len) if min_trim_len > 0 else '',
            err_rate,
            adapter, adapter1, adapter2, adapter3, adapter4,
            trimmed)     

        run_shell_cmd(cmd)
        return trimmed
    else:
        # make hard link
        linked = os.path.join(out_dir,
            os.path.basename(fastq))
        os.link(fastq, linked)
        return linked        

def trim_adapter_pe(fastq1, fastq2, adapter1, adapter2, barcodes,
        min_overlap, min_trim_len, err_rate,
        trim, keep, randombarcode, nth, out_dir):
    raise NotImplementedError

# WDL glob() globs in an alphabetical order
# so R1 and R2 can be switched, which results in an
# unexpected behavior of a workflow
# so we prepend merge_fastqs_'end'_ (R1 or R2)
# to the basename of original filename
def merge_fastqs(fastqs, end, out_dir):
    basename = os.path.basename(strip_ext_fastq(fastqs[0]))
    prefix = os.path.join(out_dir,
        'merge_fastqs_{}_{}'.format(end, basename))
    merged = '{}.merged.fastq.gz'.format(prefix)

    if len(fastqs)>1:
        cmd = 'zcat -f {} | gzip -nc > {}'.format(
            ' '.join(fastqs),
            merged)
        run_shell_cmd(cmd)
        return merged
    else:
        return hard_link(fastqs[0], merged)

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
                args.min_overlap,
                args.min_trim_len,
                args.err_rate,
                args.trim, args.keep, args.randombarcode,
                args.nth,
                args.out_dir)
            trimmed_fastqs_R1.append(fastqs[0])
            trimmed_fastqs_R2.append(fastqs[1])
        else:
            fastq = trim_adapter_se(
                fastqs[0],
                adapters_[0],
                args.barcodes,
                args.min_overlap,
                args.min_trim_len,
                args.err_rate,
                args.trim, args.keep, args.randombarcode,
                args.nth,
                args.out_dir)
            trimmed_fastqs_R1.append(fastq)

    log.info('R1 to be merged: {}'.format(trimmed_fastqs_R1))
    R1_merged = merge_fastqs(trimmed_fastqs_R1, 'R1', args.out_dir)
    temp_files.extend(trimmed_fastqs_R1)

    if args.paired_end:
        log.info('R2 to be merged: {}'.format(trimmed_fastqs_R2))
        R2_merged = merge_fastqs(trimmed_fastqs_R2, 'R2', args.out_dir)
        temp_files.extend(trimmed_fastqs_R2)

    log.info('Removing temporary files...')
    rm_f(temp_files)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('Checking if output is empty...')
    assert_file_not_empty(R1_merged)

    log.info('All done.')

if __name__=='__main__':
    main()