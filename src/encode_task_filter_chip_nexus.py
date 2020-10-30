#!/usr/bin/env python

# ENCODE DCC filter wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
import multiprocessing
from encode_lib_common import (
    copy_f_to_dir, log, ls_l, mkdir_p, rm_f, run_shell_cmd, strip_ext,
    strip_ext_bam)
from encode_lib_genomic import (
    locate_picard, remove_chrs_from_bam, samstat, samtools_index,
    samtools_name_sort, bam_is_empty,
    get_samtools_res_param)
from encode_task_filter import (
    rm_unmapped_lowq_reads_se,
)

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE DCC filter.',
                                        description='')
    parser.add_argument('bam', type=str,
                        help='Path for raw BAM file.')
    parser.add_argument('--mapq-thresh', default=30, type=int,
                        help='Threshold for low MAPQ reads removal.')
    parser.add_argument('--no-dup-removal', action="store_true",
                        help='No dupe reads removal when filtering BAM.')
    parser.add_argument('--paired-end', action="store_true",
                        help='Paired-end BAM.')
    parser.add_argument('--multimapping', default=0, type=int,
                        help='Multimapping reads.')
    parser.add_argument(
        '--filter-chrs', nargs='*',
        help='Chromosomes to be filtered for final (nodup/filt) BAM.')
    parser.add_argument('--chrsz', type=str,
                        help='2-col chromosome sizes file.')
    parser.add_argument('--mito-chr-name', default='chrM',
                        help='Mito chromosome name.')
    parser.add_argument('--nth', type=int, default=1,
                        help='Number of threads to parallelize.')
    parser.add_argument('--mem-gb', type=float,
                        help='Max. memory for samtools sort in GB. '
                        'It should be total memory for this task (not memory per thread).')
    parser.add_argument('--picard-java-heap',
                        help='Picard\'s Java max. heap: java -jar picard.jar '
                             '-Xmx[MAX_HEAP]')
    parser.add_argument('--out-dir', default='', type=str,
                            help='Output directory.')
    parser.add_argument('--log-level', default='INFO', 
                        choices=['NOTSET','DEBUG','INFO',
                            'WARNING','CRITICAL','ERROR','CRITICAL'],
                        help='Log level')
    args = parser.parse_args()

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args


def rm_dup_se(filt_bam, nth, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_bam(filt_bam)))
    # strip extension appended in the previous step
    prefix = strip_ext(prefix,'filt') 
    nodup_bam = '{}.nodup.bam'.format(prefix)

    cmd1 = 'nimnexus dedup -t {} {} | samtools view - -b > {}'
    cmd1 = cmd1.format(
        nth,
        filt_bam,
        nodup_bam)
    run_shell_cmd(cmd1)
    return nodup_bam


def pbc_qc_se(filt_bam, nodup_bam, mito_chr_name, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_bam(filt_bam)))
    pbc_qc = '{}.pbc.qc'.format(prefix)

    mt = int(run_shell_cmd('samtools view {} | grep -v "\\b{}\\b" | wc -l'.format(filt_bam, mito_chr_name)))
    m0 = int(run_shell_cmd('samtools view {} | grep -v "\\b{}\\b" | wc -l'.format(nodup_bam, mito_chr_name)))
    m0_mt = 0 if mt==0.0 else m0/float(mt)
    with open(pbc_qc, 'w') as fp:
        fp.write("{}\t{}\tN/A\tN/A\t{}\tN/A\tN/A\n".format(mt, m0, m0_mt))

    return pbc_qc


def pbc_qc_pe(filt_bam, nodup_bam, mito_chr_name, nth, out_dir):
    raise NotImplementedError


def main():
    # filt_bam - dupmark_bam - nodup_bam
    #          \ dup_qc      \ pbc_qc

    # read params
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # declare temp arrays
    temp_files = [] # files to deleted later at the end

    log.info('Removing unmapped/low-quality reads...')
    if args.paired_end:
        raise NotImplementedError('PE is not supported.')
    else:
        filt_bam = rm_unmapped_lowq_reads_se(
                args.bam, args.multimapping, args.mapq_thresh, 
                args.nth, args.mem_gb, args.out_dir)

    log.info('Checking if filtered BAM file is empty...')

    if bam_is_empty(filt_bam, args.nth):
        help_msg = (
            'No reads found in filtered BAM. '
            'Low quality sample? '
            'Or no reads passing criteria "samtools view -F 1804"? '
            'Check samtools flags at '
            'https://broadinstitute.github.io/picard/explain-flags.html. '
        )
        if args.paired_end:
            help_msg += (
                'Or is this truely PE BAM? '
                'All unpaired SE reads could be removed by "samtools view -f 2". '
            )
        raise ValueError(help_msg)

    if args.no_dup_removal:
        nodup_bam = filt_bam        
    else:
        temp_files.append(filt_bam)        
        log.info('Removing dupes...')
        if args.paired_end:
            raise NotImplementedError('PE is not supported.')
        else:
            nodup_bam = rm_dup_se(
                        filt_bam, args.nth, args.out_dir)
        temp_files.append(filt_bam)

    if len(args.filter_chrs) > 0:
        final_bam = remove_chrs_from_bam(nodup_bam, args.filter_chrs,
                                         args.chrsz, args.nth,
                                         args.out_dir)
        temp_files.append(nodup_bam)
    else:
        final_bam = nodup_bam

    log.info('Checking if final BAM file is empty...')
    if bam_is_empty(final_bam, args.nth):
        raise ValueError(
            'No reads found in final (filtered/deduped) BAM. '
            'Low quality sample? '
            'Or BAM with duplicates only? '
        )

    log.info('samtools index (final_bam)...')
    samtools_index(final_bam, args.nth, args.out_dir)

    log.info('samstat...')
    samstat(final_bam, args.nth, args.mem_gb, args.out_dir)

    log.info('Generating PBC QC log...')
    if args.paired_end:
        raise NotImplementedError('PE is not supported.')
    else:
        pbc_qc_se(filt_bam, nodup_bam, args.mito_chr_name, args.out_dir)

    log.info('samtools index (raw bam)...')
    bam = copy_f_to_dir(args.bam, args.out_dir)
    bai = samtools_index(bam, args.nth, args.out_dir)
    temp_files.extend([bam, bai])

    log.info('Removing temporary files...')
    rm_f(temp_files)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')


if __name__=='__main__':
    main()
