#!/usr/bin/env python

# ENCODE DCC bowtie wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import re
import argparse
import multiprocessing
from encode_common_genomic import *

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

# WDL glob() globs in an alphabetical order
# so R1 and R2 can be switched, which results in an
# unexpected behavior of a workflow.
# so we already prepended merge_fastqs_'end'_ (R1 or R2) 
# to the basename of original filename in 'trim_adapter' task.
# now it's time to strip it.
def strip_merge_fastqs_prefix(fastq):
    return re.sub(r'^merge\_fastqs\_R\d\_','',str(fastq))

def make_read_length_file(fastq, out_dir):
    basename = os.path.basename(strip_ext_fastq(fastq))
    prefix = os.path.join(out_dir,
        strip_merge_fastqs_prefix(basename))
    txt = '{}.read_length.txt'.format(prefix)
    read_length = get_read_length(fastq)
    with open(txt,'w') as fp:
        fp.write(str(read_length))
    return txt

def bowtie_se(fastq, param, ref_index_prefix, nth, out_dir):
    basename = os.path.basename(strip_ext_fastq(fastq))
    prefix = os.path.join(out_dir,
        strip_merge_fastqs_prefix(basename))        
    unzipped_fq = '{}.tmp'.format(prefix)
    bam = '{}.bam'.format(prefix)
    cmd = 'zcat -f {} > {}'.format(fastq, unzipped_fq)
    run_shell_cmd(cmd)

    cmd2 = 'bowtie -S -p {} {} {} {} | '
    cmd2 += 'samtools view -F 4 -Su - | samtools sort - {}'
    cmd2 = cmd2.format(
        nth,
        param,
        ref_index_prefix,
        unzipped_fq,
        prefix)
    run_shell_cmd(cmd2)

    rm_f(unzipped_fq)

    return bam

def bowtie_pe(fastq1, fastq2, param, ref_index_prefix, nth, out_dir):
    raise NotImplementedError

def chk_bowtie_index(prefix):    
    index_sa = '{}.1.ebwt'.format(prefix)
    if not os.path.exists(index_sa):
        raise Exception("bowtie index does not exists. "+
            "Prefix = {}".format(prefix))

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
    if args.bowtie_index_prefix_or_tar.endswith('.tar'):
        log.info('Unpacking bowtie index tar...')
        tar = args.bowtie_index_prefix_or_tar
        # untar
        untar(tar, args.out_dir)
        bowtie_index_prefix = os.path.join(
            args.out_dir, os.path.basename(strip_ext_tar(tar)))
        temp_files.append('{}.*'.format(
            bowtie_index_prefix))
    else:
        bowtie_index_prefix = args.bowtie_index_prefix_or_tar

    # check if bowties indices are unpacked on out_dir
    chk_bowtie_index(bowtie_index_prefix)

    # bowtie
    log.info('Running bowtie...')
    if args.paired_end:
        bam = bowtie_pe(args.fastqs[0], args.fastqs[1], args.param,
            bowtie_index_prefix, args.nth, args.out_dir)
    else:
        bam = bowtie_se(args.fastqs[0], args.param,
            bowtie_index_prefix, args.nth, args.out_dir)

    # initialize multithreading
    log.info('Initializing multi-threading...')
    num_process = min(2,args.nth)
    log.info('Number of threads={}.'.format(num_process))
    pool = multiprocessing.Pool(num_process)
    
    log.info('Running samtools index...')
    ret_val1 = pool.apply_async(
        samtools_index, (bam, args.out_dir))

    log.info('Running samtools flagstat...')
    ret_val2 = pool.apply_async(
        samtools_flagstat, (bam, args.out_dir))

    bai = ret_val1.get(BIG_INT)
    flagstat_qc = ret_val2.get(BIG_INT)

    log.info('Closing multi-threading...')
    pool.close()
    pool.join()

    log.info('Removing temporary files...')
    rm_f(temp_files)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')

if __name__=='__main__':
    main()
