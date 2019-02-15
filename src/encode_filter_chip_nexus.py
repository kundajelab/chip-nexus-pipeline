#!/usr/bin/env python

# ENCODE DCC filter wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
import multiprocessing
from encode_common_genomic import *

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
    parser.add_argument('--mito-chr-name', default='chrM',
                        help='Mito chromosome name.')
    parser.add_argument('--nth', type=int, default=1,
                        help='Number of threads to parallelize.')
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

def rm_unmapped_lowq_reads_se(bam, multimapping, mapq_thresh, nth, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_bam(bam)))
    filt_bam = '{}.filt.bam'.format(prefix)

    if multimapping:        
        # qname_sort_bam = samtools_name_sort(bam, nth, out_dir)
        qname_sort_bam = sambamba_name_sort(bam, nth, out_dir)

        cmd2 = 'samtools view -h {} | '
        cmd2 += '$(which assign_multimappers.py) -k {} | '
        cmd2 += 'samtools view -F 1804 -Su /dev/stdin | '

        # cmd2 += 'samtools sort /dev/stdin -o {} -T {} -@ {}'
        # cmd2 = cmd2.format(
        #     qname_sort_bam,
        #     multimapping,
        #     filt_bam,
        #     prefix,
        #     nth)
        cmd2 += 'sambamba sort /dev/stdin -o {} -t {}'
        cmd2 = cmd2.format(
            qname_sort_bam,
            multimapping,
            filt_bam,
            nth)
        run_shell_cmd(cmd2)
        rm_f(qname_sort_bam) # remove temporary files
    else:
        cmd = 'samtools view -F 1804 -q {} -u {} | '
        cmd += 'samtools sort /dev/stdin -o {} -T {} -@ {}'
        cmd = cmd.format(
            mapq_thresh,
            bam,
            filt_bam,
            prefix,
            nth)
        run_shell_cmd(cmd)

    return filt_bam

def rm_unmapped_lowq_reads_pe(bam, multimapping, mapq_thresh, nth, out_dir):
    raise NotImplementedError

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

def rm_dup_pe(filt_bam, nth, out_dir):
    raise NotImplementedError

def pbc_qc_se(filt_bam, nodup_bam, mito_chr_name, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_bam(filt_bam)))
    pbc_qc = '{}.pbc.qc'.format(prefix)

    mt = int(run_shell_cmd('samtools view {} | grep -v "\\b{}\\b" | wc -l'.format(filt_bam, mito_chr_name)))
    m0 = int(run_shell_cmd('samtools view {} | grep -v "\\b{}\\b" | wc -l'.format(nodup_bam, mito_chr_name)))
    m0_mt = 0 if mt==0.0 else m0/float(mt)
    with open(pbc_qc, 'w') as fp:
        fp.write("{}\t{}\tN/A\tN/A\t{}\tN/A\tN/A\n".format(mt, m0, m0_mt))

    # cmd2 = 'bedtools bamtobed -i {} | '
    # cmd2 += 'awk \'BEGIN{{OFS="\\t"}}{{print $1,$2,$3,$6}}\' | '
    # cmd2 += 'grep -v "^{}\\b" | sort | uniq -c | '
    # cmd2 += 'awk \'BEGIN{{mt=0;m0=0;m1=0;m2=0}} ($1==1){{m1=m1+1}} '
    # cmd2 += '($1==2){{m2=m2+1}} {{m0=m0+1}} {{mt=mt+$1}} END{{m1_m2=-1.0; '
    # cmd2 += 'if(m2>0) m1_m2=m1/m2; m0_mt=0; if (mt>0) m0_mt=m0/mt; m1_m0=0; if (m0>0) m1_m0=m1/m0; '
    # cmd2 += 'printf "%d\\t%d\\t%s\\t%s\\t%f\\t%s\\t%s\\n",'
    # cmd2 += 'mt,m0,"N/A","N/A",m0_mt,"N/A","N/A"}}\' > {}'
    # cmd2 = cmd2.format(
    #     bam,
    #     mito_chr_name,
    #     pbc_qc)
    # run_shell_cmd(cmd2)
    return pbc_qc

def pbc_qc_pe(filt_bam, nodup_bam, mito_chr_name, nth, out_dir):
    raise NotImplementedError

# if --no-dup-removal is on, 
# Cromwell/WDL wants to have a empty file 
# for output { File pbc_qc, File dup_qc }

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
        filt_bam = rm_unmapped_lowq_reads_pe(
                args.bam, args.multimapping, args.mapq_thresh, 
                args.nth, args.out_dir)
    else:
        filt_bam = rm_unmapped_lowq_reads_se(
                args.bam, args.multimapping, args.mapq_thresh, 
                args.nth, args.out_dir)

    if args.no_dup_removal:
        nodup_bam = filt_bam        
    else:
        log.info('Removing dupes...')
        if args.paired_end:
            nodup_bam = rm_dup_pe(
                        filt_bam, args.nth, args.out_dir)
        else:
            nodup_bam = rm_dup_se(
                        filt_bam, args.nth, args.out_dir)
        temp_files.append(filt_bam)

    # initialize multithreading
    log.info('Initializing multi-threading...')
    num_process = min(3,args.nth)
    log.info('Number of threads={}.'.format(num_process))
    pool = multiprocessing.Pool(num_process)

    # log.info('samtools index...')
    # ret_val_1 = pool.apply_async(samtools_index, 
    #                             (nodup_bam, args.out_dir))
    # log.info('samtools flagstat...')
    # ret_val_2 = pool.apply_async(samtools_flagstat,
    #                             (nodup_bam, args.out_dir))
    log.info('sambamba index...')
    ret_val_1 = pool.apply_async(sambamba_index, 
                                (nodup_bam, args.nth, args.out_dir))
    log.info('sambamba flagstat...')
    ret_val_2 = pool.apply_async(sambamba_flagstat,
                                (nodup_bam, args.nth, args.out_dir))

    log.info('Generating PBC QC log...')
    if not args.no_dup_removal:
        if args.paired_end:
            ret_val_3 = pool.apply_async(pbc_qc_pe,
                            (filt_bam, nodup_bam, args.mito_chr_name,
                                max(1,args.nth-2),
                                args.out_dir))
        else:
            ret_val_3 = pool.apply_async(pbc_qc_se,
                            (filt_bam, nodup_bam, args.mito_chr_name, args.out_dir))
            
    # gather
    nodup_bai = ret_val_1.get(BIG_INT)
    nodup_flagstat_qc = ret_val_2.get(BIG_INT)

    if not args.no_dup_removal:
        pbc_qc = ret_val_3.get(BIG_INT)

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
