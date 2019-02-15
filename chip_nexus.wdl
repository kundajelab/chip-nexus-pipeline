# ENCODE DCC TF/Histone ChIP-Seq pipeline
# Author: Jin Lee (leepc12@gmail.com)

workflow chip_nexus {
	String pipeline_ver = 'v1.1.6'
	### sample name, description
	String title = 'Untitled'
	String description = 'No description'

	### mandatory genome param
	File genome_tsv 			# reference genome data TSV file including
								# all important genome specific data file paths and parameters
	### optional but important
	String aligner = 'bowtie' 		# bowtie or bwa
	Boolean align_only = false 		# disable all post-align analysis (peak-calling, overlap, idr, ...)
	Boolean true_rep_only = false 	# disable all analyses for pseudo replicates
									# overlap and idr will also be disabled
	Boolean enable_count_signal_track = false 		# generate count signal track
	Int cutadapt_min_overlap = 4	# minimum overlap for cutadapt -O
	Int cutadapt_min_trim_len = 22	# minimum trim length for cutadapt -m
	Float cutadapt_err_rate = 0.2	# Maximum allowed adapter error rate for cutadapt -e	
	Int nimnexus_trim = 0			# Pre-trim all reads by this length before processing (nimnexus -t)
	Int nimnexus_keep = 18			# Minimum number of bases required after barcode to keep read (nimnexus -k)
	Int nimnexus_randombarcode = 5	# Number of bases at the start of each read used for random barcode (nimnexus -r)

	String bowtie_param = '--chunkmbs 512 -k 1 -m 1 -v 2 --best --strata' # Parameters for bowtie command line
	Int mapq_thresh = 30			# threshold for low MAPQ reads removal
	Boolean no_dup_removal = false	# no dupe reads removal when filtering BAM
									# dup.qc and pbc.qc will be empty files
									# and nodup_bam in the output is 
									# filtered bam with dupes

	String mito_chr_name = 'chrM' 	# name of mito chromosome. THIS IS NOT A REG-EX! you can define only one chromosome name for mito.
	String regex_filter_reads = 'chrM' 	# Perl-style regular expression pattern for chr name to filter out reads
                        			# those reads will be excluded from peak calling
	Int subsample_reads = 0			# number of reads to subsample TAGALIGN
									# 0 for no subsampling. this affects all downstream analysis
	Boolean keep_irregular_chr_in_bfilt_peak = false # when filtering with blacklist
								# do not filter peaks with irregular chr name
								# and just keep them in bfilt_peak file
								# (e.g. keep chr1_AABBCC, AABR07024382.1, ...)
								# reg-ex pattern for regular chr names: /chr[\dXY]+[ \t]/

	### task-specific variables but defined in workflow level (limit of WDL)
	Int macs2_cap_num_peak = 500000	# cap number of raw peaks called from MACS2
	Float pval_thresh = 0.01		# p.value threshold
	Float idr_thresh = 0.05			# IDR threshold

	### resources
	Int trim_adapter_cpu = 2
	Int trim_adapter_mem_mb = 12000
	Int trim_adapter_time_hr = 24
	String trim_adapter_disks = "local-disk 100 HDD"

	Int bwa_cpu = 4
	Int bwa_mem_mb = 20000
	Int bwa_time_hr = 48
	String bwa_disks = "local-disk 100 HDD"

	Int bowtie_cpu = 4
	Int bowtie_mem_mb = 20000
	Int bowtie_time_hr = 48
	String bowtie_disks = "local-disk 100 HDD"

	Int filter_cpu = 2
	Int filter_mem_mb = 20000
	Int filter_time_hr = 24
	String filter_disks = "local-disk 100 HDD"

	Int bam2ta_cpu = 2
	Int bam2ta_mem_mb = 10000
	Int bam2ta_time_hr = 6
	String bam2ta_disks = "local-disk 100 HDD"

	Int spr_mem_mb = 16000

	Int fingerprint_cpu = 2
	Int fingerprint_mem_mb = 12000
	Int fingerprint_time_hr = 6
	String fingerprint_disks = "local-disk 100 HDD"

	Int macs2_mem_mb = 16000
	Int macs2_time_hr = 24
	String macs2_disks = "local-disk 100 HDD"

	#### input file definition
	# pipeline can start from any type of inputs and then leave all other types undefined
	# supported types: fastq, bam, nodup_bam (filtered bam), ta (tagAlign), peak
	# define up to 4 replicates
	# [rep_id] is for each replicate

 	### fastqs
 	# define fastqs either with DNANexus style (1-dim array) or with default one (3-dim array)
 	# [merge_id] is for pooing fastqs
 	## DNANexus UI style fastq definition
	Array[File] fastqs_rep1_R1 = []	# [merge_id]
	Array[File] fastqs_rep1_R2 = [] # do not define _R2 array if your sample is not paired end
	Array[File] fastqs_rep2_R1 = [] # do not define if you have a single replicate
	Array[File] fastqs_rep2_R2 = []	# do not define _R2 array if your sample is not paired end
	Array[File] fastqs_rep3_R1 = [] # do not define if you have <=2 replicates
	Array[File] fastqs_rep3_R2 = []	# do not define _R2 array if your sample is not paired end
	Array[File] fastqs_rep4_R1 = [] # do not define if you have <=3 replicates
	Array[File] fastqs_rep4_R2 = []	# do not define _R2 array if your sample is not paired end
	Array[File] fastqs_rep5_R1 = [] # do not define if you have <=4 replicates
	Array[File] fastqs_rep5_R2 = []	# do not define _R2 array if your sample is not paired end
	Array[File] fastqs_rep6_R1 = [] # do not define if you have <=5 replicates
	Array[File] fastqs_rep6_R2 = []	# do not define _R2 array if your sample is not paired end

	String adapter # adapter sequence
	String? barcodes # comma-delimited barcodes

 	## default style fastq definition
 	# [read_end_id] is for fastq R1 or fastq R2
	Array[Array[Array[File]]] fastqs = [] 		# [rep_id][merge_id][read_end_id]

	### other input types (bam, nodup_bam, ta)
	Array[File] bams = [] 			# [rep_id]
	Array[File] nodup_bams = [] 	# [rep_id]
	Array[File] tas = []			# [rep_id]

	### other input types (peak)
	Array[File] peaks = []		# [PAIR(rep_id1,rep_id2)]. example for 3 reps: [rep1_rep2, rep1_rep3, rep2_rep3]
	Array[File] peaks_pr1 = []	# [rep_id]. do not define if true_rep=true
	Array[File] peaks_pr2 = []	# [rep_id]. do not define if true_rep=true
	File? peak_ppr1				# do not define if you have a single replicate or true_rep=true
	File? peak_ppr2				# do not define if you have a single replicate or true_rep=true
	File? peak_pooled			# do not define if you have a single replicate or true_rep=true

	### other inputs used for resuming pipelines (QC/txt/log/png files, ...)
	File? ta_pooled
	Array[File] flagstat_qcs = []
	Array[File] pbc_qcs = []
	Array[File] nodup_flagstat_qcs = []
	Array[File] sig_pvals = []

	Array[File] macs2_frip_qcs = []
	Array[File] macs2_pr1_frip_qcs = []
	Array[File] macs2_pr2_frip_qcs = []
	File? macs2_pooled_frip_qc_
	File? macs2_ppr1_frip_qc_
	File? macs2_ppr2_frip_qc_

	Array[File] count_signal_track_pos_bws = []
	Array[File] count_signal_track_neg_bws = []
	File? count_signal_track_pooled_pos_bw_ 
	File? count_signal_track_pooled_neg_bw_ 

	### temp vars (do not define these)
	Boolean enable_idr = true
	String peak_type = 'narrowPeak' # peak type for IDR and overlap
	String idr_rank = 'p.value'

	### read genome data and paths
	call read_genome_tsv { input:genome_tsv = genome_tsv }
	File bwa_idx_tar = read_genome_tsv.genome['bwa_idx_tar']
	File bowtie_idx_tar = read_genome_tsv.genome['bowtie_idx_tar']
	File blacklist = read_genome_tsv.genome['blacklist']
	File chrsz = read_genome_tsv.genome['chrsz']
	String gensz = read_genome_tsv.genome['gensz']

	### pipeline starts here
	# temporary 2-dim arrays for DNANexus style fastqs
	Array[Array[File]] fastqs_rep1 = if length(fastqs_rep1_R2)>0 then transpose([fastqs_rep1_R1,fastqs_rep1_R2])
									else transpose([fastqs_rep1_R1])
	Array[Array[File]] fastqs_rep2 = if length(fastqs_rep2_R2)>0 then transpose([fastqs_rep2_R1,fastqs_rep2_R2])
									else transpose([fastqs_rep2_R1])
	Array[Array[File]] fastqs_rep3 = if length(fastqs_rep3_R2)>0 then transpose([fastqs_rep3_R1,fastqs_rep3_R2])
									else transpose([fastqs_rep3_R1])
	Array[Array[File]] fastqs_rep4 = if length(fastqs_rep4_R2)>0 then transpose([fastqs_rep4_R1,fastqs_rep4_R2])
									else transpose([fastqs_rep4_R1])
	Array[Array[File]] fastqs_rep5 = if length(fastqs_rep5_R2)>0 then transpose([fastqs_rep5_R1,fastqs_rep5_R2])
									else transpose([fastqs_rep5_R1])
	Array[Array[File]] fastqs_rep6 = if length(fastqs_rep6_R2)>0 then transpose([fastqs_rep6_R1,fastqs_rep6_R2])
									else transpose([fastqs_rep6_R1])
	Array[Array[Array[File]]] fastqs_ = if length(fastqs_rep1)<1 then fastqs
		else if length(fastqs_rep2)<1 then [fastqs_rep1]
		else if length(fastqs_rep3)<1 then [fastqs_rep1,fastqs_rep2]
		else if length(fastqs_rep4)<1 then [fastqs_rep1,fastqs_rep2,fastqs_rep3]
		else if length(fastqs_rep5)<1 then [fastqs_rep1,fastqs_rep2,fastqs_rep3,fastqs_rep4]
		else if length(fastqs_rep6)<1 then [fastqs_rep1,fastqs_rep2,fastqs_rep3,fastqs_rep4,fastqs_rep5]
		else [fastqs_rep1,fastqs_rep2,fastqs_rep3,fastqs_rep4,fastqs_rep5,fastqs_rep6]

	## temp vars for resuming pipelines
	Boolean need_to_process_ta = length(peaks_pr1)==0 && length(peaks)==0
	Boolean need_to_process_nodup_bam = need_to_process_ta && length(tas)==0
	Boolean need_to_process_bam = need_to_process_nodup_bam && length(nodup_bams)==0
	Boolean need_to_process_fastq = need_to_process_bam && length(bams)==0

	scatter(fastq_set in if need_to_process_fastq then fastqs_ else []) {
		# merge fastqs
		call trim_adapter { input :
			fastqs = fastq_set,
			adapter = adapter,
			barcodes = barcodes,
			trim = nimnexus_trim,
			keep = nimnexus_keep,
			randombarcode = nimnexus_randombarcode,
			min_trim_len = cutadapt_min_trim_len,
			err_rate = cutadapt_err_rate,
			min_overlap = cutadapt_min_overlap,

			cpu = trim_adapter_cpu,
			mem_mb = trim_adapter_mem_mb,
			time_hr = trim_adapter_time_hr,
			disks = trim_adapter_disks,
		}
		if ( aligner=='bwa' ) { # align merged fastqs with bwa		
			call bwa { input :
				idx_tar = bwa_idx_tar,
				fastqs = trim_adapter.trimmed_merged_fastqs, #[R1,R2]
				paired_end = false,
				cpu = bwa_cpu,
				mem_mb = bwa_mem_mb,
				time_hr = bwa_time_hr,
				disks = bwa_disks,
			}
		}
		if ( aligner=='bowtie' ) { # align merged fastqs with bowtie
			call bowtie { input :
				idx_tar = bowtie_idx_tar,
				fastqs = trim_adapter.trimmed_merged_fastqs, #[R1,R2]
				paired_end = false,
				param = bowtie_param,
				cpu = bowtie_cpu,
				mem_mb = bowtie_mem_mb,
				time_hr = bowtie_time_hr,
				disks = bowtie_disks,
			}
		}
	}

	Array[File?] bams_ = if (aligner=='bwa') then flatten([bwa.bam, bams])
		else flatten([bowtie.bam, bams])
	scatter(bam in if need_to_process_bam then bams_ else []) {
		# filter/dedup bam		
		call filter { input :
			bam = bam,
			paired_end = false,
			mapq_thresh = mapq_thresh,
			no_dup_removal = no_dup_removal,
			mito_chr_name = mito_chr_name,

			cpu = filter_cpu,
			mem_mb = filter_mem_mb,
			time_hr = filter_time_hr,
			disks = filter_disks,
		}
	}

	Array[File] nodup_bams_ = flatten([filter.nodup_bam, nodup_bams])
	scatter(bam in if need_to_process_nodup_bam then nodup_bams_ else []) {
		# convert bam to tagalign and subsample it if necessary
		call bam2ta { input :
			bam = bam,
			paired_end = false,
			subsample = subsample_reads,
			regex_grep_v_ta = regex_filter_reads,
			mito_chr_name = mito_chr_name,

			cpu = bam2ta_cpu,
			mem_mb = bam2ta_mem_mb,
			time_hr = bam2ta_time_hr,
			disks = bam2ta_disks,
		}
	}

	Array[File] tas_ = if align_only then [] else flatten([bam2ta.ta, tas])	
	Array[File] tas__ = if need_to_process_ta then tas_ else []
	if ( length(tas__)>1 )  {
		# pool tagaligns from true replicates
		call pool_ta { input :
			tas = tas__,
		}
	}

	if ( !true_rep_only ) {
		scatter( ta in tas__ ) {
			# make two self pseudo replicates per true replicate
			call spr { input :
				ta = ta,
				paired_end = false,
				mem_mb = spr_mem_mb,
			}
		}
	}
	if ( !true_rep_only && length(tas__)>1 ) {
		# pool tagaligns from pseudo replicates
		call pool_ta as pool_ta_pr1 { input :
			tas = spr.ta_pr1,
		}
		call pool_ta as pool_ta_pr2 { input :
			tas = spr.ta_pr2,
		}
	}

	# generate count signal track
	Array[File] tas_count_signal_track = if length(count_signal_track_pos_bws)>0 then []
		else if enable_count_signal_track then tas_
		else []
	scatter(i in range(length(tas_count_signal_track))) {
		call count_signal_track { input :
			ta = tas_count_signal_track[i],
			chrsz = chrsz,
		}
	}

	if ( !defined(count_signal_track_pooled_pos_bw_) && enable_count_signal_track && length(tas_)>0 ) {
		call count_signal_track as count_signal_track_pooled { input :
			ta = select_first([pool_ta.ta_pooled, ta_pooled]),
			chrsz = chrsz,
		}
	}

	# call peaks
	scatter(i in range(length(tas__))) {
		# always call MACS2 peaks for true replicates to get signal tracks
		# call peaks on tagalign
		call macs2 { input :
			ta = tas__[i],
			gensz = gensz,
			chrsz = chrsz,
			cap_num_peak = macs2_cap_num_peak,
			pval_thresh = pval_thresh,
			make_signal = true,
			blacklist = blacklist,
			keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,

			mem_mb = macs2_mem_mb,
			disks = macs2_disks,
			time_hr = macs2_time_hr,
		}
	}

	if ( !true_rep_only ) {
		scatter(i in range(length(tas__))) {
			# call peaks on 1st pseudo replicated tagalign 
			call macs2 as macs2_pr1 { input :
				ta = select_first([spr.ta_pr1])[i],
				gensz = gensz,
				chrsz = chrsz,
				cap_num_peak = macs2_cap_num_peak,
				pval_thresh = pval_thresh,
				blacklist = blacklist,
				make_signal = false,
				keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
	
				mem_mb = macs2_mem_mb,
				disks = macs2_disks,
				time_hr = macs2_time_hr,
			}
			call macs2 as macs2_pr2 { input :
				ta = select_first([spr.ta_pr2])[i],
				gensz = gensz,
				chrsz = chrsz,
				cap_num_peak = macs2_cap_num_peak,
				pval_thresh = pval_thresh,
				blacklist = blacklist,
				make_signal = false,
				keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
	
				mem_mb = macs2_mem_mb,
				disks = macs2_disks,
				time_hr = macs2_time_hr,
			}					
		}
	}

	if ( length(tas__)>1 ) {
		# call peaks on pooled replicate
		# always call MACS2 peaks for pooled replicate to get signal tracks
		call macs2 as macs2_pooled { input :
			ta = select_first([pool_ta.ta_pooled]),
			gensz = gensz,
			chrsz = chrsz,
			cap_num_peak = macs2_cap_num_peak,
			pval_thresh = pval_thresh,
			make_signal = true,
			blacklist = blacklist,
			keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,

			mem_mb = macs2_mem_mb,
			disks = macs2_disks,
			time_hr = macs2_time_hr,
		}
	}
	if ( !true_rep_only && length(tas__)>1 ) {
		# call peaks on 1st pooled pseudo replicates
		call macs2 as macs2_ppr1 { input :
			ta = select_first([pool_ta_pr1.ta_pooled]),
			gensz = gensz,
			chrsz = chrsz,
			cap_num_peak = macs2_cap_num_peak,
			pval_thresh = pval_thresh,
			blacklist = blacklist,
			make_signal = false,
			keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,

			mem_mb = macs2_mem_mb,
			disks = macs2_disks,
			time_hr = macs2_time_hr,
		}
	}
	if ( !true_rep_only && length(tas__)>1 ) {
		# call peaks on 2nd pooled pseudo replicates
		call macs2 as macs2_ppr2 { input :
			ta = select_first([pool_ta_pr2.ta_pooled]),
			gensz = gensz,
			chrsz = chrsz,
			cap_num_peak = macs2_cap_num_peak,
			pval_thresh = pval_thresh,
			blacklist = blacklist,
			make_signal = false,
			keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,

			mem_mb = macs2_mem_mb,
			disks = macs2_disks,
			time_hr = macs2_time_hr,
		}
	}

	# make peak arrays
	Array[File] peaks_ = if align_only then [] else flatten([macs2.npeak, peaks])

	# generate all possible pairs of true replicates (pair: left=prefix, right=[peak1,peak2])
	Array[Pair[String,Array[File]]] peak_pairs =  
		if length(peaks_)<=1 then [] # 1 rep
		else if length(peaks_)<=2 then # 2 reps
			 [('rep1-rep2',[peaks_[0],peaks_[1]])]
		else if length(peaks_)<=3 then # 3 reps
			 [('rep1-rep2',[peaks_[0],peaks_[1]]), ('rep1-rep3',[peaks_[0],peaks_[2]]),
			  ('rep2-rep3',[peaks_[1],peaks_[2]])]
		else if length(peaks_)<=4 then # 4 reps
			 [('rep1-rep2',[peaks_[0],peaks_[1]]), ('rep1-rep3',[peaks_[0],peaks_[2]]), ('rep1-rep4',[peaks_[0],peaks_[3]]),
			  ('rep2-rep3',[peaks_[1],peaks_[2]]), ('rep2-rep4',[peaks_[1],peaks_[3]]),
			  ('rep3-rep4',[peaks_[2],peaks_[3]])]
		else if length(peaks_)<=5 then # 5 reps
			 [('rep1-rep2',[peaks_[0],peaks_[1]]), ('rep1-rep3',[peaks_[0],peaks_[2]]), ('rep1-rep4',[peaks_[0],peaks_[3]]), ('rep1-rep5',[peaks_[0],peaks_[4]]),
			  ('rep2-rep3',[peaks_[1],peaks_[2]]), ('rep2-rep4',[peaks_[1],peaks_[3]]), ('rep2-rep5',[peaks_[1],peaks_[4]]),
			  ('rep3-rep4',[peaks_[2],peaks_[3]]), ('rep3-rep5',[peaks_[2],peaks_[4]]),
			  ('rep4-rep5',[peaks_[3],peaks_[4]])]
		else # 6 reps
			 [('rep1-rep2',[peaks_[0],peaks_[1]]), ('rep1-rep3',[peaks_[0],peaks_[2]]), ('rep1-rep4',[peaks_[0],peaks_[3]]), ('rep1-rep5',[peaks_[0],peaks_[4]]), ('rep1-rep6',[peaks_[0],peaks_[5]]),
			  ('rep2-rep3',[peaks_[1],peaks_[2]]), ('rep2-rep4',[peaks_[1],peaks_[3]]), ('rep2-rep5',[peaks_[1],peaks_[4]]), ('rep2-rep6',[peaks_[1],peaks_[5]]),
			  ('rep3-rep4',[peaks_[2],peaks_[3]]), ('rep3-rep5',[peaks_[2],peaks_[4]]), ('rep3-rep6',[peaks_[2],peaks_[5]]),
			  ('rep4-rep5',[peaks_[3],peaks_[4]]), ('rep4-rep6',[peaks_[3],peaks_[5]]),
			  ('rep5-rep6',[peaks_[4],peaks_[5]])]
	if ( length(peaks_)>0 ) {
		scatter( pair in peak_pairs ) {
			# Naive overlap on every pair of true replicates
			call overlap { input :
				prefix = pair.left,
				peak1 = pair.right[0],
				peak2 = pair.right[1],
				peak_pooled = select_first([macs2_pooled.npeak, peak_pooled]),
				peak_type = peak_type,
				blacklist = blacklist,
				chrsz = chrsz,
				keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
				ta = if defined(ta_pooled) then ta_pooled else pool_ta.ta_pooled,
			}
		}
	}
	if ( length(peaks_)>0 && enable_idr ) {
		scatter( pair in peak_pairs ) {
			# IDR on every pair of true replicates
			call idr { input : 
				prefix = pair.left,
				peak1 = pair.right[0],
				peak2 = pair.right[1],
				peak_pooled = select_first([macs2_pooled.npeak, peak_pooled]),
				idr_thresh = idr_thresh,
				peak_type = peak_type,
				rank = idr_rank,
				blacklist = blacklist,
				chrsz = chrsz,
				keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
				ta = if defined(ta_pooled) then ta_pooled else pool_ta.ta_pooled,
			}
		}
	}

	Array[File] peaks_pr1_ = flatten(select_all([macs2_pr1.npeak, peaks_pr1]))
	Array[File] peaks_pr2_ = flatten(select_all([macs2_pr2.npeak, peaks_pr2]))

	scatter( i in range(length(peaks_pr1_)) ) {
		# Naive overlap on pseduo replicates
		call overlap as overlap_pr { input : 
			prefix = "rep"+(i+1)+"-pr",
			peak1 = peaks_pr1_[i],
			peak2 = peaks_pr2_[i],
			peak_pooled = peaks_[i],
			peak_type = peak_type,
			blacklist = blacklist,
			chrsz = chrsz,
			keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
			ta = if length(tas_)>0 then tas_[i] else if defined(ta_pooled) then ta_pooled else pool_ta.ta_pooled,
		}
	}
	if ( enable_idr ) {
		scatter( i in range(length(peaks_pr1_)) ) {
			# IDR on pseduo replicates
			call idr as idr_pr { input : 
				prefix = "rep"+(i+1)+"-pr",
				peak1 = peaks_pr1_[i],
				peak2 = peaks_pr2_[i],
				peak_pooled = peaks_[i],
				idr_thresh = idr_thresh,
				peak_type = peak_type,
				rank = idr_rank,
				blacklist = blacklist,
				chrsz = chrsz,
				keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
				ta = if length(tas_)>0 then tas_[i] else if defined(ta_pooled) then ta_pooled else pool_ta.ta_pooled,
			}
		}
	}
	if ( length(peaks_pr1_)>1 ) {
		# Naive overlap on pooled pseudo replicates
		call overlap as overlap_ppr { input : 
			prefix = "ppr",
			peak1 = select_first([macs2_ppr1.npeak, peak_ppr1]), #peak_ppr1_[0],
			peak2 = select_first([macs2_ppr2.npeak, peak_ppr2]), #peak_ppr2_[0],
			peak_pooled = select_first([macs2_pooled.npeak, peak_pooled]),
			peak_type = peak_type,
			blacklist = blacklist,
			chrsz = chrsz,
			keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
			ta = if defined(ta_pooled) then ta_pooled else pool_ta.ta_pooled,
		}
	}
	if ( enable_idr && length(peaks_pr1_)>1 ) {
		# IDR on pooled pseduo replicates
		call idr as idr_ppr { input : 
			prefix = "ppr",
			peak1 = select_first([macs2_ppr1.npeak, peak_ppr1]), #peak_ppr1_[0],
			peak2 = select_first([macs2_ppr2.npeak, peak_ppr2]), #peak_ppr2_[0],
			peak_pooled = select_first([macs2_pooled.npeak, peak_pooled]),
			idr_thresh = idr_thresh,
			peak_type = peak_type,
			rank = idr_rank,
			blacklist = blacklist,
			chrsz = chrsz,
			keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
			ta = if defined(ta_pooled) then ta_pooled else pool_ta.ta_pooled,
		}
	}

	if ( !align_only && !true_rep_only ) {
		# reproducibility QC for overlapping peaks
		call reproducibility as reproducibility_overlap { input :
			prefix = 'overlap',
			peaks = overlap.bfilt_overlap_peak,
			peaks_pr = overlap_pr.bfilt_overlap_peak,
			peak_ppr = overlap_ppr.bfilt_overlap_peak,
			peak_type = peak_type,
			chrsz = chrsz,
			keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
		}
	}
	if ( !align_only && !true_rep_only && enable_idr ) {
		# reproducibility QC for IDR peaks
		call reproducibility as reproducibility_idr { input :
			prefix = 'idr',
			peaks = idr.bfilt_idr_peak,
			peaks_pr = idr_pr.bfilt_idr_peak,
			peak_ppr = idr_ppr.bfilt_idr_peak,
			peak_type = peak_type,
			chrsz = chrsz,
			keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
		}
	}

	#Array[File?] flagstat_qcs_ = flatten([flagstat_qcs, bwa.flagstat_qc, bowtie.flagstat_qc])
	Array[File?] flagstat_qcs_ = if (aligner=='bwa') then flatten([flagstat_qcs, bwa.flagstat_qc])
		else flatten([flagstat_qcs, bowtie.flagstat_qc])

	Array[File] pbc_qcs_ = flatten([pbc_qcs, filter.pbc_qc])
	Array[File] nodup_flagstat_qcs_ = flatten([nodup_flagstat_qcs, filter.flagstat_qc])

	Array[File] sig_pvals_ = flatten([sig_pvals, macs2.sig_pval])
	
	Array[File] macs2_frip_qcs_ = flatten([macs2_frip_qcs, macs2.frip_qc])
	Array[File] macs2_pr1_frip_qcs_ = flatten(select_all([macs2_pr1_frip_qcs, macs2_pr1.frip_qc]))
	Array[File] macs2_pr2_frip_qcs_ = flatten(select_all([macs2_pr2_frip_qcs, macs2_pr2.frip_qc]))

	# Generate final QC report and JSON
	call qc_report { input :
		pipeline_ver = pipeline_ver,
		title = title,
		description = description,
		genome = basename(genome_tsv),
		paired_end = false,
		pipeline_type = 'ChIP-nexus',
		peak_caller = 'macs2',
		macs2_cap_num_peak = macs2_cap_num_peak,
		idr_thresh = idr_thresh,

		flagstat_qcs = flagstat_qcs_,
		nodup_flagstat_qcs = nodup_flagstat_qcs_,
		pbc_qcs = pbc_qcs_,

		frip_macs2_qcs = macs2_frip_qcs_,
		frip_macs2_qcs_pr1 = macs2_pr1_frip_qcs_,
		frip_macs2_qcs_pr2 = macs2_pr2_frip_qcs_,
		frip_macs2_qc_pooled = if defined(macs2_pooled_frip_qc_) then macs2_pooled_frip_qc_ else macs2_pooled.frip_qc,
		frip_macs2_qc_ppr1 = if defined(macs2_ppr1_frip_qc_) then macs2_ppr1_frip_qc_ else macs2_ppr1.frip_qc,
		frip_macs2_qc_ppr2 = if defined(macs2_ppr2_frip_qc_) then macs2_ppr2_frip_qc_ else macs2_ppr2.frip_qc,

		idr_plots = idr.idr_plot,
		idr_plots_pr = idr_pr.idr_plot,
		idr_plot_ppr = idr_ppr.idr_plot,
		frip_idr_qcs = idr.frip_qc,
		frip_idr_qcs_pr = idr_pr.frip_qc,
		frip_idr_qc_ppr = idr_ppr.frip_qc,
		frip_overlap_qcs = overlap.frip_qc,
		frip_overlap_qcs_pr = overlap_pr.frip_qc,
		frip_overlap_qc_ppr = overlap_ppr.frip_qc,
		idr_reproducibility_qc = reproducibility_idr.reproducibility_qc,
		overlap_reproducibility_qc = reproducibility_overlap.reproducibility_qc,
	}

	output {
		File report = qc_report.report
		File qc_json = qc_report.qc_json
		Boolean qc_json_ref_match = qc_report.qc_json_ref_match
	}	
}

task trim_adapter { # trim adapters and merge trimmed fastqs
	Array[Array[File]] fastqs 		# [merge_id][read_end_id]
	String adapter
	String? barcodes
	# optional
	Int min_overlap 		# minimum overlap for cutadapt -O
	Int min_trim_len 		# minimum trim length for cutadapt -m
	Float err_rate			# Maximum allowed adapter error rate 
							# for cutadapt -e
	Int trim				# Pre-trim all reads by this length before processing (nimnexus -t)
	Int keep				# Minimum number of bases required after barcode to keep read (nimnexus -k)
	Int randombarcode		# Number of bases at the start of each read used for random barcode (nimnexus -r)

	Int cpu
	Int mem_mb
	Int time_hr
	String disks

	command {
		python $(which encode_trim_adapter_chip_nexus.py) \
			${write_tsv(fastqs)} \
			${"--adapter " + adapter} \
			${"--barcodes " + barcodes} \
			${"--min-overlap " + min_overlap} \
			${"--min-trim-len " + min_trim_len} \
			${"--err-rate " + err_rate} \
			${"--trim " + trim} \
			${"--keep " + keep} \
			${"--randombarcode " + randombarcode} \
			${"--nth " + cpu}
	}
	output {
		# WDL glob() globs in an alphabetical order
		# so R1 and R2 can be switched, which results in an
		# unexpected behavior of a workflow
		# so we prepend merge_fastqs_'end'_ (R1 or R2)
		# to the basename of original filename
		# this prefix will be later stripped in bwa/bowtie task
		Array[File] trimmed_merged_fastqs = glob("merge_fastqs_R?_*.fastq.gz")
	}
	runtime {
		cpu : cpu
		memory : "${mem_mb} MB"
		time : time_hr
		disks : disks
	}
}

task bwa {
	File idx_tar 		# reference bwa index tar
	Array[File] fastqs 	# [read_end_id]
	Boolean paired_end

	Int cpu
	Int mem_mb
	Int time_hr
	String disks

	command {
		python $(which encode_bwa.py) \
			${idx_tar} \
			${sep=' ' fastqs} \
			${if paired_end then "--paired-end" else ""} \
			${"--nth " + cpu}
	}
	output {
		File bam = glob("*.bam")[0]
		File bai = glob("*.bai")[0]
		File flagstat_qc = glob("*.flagstat.qc")[0]
	}
	runtime {
		cpu : cpu
		memory : "${mem_mb} MB"
		time : time_hr
		disks : disks
		preemptible: 0
	}
}

task bowtie {
	File idx_tar 		# reference bowtie index tar
	Array[File] fastqs 	# [read_end_id]
	String param
	Boolean paired_end

	Int cpu
	Int mem_mb
	Int time_hr
	String disks

	command {
		python $(which encode_bowtie.py) \
			${idx_tar} \
			${sep=' ' fastqs} \
			${if paired_end then "--paired-end" else ""} \
			${"--param ' " + param + "'"} \
			${"--nth " + cpu}
	}
	output {
		File bam = glob("*.bam")[0]
		File bai = glob("*.bai")[0]
		File flagstat_qc = glob("*.flagstat.qc")[0]
	}
	runtime {
		cpu : cpu
		memory : "${mem_mb} MB"
		time : time_hr
		disks : disks
		preemptible: 0
	}
}

task filter {
	File bam
	Boolean paired_end
	Int mapq_thresh				# threshold for low MAPQ reads removal
	Boolean no_dup_removal 		# no dupe reads removal when filtering BAM
									# dup.qc and pbc.qc will be empty files
									# and nodup_bam in the output is 
									# filtered bam with dupes	
	String mito_chr_name
	Int cpu
	Int mem_mb
	Int time_hr
	String disks

	command {
		${if no_dup_removal then "touch null.dup.qc null.pbc.qc; " else ""}
		touch null
		python $(which encode_filter_chip_nexus.py) \
			${bam} \
			${if paired_end then "--paired-end" else ""} \
			${"--multimapping " + 0} \
			${"--mapq-thresh " + mapq_thresh} \
			${if no_dup_removal then "--no-dup-removal" else ""} \
			${"--mito-chr-name " + mito_chr_name} \
			${"--nth " + cpu}
	}
	output {
		File nodup_bam = glob("*.bam")[0]
		File nodup_bai = glob("*.bai")[0]
		File flagstat_qc = glob("*.flagstat.qc")[0]
		File pbc_qc = if no_dup_removal then glob("null")[0] else glob("*.pbc.qc")[0]
	}
	runtime {
		cpu : cpu
		memory : "${mem_mb} MB"
		time : time_hr
		disks : disks
	}
}

task bam2ta {
	File bam
	Boolean paired_end
	String regex_grep_v_ta   	# Perl-style regular expression pattern 
                        		# to remove matching reads from TAGALIGN
	String mito_chr_name 		# mito chromosome name
	Int subsample 				# number of reads to subsample TAGALIGN
								# this affects all downstream analysis
	Int cpu
	Int mem_mb
	Int time_hr
	String disks

	command {
		python $(which encode_bam2ta.py) \
			${bam} \
			--disable-tn5-shift \
			${if paired_end then "--paired-end" else ""} \
			${if regex_grep_v_ta!="" then "--regex-grep-v-ta '"+regex_grep_v_ta+"'" else ""} \
			${"--mito-chr-name " + mito_chr_name} \
			${"--subsample " + subsample} \
			${"--nth " + cpu}
	}
	output {
		File ta = glob("*.tagAlign.gz")[0]
	}
	runtime {
		cpu : cpu
		memory : "${mem_mb} MB"
		time : time_hr
		disks : disks
	}
}

task spr { # make two self pseudo replicates
	File ta
	Boolean paired_end

	Int mem_mb

	command {
		python $(which encode_spr.py) \
			${ta} \
			${if paired_end then "--paired-end" else ""}
	}
	output {
		File ta_pr1 = glob("*.pr1.tagAlign.gz")[0]
		File ta_pr2 = glob("*.pr2.tagAlign.gz")[0]
	}
	runtime {
		cpu : 1
		memory : "${mem_mb} MB"
		time : 1
		disks : "local-disk 50 HDD"
	}
}

task pool_ta {
	Array[File] tas

	command {
		python $(which encode_pool_ta.py) \
			${sep=' ' tas}
	}
	output {
		File ta_pooled = glob("*.tagAlign.gz")[0]
	}
	runtime {
		cpu : 1
		memory : "4000 MB"
		time : 1
		disks : "local-disk 50 HDD"
	}
}

task count_signal_track {
	File ta 			# tag-align
	File chrsz			# 2-col chromosome sizes file

	command {
		python $(which encode_count_signal_track.py) \
			${ta} \
			${"--chrsz " + chrsz}
	}
	output {
		File pos_bw = glob("*.positive.bigwig")[0]
		File neg_bw = glob("*.negative.bigwig")[0]
	}
	runtime {
		cpu : 1
		memory : "8000 MB"
		time : 4
		disks : "local-disk 50 HDD"
	}
}

task macs2 {
	File ta 			# tag-align
	String gensz		# Genome size (sum of entries in 2nd column of 
                        # chr. sizes file, or hs for human, ms for mouse)
	File chrsz			# 2-col chromosome sizes file
	Int cap_num_peak	# cap number of raw peaks called from MACS2
	Float pval_thresh 	# p.value threshold
	Boolean make_signal
	File blacklist 		# blacklist BED to filter raw peaks
	Boolean	keep_irregular_chr_in_bfilt_peak

	Int mem_mb
	Int time_hr
	String disks

	command {
		${if make_signal then "" 
			else "touch null.pval.signal.bigwig null.fc.signal.bigwig"}
		touch null
		python $(which encode_macs2_chip.py) \
			${ta} \
			${"--gensz "+ gensz} \
			${"--chrsz " + chrsz} \
			--fraglen 150 --shift -75 \
			${"--cap-num-peak " + cap_num_peak} \
			${"--pval-thresh "+ pval_thresh} \
			${if make_signal then "--make-signal" else ""} \
			${if keep_irregular_chr_in_bfilt_peak then "--keep-irregular-chr" else ""} \
			${"--blacklist "+ blacklist}
	}
	output {
		File npeak = glob("*[!.][!b][!f][!i][!l][!t].narrowPeak.gz")[0]
		File bfilt_npeak = glob("*.bfilt.narrowPeak.gz")[0]
		File bfilt_npeak_bb = glob("*.bfilt.narrowPeak.bb")[0]
		Array[File] bfilt_npeak_hammock = glob("*.bfilt.narrowPeak.hammock.gz*")
		File sig_pval = if make_signal then glob("*.pval.signal.bigwig")[0] else glob("null")[0]
		File sig_fc = if make_signal then glob("*.fc.signal.bigwig")[0] else glob("null")[0]
		File frip_qc = glob("*.frip.qc")[0]
	}
	runtime {
		cpu : 1
		memory : "${mem_mb} MB"
		time : time_hr
		disks : disks
	}
}

task idr {
	String prefix 		# prefix for IDR output file
	File peak1 			
	File peak2
	File peak_pooled
	Float idr_thresh
	File blacklist 		# blacklist BED to filter raw peaks
	Boolean	keep_irregular_chr_in_bfilt_peak
	# parameters to compute FRiP
	File? ta			# to calculate FRiP
	File chrsz			# 2-col chromosome sizes file
	String peak_type
	String rank

	command {
		${if defined(ta) then "" else "touch null.frip.qc"}			
		touch null 
		python $(which encode_idr.py) \
			${peak1} ${peak2} ${peak_pooled} \
			${"--prefix " + prefix} \
			${"--idr-thresh " + idr_thresh} \
			${"--peak-type " + peak_type} \
			--idr-rank ${rank} \
			--fraglen 150 \
			${"--chrsz " + chrsz} \
			${"--blacklist "+ blacklist} \
			${if keep_irregular_chr_in_bfilt_peak then "--keep-irregular-chr" else ""} \
			${"--ta " + ta}
	}
	output {
		File idr_peak = glob("*[!.][!b][!f][!i][!l][!t]."+peak_type+".gz")[0]
		File bfilt_idr_peak = glob("*.bfilt."+peak_type+".gz")[0]
		File bfilt_idr_peak_bb = glob("*.bfilt."+peak_type+".bb")[0]
		Array[File] bfilt_idr_peak_hammock = glob("*.bfilt."+peak_type+".hammock.gz*")
		File idr_plot = glob("*.txt.png")[0]
		File idr_unthresholded_peak = glob("*.txt.gz")[0]
		File idr_log = glob("*.log")[0]
		File frip_qc = if defined(ta) then glob("*.frip.qc")[0] else glob("null")[0]
	}
	runtime {
		cpu : 1
		memory : "4000 MB"
		time : 1
		disks : "local-disk 50 HDD"
	}	
}

task overlap {
	String prefix 		# prefix for IDR output file
	File peak1
	File peak2
	File peak_pooled
	File blacklist 	# blacklist BED to filter raw peaks
	Boolean	keep_irregular_chr_in_bfilt_peak
	# parameters to compute FRiP
	File? ta		# to calculate FRiP
	File chrsz		# 2-col chromosome sizes file
	String peak_type

	command {
		${if defined(ta) then "" else "touch null.frip.qc"}			
		touch null 
		python $(which encode_naive_overlap.py) \
			${peak1} ${peak2} ${peak_pooled} \
			${"--prefix " + prefix} \
			${"--peak-type " + peak_type} \
			--fraglen 150 \
			${"--chrsz " + chrsz} \
			${"--blacklist "+ blacklist} \
			--nonamecheck \
			${if keep_irregular_chr_in_bfilt_peak then "--keep-irregular-chr" else ""} \
			${"--ta " + ta}
	}
	output {
		File overlap_peak = glob("*[!.][!b][!f][!i][!l][!t]."+peak_type+".gz")[0]
		File bfilt_overlap_peak = glob("*.bfilt."+peak_type+".gz")[0]
		File bfilt_overlap_peak_bb = glob("*.bfilt."+peak_type+".bb")[0]
		Array[File] bfilt_overlap_peak_hammock = glob("*.bfilt."+peak_type+".hammock.gz*")
		File frip_qc = if defined(ta) then glob("*.frip.qc")[0] else glob("null")[0]
	}
	runtime {
		cpu : 1
		memory : "4000 MB"
		time : 1
		disks : "local-disk 50 HDD"
	}	
}

task reproducibility {
	String prefix
	Array[File] peaks # peak files from pair of true replicates
						# in a sorted order. for example of 4 replicates,
						# 1,2 1,3 1,4 2,3 2,4 3,4.
                        # x,y means peak file from rep-x vs rep-y
	Array[File] peaks_pr	# peak files from pseudo replicates
	File? peak_ppr			# Peak file from pooled pseudo replicate.
	String peak_type
	File chrsz			# 2-col chromosome sizes file
	Boolean	keep_irregular_chr_in_bfilt_peak

	command {
		python $(which encode_reproducibility_qc.py) \
			${sep=' ' peaks} \
			--peaks-pr ${sep=' ' peaks_pr} \
			${"--peak-ppr "+ peak_ppr} \
			--prefix ${prefix} \
			${"--peak-type " + peak_type} \
			${if keep_irregular_chr_in_bfilt_peak then "--keep-irregular-chr" else ""} \
			${"--chrsz " + chrsz}
	}
	output {
		File optimal_peak = glob("optimal_peak.*.gz")[0]
		File conservative_peak = glob("conservative_peak.*.gz")[0]
		File optimal_peak_bb = glob("optimal_peak.*.bb")[0]
		File conservative_peak_bb = glob("conservative_peak.*.bb")[0]
		Array[File] optimal_peak_hammock = glob("optimal_peak.*.hammock.gz*")
		Array[File] conservative_peak_hammock = glob("conservative_peak.*.hammock.gz*")
		File reproducibility_qc = glob("*reproducibility.qc")[0]
	}
	runtime {
		cpu : 1
		memory : "4000 MB"
		time : 1
		disks : "local-disk 50 HDD"
	}	
}

# gather all outputs and generate 
# - qc.html		: organized final HTML report
# - qc.json		: all QCs
task qc_report {
	# optional metadata
	String pipeline_ver
 	String title # name of sample
	String description # description for sample
	String? genome
	#String? encode_accession_id	# ENCODE accession ID of sample
	# workflow params
	Boolean paired_end
	String pipeline_type
	String peak_caller
	Int? macs2_cap_num_peak
	Float idr_thresh
	# QCs
	Array[File?]? flagstat_qcs
	Array[File]? nodup_flagstat_qcs
	Array[File]? pbc_qcs
	File? jsd_plot 
	Array[File]? jsd_qcs
	Array[File]? idr_plots
	Array[File]? idr_plots_pr
	File? idr_plot_ppr
	Array[File]? frip_macs2_qcs
	Array[File]? frip_macs2_qcs_pr1
	Array[File]? frip_macs2_qcs_pr2
	File? frip_macs2_qc_pooled
	File? frip_macs2_qc_ppr1 
	File? frip_macs2_qc_ppr2 
	Array[File]? frip_idr_qcs
	Array[File]? frip_idr_qcs_pr
	File? frip_idr_qc_ppr 
	Array[File]? frip_overlap_qcs
	Array[File]? frip_overlap_qcs_pr
	File? frip_overlap_qc_ppr
	File? idr_reproducibility_qc
	File? overlap_reproducibility_qc

	File? qc_json_ref

	command {
		python $(which encode_qc_report.py) \
			${"--pipeline-ver " + pipeline_ver} \
			${"--title '" + sub(title,"'","_") + "'"} \
			${"--desc '" + sub(description,"'","_") + "'"} \
			${"--genome " + genome} \
			${"--multimapping " + 0} \
			${if paired_end then "--paired-end" else ""} \
			--pipeline-type ${pipeline_type} \
			--peak-caller ${peak_caller} \
			${"--macs2-cap-num-peak " + macs2_cap_num_peak} \
			--idr-thresh ${idr_thresh} \
			--flagstat-qcs ${sep=' ' flagstat_qcs} \
			--nodup-flagstat-qcs ${sep=' ' nodup_flagstat_qcs} \
			--pbc-qcs ${sep=' ' pbc_qcs} \
			${"--jsd-plot " + jsd_plot} \
			--jsd-qcs ${sep=' ' jsd_qcs} \
			--idr-plots ${sep=' ' idr_plots} \
			--idr-plots-pr ${sep=' ' idr_plots_pr} \
			${"--idr-plot-ppr " + idr_plot_ppr} \
			--frip-macs2-qcs ${sep=' ' frip_macs2_qcs} \
			--frip-macs2-qcs-pr1 ${sep=' ' frip_macs2_qcs_pr1} \
			--frip-macs2-qcs-pr2 ${sep=' ' frip_macs2_qcs_pr2} \
			${"--frip-macs2-qc-pooled " + frip_macs2_qc_pooled} \
			${"--frip-macs2-qc-ppr1 " + frip_macs2_qc_ppr1} \
			${"--frip-macs2-qc-ppr2 " + frip_macs2_qc_ppr2} \
			--frip-idr-qcs ${sep=' ' frip_idr_qcs} \
			--frip-idr-qcs-pr ${sep=' ' frip_idr_qcs_pr} \
			${"--frip-idr-qc-ppr " + frip_idr_qc_ppr} \
			--frip-overlap-qcs ${sep=' ' frip_overlap_qcs} \
			--frip-overlap-qcs-pr ${sep=' ' frip_overlap_qcs_pr} \
			${"--frip-overlap-qc-ppr " + frip_overlap_qc_ppr} \
			${"--idr-reproducibility-qc " + idr_reproducibility_qc} \
			${"--overlap-reproducibility-qc " + overlap_reproducibility_qc} \
			--out-qc-html qc.html \
			--out-qc-json qc.json \
			${"--qc-json-ref " + qc_json_ref}
	}
	output {
		File report = glob('*qc.html')[0]
		File qc_json = glob('*qc.json')[0]
		Boolean qc_json_ref_match = read_string("qc_json_ref_match.txt")=="True"
	}
	runtime {
		cpu : 1
		memory : "4000 MB"
		time : 1
		disks : "local-disk 50 HDD"
	}	
}

### workflow system tasks

task read_genome_tsv {
	File genome_tsv
	command {
		cat ${genome_tsv} > 'tmp.tsv'
	}
	output {
		Map[String,String] genome = read_map('tmp.tsv')
	}
	runtime {
		cpu : 1
		memory : "4000 MB"
		time : 1
		disks : "local-disk 50 HDD"		
	}
}

task rounded_mean {
	Array[Int] ints
	command <<<
		python <<CODE
		arr = [${sep=',' ints}]
		with open('tmp.txt','w') as fp:
			if len(arr):
			    sum_ = sum(arr)
			    mean_ = sum(arr)/float(len(arr))
			    fp.write('{}'.format(int(round(mean_))))
			else:
			    fp.write('0')
		CODE
	>>>
	output {
		Int rounded_mean = read_int('tmp.txt')
	}
	runtime {
		cpu : 1
		memory : "4000 MB"
		time : 1
		disks : "local-disk 50 HDD"
	}	
}

task compare_md5sum {
	Array[String] labels
	Array[File] files
	Array[File] ref_files

	command <<<
		python <<CODE	
		from collections import OrderedDict
		import os
		import json
		import hashlib

		def md5sum(filename, blocksize=65536):
		    hash = hashlib.md5()
		    with open(filename, 'rb') as f:
		        for block in iter(lambda: f.read(blocksize), b""):
		            hash.update(block)
		    return hash.hexdigest()

		with open('${write_lines(labels)}','r') as fp:
			labels = fp.read().splitlines()
		with open('${write_lines(files)}','r') as fp:
			files = fp.read().splitlines()
		with open('${write_lines(ref_files)}','r') as fp:
			ref_files = fp.read().splitlines()

		result = OrderedDict()
		match = OrderedDict()
		match_overall = True

		result['tasks'] = []
		result['failed_task_labels'] = []
		result['succeeded_task_labels'] = []
		for i, label in enumerate(labels):
			f = files[i]
			ref_f = ref_files[i]
			md5 = md5sum(f)
			ref_md5 = md5sum(ref_f)
			# if text file, read in contents
			if f.endswith('.qc') or f.endswith('.txt') or \
				f.endswith('.log') or f.endswith('.out'):
				with open(f,'r') as fp:
					contents = fp.read()
				with open(ref_f,'r') as fp:
					ref_contents = fp.read()
			else:
				contents = ''
				ref_contents = ''
			matched = md5==ref_md5
			result['tasks'].append(OrderedDict([
				('label', label),
				('match', matched),
				('md5sum', md5),
				('ref_md5sum', ref_md5),
				('basename', os.path.basename(f)),
				('ref_basename', os.path.basename(ref_f)),
				('contents', contents),
				('ref_contents', ref_contents),
				]))
			match[label] = matched
			match_overall &= matched
			if matched:
				result['succeeded_task_labels'].append(label)
			else:
				result['failed_task_labels'].append(label)		
		result['match_overall'] = match_overall

		with open('result.json','w') as fp:
			fp.write(json.dumps(result, indent=4))
		match_tmp = []
		for key in match:
			val = match[key]
			match_tmp.append('{}\t{}'.format(key, val))
		with open('match.tsv','w') as fp:
			fp.writelines('\n'.join(match_tmp))
		with open('match_overall.txt','w') as fp:
			fp.write(str(match_overall))
		CODE
	>>>
	output {
		Map[String,String] match = read_map('match.tsv') # key:label, val:match
		Boolean match_overall = read_boolean('match_overall.txt')
		File json = glob('result.json')[0] # details (json file)
		String json_str = read_string('result.json') # details (string)
	}
	runtime {
		cpu : 1
		memory : "4000 MB"
		time : 1
		disks : "local-disk 50 HDD"		
	}
}
