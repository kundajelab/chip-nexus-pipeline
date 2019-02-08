# ENCODE DCC ChIP-Seq pipeline tester for task filter
# Author: Jin Lee (leepc12@gmail.com)
import "../../chip.wdl" as chip

workflow test_filter {
	String dup_marker = 'picard'	# picard.jar MarkDuplicates (picard) or 
									# sambamba markdup (sambamba)
	Int mapq_thresh = 30			# threshold for low MAPQ reads removal
	Boolean no_dup_removal = false	# no dupe reads removal when filtering BAM
									# dup.qc and pbc.qc will be empty files
									# and nodup_bam in the output is 
									# filtered bam with dupes	
	String pe_bam
	String se_bam

	String ref_pe_nodup_bam
	String ref_pe_filt_bam
	String ref_se_nodup_bam
	String ref_se_filt_bam
	String mito_chr_name = 'chrM'

	Int filter_cpu = 1
	Int filter_mem_mb = 20000
	Int filter_time_hr = 24
	String filter_disks = "local-disk 100 HDD"

	call chip.filter as pe_filter { input :
		bam = pe_bam,
		no_dup_removal = false,
		paired_end = true,
		dup_marker = dup_marker,
		mapq_thresh = mapq_thresh,
		mito_chr_name = mito_chr_name,

		cpu = filter_cpu,
		mem_mb = filter_mem_mb,
		time_hr = filter_time_hr,
		disks = filter_disks,
	}
	call chip.filter as pe_filter_no_dup_removal { input :
		bam = pe_bam,
		no_dup_removal = true,
		paired_end = true,
		dup_marker = dup_marker,
		mapq_thresh = mapq_thresh,
		mito_chr_name = mito_chr_name,

		cpu = filter_cpu,
		mem_mb = filter_mem_mb,
		time_hr = filter_time_hr,
		disks = filter_disks,
	}
	call chip.filter as se_filter { input :
		bam = se_bam,
		no_dup_removal = false,
		paired_end = false,
		dup_marker = dup_marker,
		mapq_thresh = mapq_thresh,
		mito_chr_name = mito_chr_name,

		cpu = filter_cpu,
		mem_mb = filter_mem_mb,
		time_hr = filter_time_hr,
		disks = filter_disks,
	}
	call chip.filter as se_filter_no_dup_removal { input :
		bam = se_bam,
		no_dup_removal = true,
		paired_end = false,
		dup_marker = dup_marker,
		mapq_thresh = mapq_thresh,
		mito_chr_name = mito_chr_name,

		cpu = filter_cpu,
		mem_mb = filter_mem_mb,
		time_hr = filter_time_hr,
		disks = filter_disks,
	}

	call chip.compare_md5sum { input :
		labels = [
			'pe_filter',
			'pe_filter_no_dup_removal',
			'se_filter',
			'se_filter_no_dup_removal',
		],
		files = [
			pe_filter.nodup_bam,
			pe_filter_no_dup_removal.nodup_bam,
			se_filter.nodup_bam,
			se_filter_no_dup_removal.nodup_bam,
		],
		ref_files = [
			ref_pe_nodup_bam,
			ref_pe_filt_bam,
			ref_se_nodup_bam,
			ref_se_filt_bam,
		],
	}
}
