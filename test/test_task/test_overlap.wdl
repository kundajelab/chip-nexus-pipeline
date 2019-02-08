# ENCODE DCC ChIP-Seq pipeline tester
# Author: Jin Lee (leepc12@gmail.com)
import "../../chip.wdl" as chip

workflow test_overlap {
	String se_peak_rep1 # test overlap,idr for SE set only
	String se_peak_rep2
	String se_peak_pooled
	String se_ta_pooled

	String ref_se_overlap_peak
	String ref_se_overlap_bfilt_peak
	String ref_se_overlap_frip_qc

	String se_blacklist
	String se_chrsz

	Int fraglen

	call chip.overlap as se_overlap { input :
		prefix = "rep1-rep2",
		peak1 = se_peak_rep1,
		peak2 = se_peak_rep2,
		peak_pooled = se_peak_pooled,
		peak_type = 'regionPeak',
		blacklist = se_blacklist,
		keep_irregular_chr_in_bfilt_peak = false,
		chrsz = se_chrsz,
		fraglen = fraglen,
		ta = se_ta_pooled,
	}

	call chip.compare_md5sum { input :
		labels = [
			'se_overlap_peak',
			'se_overlap_bfilt_peak',
			'se_overlap_frip_qc',
		],
		files = [
			se_overlap.overlap_peak,
			se_overlap.bfilt_overlap_peak,
			se_overlap.frip_qc,
		],
		ref_files = [
			ref_se_overlap_peak,
			ref_se_overlap_bfilt_peak,
			ref_se_overlap_frip_qc,
		],
	}
}
