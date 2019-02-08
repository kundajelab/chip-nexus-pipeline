# ENCODE DCC ChIP-Seq pipeline tester for task spp
# Author: Jin Lee (leepc12@gmail.com)
import "../../chip.wdl" as chip

workflow test_spp {
	Int cap_num_peak

	Int fraglen
	# test spp for SE set only
	String se_ta
	String se_ctl_ta

	String ref_se_spp_rpeak # raw narrow-peak
	String ref_se_spp_bfilt_rpeak # blacklist filtered narrow-peak
	String ref_se_spp_frip_qc 

	String se_blacklist
	String se_chrsz

	Int spp_cpu = 1
	Int spp_mem_mb = 16000
	Int spp_time_hr = 72
	String spp_disks = "local-disk 100 HDD"

	call chip.spp as se_spp { input :
		tas = [se_ta, se_ctl_ta],
		chrsz = se_chrsz,
		fraglen = fraglen,
		cap_num_peak = cap_num_peak,
		blacklist = se_blacklist,
		keep_irregular_chr_in_bfilt_peak = false,

		cpu = spp_cpu,
		mem_mb = spp_mem_mb,
		time_hr = spp_time_hr,
		disks = spp_disks,
	}

	call chip.compare_md5sum { input :
		labels = [
			'se_spp_rpeak',
			'se_spp_bfilt_rpeak',
			'se_spp_frip_qc',
		],
		files = [
			se_spp.rpeak,
			se_spp.bfilt_rpeak,
			se_spp.frip_qc,
		],
		ref_files = [
			ref_se_spp_rpeak,
			ref_se_spp_bfilt_rpeak,
			ref_se_spp_frip_qc,
		],
	}
}
