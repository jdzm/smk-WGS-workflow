#merge ratios and CNV calls from control freec
rule merge_freec_control: 
	input: 
		cnvs = ['%s/%s/freec_control/freec_CNVs.p.value.txt' % (derived, t) for t in tumors],
		ratios = ['%s/%s/freec_control/freec_ratio.txt' % (derived, t) for t in tumors],
		info = ['%s/%s/freec_control/freec_info.txt' % (derived, t) for t in tumors]
	output:
		cnv_sam = '%s/cnv_calls/{m}/control-freec/cnvs_control.tsv' % (derived),
		ratio_mat = '%s/cnv_calls/{m}/control-freec/ratio_control.tsv' % (derived), 
		freec_info = '%s/cnv_calls/{m}/control-freec/freec_info_control.tsv' % (derived)
	params:
		blacklist = config["ref_genome"]["blacklist"],
		excl = '%s/hg19.delly.excl.tsv' % (metadata)
	priority: 10
	conda: '../envs/useR.yaml'
	script:
		'../scripts/mergeFreec.R'

rule merge_freec_single: 
	input: 
		cnvs = ['%s/%s/freec_single/freec_CNVs.p.value.txt' % (derived, s) for s in samples],
		ratios = ['%s/%s/freec_single/freec_ratio.txt' % (derived, s) for s in samples],
		info = ['%s/%s/freec_single/freec_info.txt' % (derived, s) for s in samples]
	output:
		cnv_sam = '%s/cnv_calls/{m}/control-freec/cnvs_single.tsv' % (derived),
		ratio_mat = '%s/cnv_calls/{m}/control-freec/ratio_single.tsv' % (derived), 
		freec_info = '%s/cnv_calls/{m}/control-freec/freec_info_single.tsv' % (derived)
	params:
		blacklist = config["ref_genome"]["blacklist"],
		excl = '%s/hg19.delly.excl.tsv' % (metadata)
	conda: '../envs/useR.yaml'
	priority: 10
	script:
		'../scripts/mergeFreec.R'	