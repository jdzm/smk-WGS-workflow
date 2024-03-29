
rule control_freec:
	input:
		bam = '%s/{s}/aligned/raw.mdups.recal.bam' % (in_data), 
		cont = lambda wc: '%s/%s/aligned/raw.mdups.recal.bam' % (in_data, samples[wc.s]['control'])
	output: 
		con = '%s/{s}/freec_AltCtrl/config_control.txt' % (derived),
		cnvs = '%s/{s}/freec_AltCtrl/freec_CNVs' % (derived),
		freec_ratio = '%s/{s}/freec_AltCtrl/freec_ratio.txt' % (derived),
		freec_info = '%s/{s}/freec_AltCtrl/freec_info.txt' % (derived),
	conda: '../envs/freec.yaml'
	threads: config['threads'] // 5
	log: '%s/{s}/05_freec_AltCtrl.log' % (logs)
	params: 
		outdir = '%s/{s}/freec_AltCtrl' % (derived), 
		chromsizes = genome['chrom_sizes']['main'],
		gchroms = genome['by_chroms'],
		gem_map = genome['gem_mapp'],
		template = '%s/freec_config/config_template_control.txt' % (metadata), 
		my_threads = config['threads'] // 5, 
		exp_ploidy = '2,3,4',
		g = genome['fasta'],
		snps = genome['dbsnp_vcf']
	priority: 10
	shell: 
		"""
		sed 's+MY_CHRSIZES+{params.chromsizes}+g' {params.template} \
		| sed 's+MY_CHR_FASTA+{params.gchroms}+g' \
		| sed 's+MY_GEM_MAP+{params.gem_map}+g' \
		| sed 's+MY_FULL_FASTA+{params.g}+g' \
		| sed 's+MY_SNPS+{params.snps}+g' \
		| sed 's+MY_THREADS+{params.my_threads}+g' \
		| sed 's+MY_INPUT+{input.bam}+g' \
		| sed 's+MY_OUTDIR+{params.outdir}+g' \
		| sed 's+MY_PLOIDY+{params.exp_ploidy}+g' \
		| sed 's+MY_CONTROL+{input.cont}+g' > {output.con}
        
		freec -conf {output.con} &> {log}
		rename raw.mdups.recal.bam_ freec_ {params.outdir}/*
        """

rule extra_files_freec:
	input:
		cnvs = '%s/{s}/freec_AltCtrl/freec_CNVs' % (derived),
		freec_ratio = '%s/{s}/freec_AltCtrl/freec_ratio.txt' % (derived)
	output: 
		ratio_plot = '%s/{s}/freec_AltCtrl/freec_ratio.png' % (derived),
		cnvsignif = '%s/{s}/freec_AltCtrl/freec_CNVs.p.value.txt' % (derived),
		report_plot = '%s/cnv_calls/RPE_TP53/control-freec/freec_plots/%s/AltCtrl_run/{s}_freec_ratio.png' % (derived, config["WGD_Exp"])
	params: 
		outdir = '%s/{s}/freec_AltCtrl/' % (derived),  
		exp_ploidy = '2,3,4',
		sam_id = '{s}',
	conda: '../envs/freec.yaml'
	log: '%s/{s}/05_freec_AltCtrl_plots.log' % (logs)
	priority: 10
	shell: 
		"""
		cat scripts/CNV_significance.R | R --slave --args {input.freec_ratio} {input.cnvs} {params.outdir} &> {log}
		
		Rscript --vanilla scripts/plotPloidyChr.R {params.sam_id} {input.freec_ratio} {params.exp_ploidy} {params.outdir} &> {log}
		
		cp {output.ratio_plot} {output.report_plot}
		"""

# rule gather_freec_plots: 
# 	input:
# 		freec_plot = '%s/{s}/freec_{m}/freec_ratio.png' % (derived),
# 	output: 
# 		new_loc = '/mnt/data3/Juan/freec_plots/{wgd}/{m}_run/{s}.png' % (derived),
# 	shell: 
# 		"cp {input.freec_plot} {output.new_loc}"

# merge ratios and CNV calls from all samples
# rule merge_freec_control: 
# 	input: 
# 		cnvs = ['%s/%s/freec_control/freec_CNVs.p.value.txt' % (derived, t) for t in tumors],
# 		ratios = ['%s/%s/freec_control/freec_ratio.txt' % (derived, t) for t in tumors],
# 		info = ['%s/%s/freec_control/freec_info.txt' % (derived, t) for t in tumors]
# 	output:
# 		cnv_sam = '%s/cnv_calls/freec_AltCtrl/control-freec/cnvs_control.tsv' % (derived),
# 		ratio_mat = '%s/cnv_calls/freec_AltCtrl/control-freec/ratio_control.tsv' % (derived), 
# 		freec_info = '%s/cnv_calls/freec_AltCtrl/control-freec/freec_info_control.tsv' % (derived)
# 	params:
# 		blacklist = genome["blacklist"],
# 		excl = genome["telocent"]
# 	priority: 10
# 	conda: '../envs/useR.yaml'
# 	script:
# 		'../scripts/mergeFreec.R'

# rule merge_freec_single: 
# 	input: 
# 		cnvs = ['%s/%s/freec_single/freec_CNVs.p.value.txt' % (derived, s) for s in samples],
# 		ratios = ['%s/%s/freec_single/freec_ratio.txt' % (derived, s) for s in samples],
# 		info = ['%s/%s/freec_single/freec_info.txt' % (derived, s) for s in samples]
# 	output:
# 		cnv_sam = '%s/cnv_calls/freec_AltCtrl/control-freec/cnvs_single.tsv' % (derived),
# 		ratio_mat = '%s/cnv_calls/freec_AltCtrl/control-freec/ratio_single.tsv' % (derived), 
# 		freec_info = '%s/cnv_calls/freec_AltCtrl/control-freec/freec_info_single.tsv' % (derived)
# 	params:
# 		blacklist = genome["blacklist"],
# 		excl = genome["telocent"]
# 	conda: '../envs/useR.yaml'
# 	priority: 10
# 	script:
# 		'../scripts/mergeFreec.R'	



