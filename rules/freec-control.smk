#lambda wc: samples[wc.s]['expectedPloidy']
rule control_freec:
	input:
		bam = '%s/{s}/aligned/filt.bam' % (in_data), 
		cont = lambda wc: '%s/%s/aligned/filt.bam' % (in_data, samples[wc.s]['control'])
	output: 
		con = '%s/{s}/freec_{m}/config_{m}.txt' % (derived),
		cnvs = '%s/{s}/freec_{m}/freec_CNVs' % (derived),
		freec_ratio = '%s/{s}/freec_{m}/freec_ratio.txt' % (derived)
	conda: '../envs/freec.yaml'
	threads: config['threads'] // 5
	log: '%s/{s}/05_freec_{m}.log' % (logs)
	params: 
		outdir = '%s/{s}/freec_{m}' % (derived), 
		gname = genome['name'],
		template = '%s/freec_config/config_template_{m}.txt' % (metadata), 
		my_threads = config['threads'] // 5, 
		exp_ploidy = '2,3,4'
	priority: 10
	shell: 
		"""
		sed 's/REF_GENOME/{params.gname}/g' {params.template} \
		| sed 's+MY_THREADS+{params.my_threads}+g' \
		| sed 's+MY_INPUT+{input.bam}+g' \
		| sed 's+MY_OUTDIR+{params.outdir}+g' \
		| sed 's+MY_PLOIDY+{params.exp_ploidy}+g' \
		| sed 's+MY_CONTROL+{input.cont}+g' > {output.con}
        
		freec -conf {output.con} &> {log}
		rename filt.bam_ freec_ {params.outdir}/*
        """

rule extra_files_freec:
	input:
		cnvs = '%s/{s}/freec_{m}/freec_CNVs' % (derived),
		freec_ratio = '%s/{s}/freec_{m}/freec_ratio.txt' % (derived)
	output: 
		ratio_plot = '%s/{s}/freec_{m}/freec_ratio.png' % (derived),
		cnvsignif = '%s/{s}/freec_{m}/freec_CNVs.p.value.txt' % (derived)
	params: 
		outdir = '%s/{s}/freec_{m}/' % (derived),  
		exp_ploidy = '2,3,4',
		sam_id = '{s}'
	conda: '../envs/freec.yaml'
	log: '%s/{s}/05_freec_{m}_plots.log' % (logs)
	priority: 10
	shell: 
		"""
		cat scripts/CNV_significance.R | R --slave --args {input.freec_ratio} {input.cnvs} {params.outdir} &> {log}
		
		Rscript --vanilla scripts/plotPloidyChr.R {params.sam_id} {input.freec_ratio} {params.exp_ploidy} {params.outdir} &> {log}
		"""






