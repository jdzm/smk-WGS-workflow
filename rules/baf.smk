rule baf_freec:
	input:
		bam = '%s/{s}/aligned/raw.mdups.recal.bam' % (in_data), 
		cont = lambda wc: '%s/%s/aligned/raw.mdups.recal.bam' % (in_data, samples[wc.s]['control'])
	output: 
		con = '%s/{s}/freec_{m}/config_{m}.txt' % (derived),
		cnvs = '%s/{s}/freec_{m}/freec_CNVs' % (derived),
		freec_ratio = '%s/{s}/freec_{m}/freec_ratio.txt' % (derived),
		freec_info = '%s/{s}/freec_{m}/freec_info.txt' % (derived),
	conda: '../envs/freec.yaml'
	threads: config['threads'] // 5
	log: '%s/{s}/05_freec_{m}.log' % (logs)
	params: 
		outdir = '%s/{s}/freec_{m}' % (derived), 
		chromsizes = genome['chrom_sizes']['main'],
		gchroms = genome['by_chroms'],
		gem_map = genome['gem_mapp'],
		template = '%s/freec_config/config_template_{m}.txt' % (metadata), 
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