rule mutect2_call:
	input: 
		bam = '%s/{s}/aligned/filt.bam' % (in_data), 
		cont = lambda wc: '%s/%s/aligned/filt.bam' % (in_data, samples[wc.s]['control'])
	output:
		vcf = '%s/{s}/mutect_calls/somatic_raw.vcf.gz' % (derived),
		orient = '%s/{s}/mutect_calls/f1r2.tar.gz' % (derived),
		romodel = '%s/{s}/mutect_calls/read-orientation-model.tar.gz' % (derived)
	threads: config['threads'] // 3
	conda:
		'../envs/bwa-gatk.yaml'
	params:
		g = config["ref_genome"]["bwa_idx"],
		dbsnp = config["ref_genome"]["dbsnp_vcf"],
		gnomadaf = config["ref_genome"]["gnomad_af_only"],
		con_name = lambda wc: '%s' % (samples[wc.s]['control'])
	log:
		'%s/{s}/07_mutect.log' % (logs)
	shell:
		"""
		gatk Mutect2 2> {log} \
			-R {params.g} \
			-I {input.bam} \
			-I {input.cont} \
			--native-pair-hmm-threads {threads} \
			--f1r2-tar-gz {output.orient} \
			--germline-resource {params.gnomadaf} \
			-normal {params.con_name} \
			-O {output.vcf}

		gatk LearnReadOrientationModel 2> {log} \
			-I {output.orient} -O {output.romodel}
		"""
		# --panel-of-normals pon.vcf.gz \
		# The -f1r2 file is needed to improve the filtering based on the read orientation model. 
		

rule mutect2_filter:
	input: 
		vcf = '%s/{s}/mutect_calls/somatic_raw.vcf.gz' % (derived),
		romodel = '%s/{s}/mutect_calls/read-orientation-model.tar.gz' % (derived)
	output:
		filt = '%s/{s}/mutect_calls/somatic.vcf.gz' % (derived)
	threads: config["threads"]
	conda:
		'../envs/bwa-gatk.yaml'
	params:
		g = config["ref_genome"]["bwa_idx"],
		dbsnp = config["ref_genome"]["dbsnp_vcf"]
	log:
		'%s/{s}/075_mutect_filter.log' % (logs)
	shell:
		"""
		
		gatk FilterMutectCalls 2> {log} \
			-R {params.g} \
			-V {input.vcf} \
			--contamination-estimate 0.0 \
			--ob-priors {input.romodel} \
			-O {output.filt}

		"""
		# add contamination part later. 


#expand (['%s/{t}/aligned/filt.bam' % (in_data)], t=tumors)
rule mutect2_joint_call:
	input: 
		['%s/%s/aligned/filt.bam' % (in_data, t) for t in tumors]
	output:
		vcf = '%s/snv_calls/{m}/mutect2/somatic_raw.vcf.gz' % (derived)
	threads: config["threads"]
	conda:
		'../envs/bwa-gatk.yaml'
	params:
		g = config["ref_genome"]["bwa_idx"],
		dbsnp = config["ref_genome"]["dbsnp_vcf"],
		con_name = "WT_Ctrl",
		con_bam = '%s/WT_Ctrl/aligned/filt.bam' % (in_data)
	log:
		'%s/joint/{m}/01_mutect_joint.log' % (logs)
	shell:
		"""
		in_files="{input}"
		deconvolve_inputs=${{in_files// /" -I "}}
		
		gatk Mutect2 2> {log} \
			-R {params.g} \
			-I $deconvolve_inputs \
			-I {params.con_bam} \
			--native-pair-hmm-threads {threads} \
			-normal {params.con_name} \
			--germline-resource {params.dbsnp} \
			-O {output.vcf}
		"""

rule mutect2_joint_filter:
	input: 
		vcf = '%s/snv_calls/{m}/mutect2/somatic_raw.vcf.gz' % (derived)
	output:
		filt = '%s/snv_calls/{m}/mutect2/somatic.vcf.gz' % (derived)
	threads: config["threads"]
	conda:
		'../envs/bwa-gatk.yaml'
	params:
		g = config["ref_genome"]["bwa_idx"]
	log:
		'%s/joint/{m}/015_mutect_joint_filter.log' % (logs)
	shell:
		"""
		
		gatk FilterMutectCalls 2> {log} \
			-R {params.g} \
			-V {input.vcf} \
			--contamination-estimate 0.0 \
			-O {output.filt} 

		"""
# 		--unique-alt-read-count 4 \


# rule mutect2_refine:
# 	input: 
# 		bam = '%s/{s}/aligned/filt.bam' % (in_data), 
# 		orient = '%s/{s}/mutect/f1r2.tar.gz' % (derived)
# 	output:
# 		romodel = '%s/{s}/mutect/read-orientation-model.tar.gz' % (derived),
# 		pilesum = '%s/{s}/mutect/pileupsummaries.table' % (derived),
# 		contam = '%s/{s}/mutect/contamination.table' % (derived),
		
# 	threads: config['threads'] // 4
# 	conda:
# 		'../envs/bwa-gatk.yaml'
# 	params:
# 		g = config["ref_genome"]["bwa_idx"],
# 		dbsnp = config["ref_genome"]["dbsnp_vcf"]
# 	log:
# 		'%s/{s}/075_mutect_refine.log' % (logs)
# 	shell:
# 		"""
# 		gatk LearnReadOrientationModel 2> {log} \
# 			-I {input.orient} -O {output.romodel}
		
# 		gatk GetPileupSummaries 2> {log} \
# 		    -I {input.bam} \
# 		    -V {params.g} \ # might need to be the biallelic version!!
# 		    -L {params.g} \
# 		    -O {output.pilesum}
		
# 		gatk CalculateContamination 2> {log} \
# 	        -I getpileupsummaries.table \
# 	        -tumor-segmentation segments.table \
# 	        -O calculatecontamination.table

# 		"""










