rule mutect2_call:
	input: 
		bam = '%s/{s}/aligned/raw.mdups.recal.bam' % (in_data), 
		cont = lambda wc: '%s/%s/aligned/raw.mdups.recal.bam' % (in_data, samples[wc.s]['control'])
	output:
		vcf = '%s/{s}/mutect_calls/somatic_raw.vcf.gz' % (derived),
		orient = '%s/{s}/mutect_calls/f1r2.tar.gz' % (derived),
		romodel = '%s/{s}/mutect_calls/read-orientation-model.tar.gz' % (derived)
	threads: config['threads'] // 3
	conda:
		'../envs/bwa-gatk.yaml'
	params:
		g = genome["fasta"],
		gnomadaf = genome["gnomad_af_only"],
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
		g = genome["fasta"],
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

rule mutect2_joint_call:
	input: 
		bam = ['%s/%s/aligned/raw.mdups.recal.bam' % (in_data, t) for t in tumors]
	output:
		vcf = '%s/snv_calls/{m}/mutect2-gnomAD-joint/somatic_raw.vcf.gz' % (derived),
		orient = '%s/snv_calls/{m}/mutect2-gnomAD-joint/f1r2.tar.gz' % (derived)
	threads: config["threads"]
	conda:
		'../envs/bwa-gatk.yaml'
	params:
		g = genome["fasta"],
		gnomadaf = genome["gnomad_af_only"],
		con_name = "WT_Ctrl",
		con_bam = '%s/WT_Ctrl/aligned/raw.mdups.recal.bam' % (in_data)
	log:
		'%s/joint/{m}/01_mutect_joint_gnomad.log' % (logs)
	shell:
		"""
		in_files="{input.bam}"
		deconvolve_inputs=${{in_files// /" -I "}}

		gatk Mutect2 2> {log} \
			-R {params.g} \
			-I $deconvolve_inputs \
			-I {params.con_bam} \
			--native-pair-hmm-threads {threads} \
			--f1r2-tar-gz {output.orient} \
			--germline-resource {params.gnomadaf} \
			-normal {params.con_name} \
			-O {output.vcf}
		"""

rule mutect2_joint_filter:
	input: 
		vcf = '%s/snv_calls/{m}/mutect2-gnomAD-joint/somatic_raw.vcf.gz' % (derived),
		orient = '%s/snv_calls/{m}/mutect2-gnomAD-joint/f1r2.tar.gz' % (derived)
	output:
		filt = '%s/snv_calls/{m}/mutect2-gnomAD-joint/somatic.vcf.gz' % (derived),
		romodel = '%s/snv_calls/{m}/mutect2-gnomAD-joint/read-orientation-model.tar.gz' % (derived)
	threads: config["threads"]
	conda:
		'../envs/bwa-gatk.yaml'
	params:
		g = genome["fasta"]
	log:
		'%s/joint/{m}/015_mutect_joint_filter_gnomad.log' % (logs)
	shell:
		"""
		gatk LearnReadOrientationModel 2> {log} \
			-I {input.orient} -O {output.romodel}

		gatk FilterMutectCalls 2> {log} \
			-R {params.g} \
			-V {input.vcf} \
			--contamination-estimate 0.0 \
			--ob-priors {output.romodel} \
			-O {output.filt}

		"""
# 		--unique-alt-read-count 4 \


# rule mutect2_refine:
# When I have other sort of samples I will need to 
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










