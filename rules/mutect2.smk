## Runs mutect on all control-tumor pairs. Controls are defined in the config file.
## Currently using genomAD allele frequencies as germline resource, and excludion low comp
## regions (telocentromeric+blacklisted). f1r2 file is needed to improve filtering at later
## step. Running without panel of normals. 
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
		intervals = genome['mutect_exclude'],
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
			-XL {params.intervals} \
			--native-pair-hmm-threads {threads} \
			--f1r2-tar-gz {output.orient} \
			--germline-resource {params.gnomadaf} \
			-normal {params.con_name} \
			-O {output.vcf}

		gatk LearnReadOrientationModel 2>> {log} \
			-I {output.orient} -O {output.romodel}
		"""

# Adds filter fields to mutect result. Does not remove variants.
rule mutect2_filter:
	input: 
		vcf = '%s/{s}/mutect_calls/somatic_raw.vcf.gz' % (derived),
		romodel = '%s/{s}/mutect_calls/read-orientation-model.tar.gz' % (derived)
	output:
		filt = '%s/{s}/mutect_calls/somatic.vcf.gz' % (derived)
	threads: config["threads"] // 3
	conda:
		'../envs/bwa-gatk.yaml'
	params:
		g = genome["fasta"],
	log:
		'%s/{s}/07_mutect_filter.log' % (logs)
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

# Joint call of all the samples labelled as 'tumor'. Control is hardcoded (only one in our)
# experimental setup. It will require minor tweaking to include other controls. 
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
		intervals = genome['mutect_exclude'],
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
			-XL {params.intervals} \
			--native-pair-hmm-threads {threads} \
			--f1r2-tar-gz {output.orient} \
			--germline-resource {params.gnomadaf} \
			-normal {params.con_name} \
			-O {output.vcf}
		"""

rule mutect2_joint_filter:
	input: 
		vcf = '%s/snv_calls/{m}/mutect2-joint/somatic_raw.vcf.gz' % (derived),
		orient = '%s/snv_calls/{m}/mutect2/f1r2.tar.gz' % (derived)
	output:
		filt = '%s/snv_calls/{m}/mutect2-joint/somatic.vcf.gz' % (derived),
		romodel = '%s/snv_calls/{m}/mutect2-joint/read-orientation-model.tar.gz' % (derived)
	threads: config["threads"]
	conda:
		'../envs/bwa-gatk.yaml'
	params:
		g = genome["fasta"]
	log:
		'%s/joint/{m}/01_mutect_joint_filt.log' % (logs)
	shell:
		"""
		gatk LearnReadOrientationModel 2> {log} \
			-I {input.orient} -O {output.romodel}

		gatk FilterMutectCalls 2>> {log} \
			-R {params.g} \
			-V {input.vcf} \
			--contamination-estimate 0.0 \
			--ob-priors {output.romodel} \
			-O {output.filt}

		"""



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










