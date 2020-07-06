# CHRS = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X']

## Runs mutect on all control-tumor pairs. Controls are defined in the config file.
## Currently using genomAD allele frequencies as germline resource, and excludion low comp
## regions (telocentromeric+blacklisted). f1r2 file is needed to improve filtering at later
## step. Running without panel of normals. 
rule def_output:
	input: 
		expand ('%s/{s}/mutect_calls/chrs/chr{chr}_somatic_raw.vcf.gz' % (derived), s=tumors, chr=CHRS)

rule mutect2_call:
	input: 
		bam = '%s/{s}/aligned/raw.mdups.recal.bam' % (in_data), 
		cont = lambda wc: '%s/%s/aligned/raw.mdups.recal.bam' % (in_data, samples[wc.s]['control'])
	output:
		vcf = '%s/{s}/mutect_calls/chrs/chr{chr}_somatic_raw.vcf.gz' % (derived),
		orient = '%s/{s}/mutect_calls/chrs/chr{chr}_f1r2.tar.gz' % (derived),
		romodel = '%s/{s}/mutect_calls/chrs/chr{chr}_read-orientation-model.tar.gz' % (derived)
	threads: 1
	conda:
		'../envs/bwa-gatk.yaml'
	params:
		g = genome["fasta"],
		intervals = genome['mutect_exclude'],
		gnomadaf = genome["gnomad_af_only"],
		con_name = lambda wc: '%s' % (samples[wc.s]['control']),
		chosen_chr = '{chr}'
	log:
		'%s/{s}/07_mutect_chr{chr}.log' % (logs)
	shell:
		"""
		gatk Mutect2 2> {log} \
			-R {params.g} \
			-I {input.bam} \
			-I {input.cont} \
			-L {params.chosen_chr} \
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
		vcf = '%s/{s}/mutect_calls/chrs/chr{chr}_somatic_raw.vcf.gz' % (derived),
		romodel = '%s/{s}/mutect_calls/chrs/chr{chr}_read-orientation-model.tar.gz' % (derived)
	output:
		filt = '%s/{s}/mutect_calls/chrs/chr{chr}_somatic.vcf.gz' % (derived)
	threads: 1
	conda:
		'../envs/bwa-gatk.yaml'
	params:
		g = genome["fasta"],
	log:
		'%s/{s}/07_mutect_filter_chr{chr}.log' % (logs)
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

rule merge_mutect2:
	input:
		expand ('%s/{{s}}/mutect_calls/chrs/chr{chr}_somatic.vcf.gz' % (derived), chr=CHRS),
	output:
		'%s/{s}/mutect_calls/somatic.vcf.gz' % (derived), 
	conda:
		'../envs/bwa-gatk.yaml'
	log:
		'%s/{s}/07_mutect_merge.log' % (logs)
	shell:
		"""
		in_files="{input}"
		deconvolve_inputs=${{in_files// /" -I "}}

		gatk MergeVcfs 2> {log} -I $deconvolve_inputs -O {output}
		"""

