rule delly_call:
	input: 
		bam = '%s/{s}/aligned/filt.bam' % (in_data), 
		control_bam = lambda wc: '%s/%s/aligned/filt.bam' % (in_data, samples[wc.s]['control'])
	output:
		bcf = '%s/{s}/delly/raw_calls.bcf' % (derived),
		svprops = '%s/{s}/delly/raw_calls.tab' % (derived),
		sampleprops = '%s/{s}/delly/raw_calls.sampleprops.tsv' % (derived)
	threads: config['threads'] // 4
	conda:
		'../envs/trans.yaml'
	params:
		g = config["ref_genome"]["bwa_idx"],
		excl = '%s/hg19.delly.excl.tsv' % (metadata), 
		sv_root = config["svprops"]
	log:
		'%s/{s}/06_delly_call.log' % (logs)
	shell:
		"""
		touch {input.bam}.bai
		delly call -x {params.excl} -o {output.bcf} -g {params.g} \
			{input.bam} {input.control_bam} &> {log}
		{params.sv_root}/svprops {output.bcf} > {output.svprops}
		{params.sv_root}/sampleprops {output.bcf} > {output.sampleprops}
		"""
		
rule delly_filter:
	input: 
		bcf_raw = '%s/{s}/delly/raw_calls.bcf' % (derived),
	output:
		bcf = '%s/{s}/delly/filtered_calls.bcf' % (derived),
		vcf = '%s/{s}/delly/filtered_calls.vcf' % (derived),
		svprops = '%s/{s}/delly/filtered_calls.tab' % (derived),
		sampleprops = '%s/{s}/delly/filtered_calls.sampleprops.tsv' % (derived)
	threads: config['threads'] // 6
	conda:
		'../envs/trans.yaml'
	params:
		sam_info = '%s/delly_samples_RPE.tsv' % (metadata),
		sv_root = config["svprops"]
	log:
		'%s/{s}/06_delly_filter.log' % (logs)
	shell:
		"""
		# no support in the matched normal and an overall confident VCF filter equal to PASS.
		# Not in place. just keeping code > -a 0.25 # equire a minimum variant allele frequency of 25%
		
		delly filter -p -f somatic -o {output.bcf} \
			-s {params.sam_info} {input.bcf_raw} &> {log}
		{params.sv_root}/svprops {output.bcf} > {output.svprops}
		{params.sv_root}/sampleprops {output.bcf} > {output.sampleprops}

		bcftools convert -O v -o {output.vcf} {output.bcf} &>> {log}
		"""

# one might need to reindex all bam files before running delly 
rule delly_somatic_joint:
	input: 
		bam = ['%s/%s/aligned/filt.bam' % (in_data, t) for t in tumors],
		control_bam = '%s/%s/aligned/filt.bam' % (in_data, "WT_Ctrl")
	output:
		bcf = '%s/sv_calls/RPE_TP53/delly/raw_calls.bcf' % (derived),
		svprops = '%s/sv_calls/RPE_TP53/delly/raw_calls.tab' % (derived),
		sampleprops = '%s/sv_calls/RPE_TP53/delly/raw_calls.sampleprops.tsv' % (derived)
	threads: config['threads']
	conda:
		'../envs/trans.yaml'
	params:
		g = config["ref_genome"]["bwa_idx"],
		excl = '%s/hg19.delly.excl.tsv' % (metadata), 
		sv_root = config["svprops"]
	log:
		'%s/joint/RPE_TP53/02_delly_joint.log' % (logs)
	shell:
		"""
		delly call -x {params.excl} -o {output.bcf} -g {params.g} \
			{input.bam} {input.control_bam} &> {log}
		{params.sv_root}/svprops {output.bcf} > {output.svprops}
		{params.sv_root}/sampleprops {output.bcf} > {output.sampleprops}
		"""

# Watch out for the sam_info in params. It is hardcoded
rule delly_filter_joint:
	input: 
		bcf_raw = '%s/sv_calls/RPE_TP53/delly/raw_calls.bcf' % (derived)
	output:
		bcf = '%s/sv_calls/RPE_TP53/delly/somatic_svs.bcf' % (derived),
		svprops = '%s/sv_calls/RPE_TP53/delly/somatic_svs.tab' % (derived),
		sampleprops = '%s/sv_calls/RPE_TP53/delly/somatic_svs.sampleprops.tsv' % (derived),
		vcf = '%s/sv_calls/RPE_TP53/delly/somatic_svs.vcf' % (derived),
	threads: 1
	conda:
		'../envs/trans.yaml'
	params:
		sam_info = '%s/delly_samples_RPE_oldnames.tsv' % (metadata),
		sv_root = config["svprops"]
	log:
		'%s/joint/RPE_TP53/03_delly_filter.log' % (logs)
	shell:
		"""
		touch {input.bcf_raw}.csi
		# no support in the matched normal and an overall confident VCF filter equal to PASS.
		# Not in place. just keeping code > -a 0.25 # equire a minimum variant allele frequency of 25%
		
		delly filter -p -f somatic -o {output.bcf} \
			-s {params.sam_info} {input.bcf_raw} &> {log}
		{params.sv_root}/svprops {output.bcf} > {output.svprops}
		{params.sv_root}/sampleprops {output.bcf} > {output.sampleprops}

		bcftools convert -O v -o {output.vcf} {output.bcf} &>> {log}
		"""

# rule to_vcf: 
# 	input: '%s/bcf_file.bcf' % (derived)
# 	output: '%s/vcf_file.vcf' % (derived)
# 	conda: '../envs/trans.yaml'
# 	threads: 1
# 	shell:
# 		"""
# 		# outputs uncompressed vcf file
# 		bcftools convert -O v -o {output} {input}
# 		"""



