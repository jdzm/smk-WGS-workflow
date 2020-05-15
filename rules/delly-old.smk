rule delly_call:
	input: 
		bam = '%s/{s}/aligned/raw.mdups.recal.bam' % (in_data), 
		control_bam = lambda wc: '%s/%s/aligned/raw.mdups.recal.bam' % (in_data, samples[wc.s]['control'])
	output:
		bcf = '%s/{s}/delly/raw_calls.bcf' % (derived),
		csi = '%s/{s}/delly/raw_calls.bcf.csi' % (derived),
		vcf = '%s/{s}/delly/raw_calls.vcf.gz' % (derived),
		tbi = '%s/{s}/delly/raw_calls.vcf.gz.tbi' % (derived),
		svprops = '%s/{s}/delly/raw_calls.tab' % (derived),
		sampleprops = '%s/{s}/delly/raw_calls.sampleprops.tsv' % (derived)
	threads: 2
	conda:
		'../envs/trans.yaml'
	params:
		g = genome["fasta"],
		excl = '%s/delly.excl.tsv' % (metadata), 
		intervals = genome['mutect_exclude'],
		sv_root = config["svprops"]
	log:
		'%s/{s}/06_delly_call.log' % (logs)
	shell:
		"""
		export OMP_NUM_THREADS=2
		delly_precomp call -x {params.intervals} -o {output.bcf} -g {params.g} \
			{input.bam} {input.control_bam} &> {log}
		
		bcftools convert -O z -o {output.vcf} {output.bcf} &>> {log}
		tabix {output.vcf}

		{params.sv_root}/svprops {output.bcf} > {output.svprops}
		{params.sv_root}/sampleprops {output.bcf} > {output.sampleprops}
		"""
		
rule delly_filter:
	input: 
		bcf_raw = '%s/{s}/delly/raw_calls.bcf' % (derived),
		csi_raw = '%s/{s}/delly/raw_calls.bcf.csi' % (derived)
	output:
		bcf = temp('%s/{s}/delly/somatic_svs.bcf' % (derived)),
		csi = temp('%s/{s}/delly/somatic_svs.bcf.csi' % (derived)),
		vcf = '%s/{s}/delly/somatic_svs.vcf.gz' % (derived),
		tbi = '%s/{s}/delly/somatic_svs.vcf.gz.tbi' % (derived),
		svprops = '%s/{s}/delly/somatic_svs.tab' % (derived),
		sampleprops = '%s/{s}/delly/somatic_svs.sampleprops.tsv' % (derived),
	threads: 2
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
		# Not in place. just keeping code > -a 0.01 # require a minimum variant allele frequency of 25%
		
		delly filter -a 0.01 -p -f somatic -o {output.bcf} \
			-s {params.sam_info} {input.bcf_raw} &> {log}
		{params.sv_root}/svprops {output.bcf} > {output.svprops}
		{params.sv_root}/sampleprops {output.bcf} > {output.sampleprops}

		bcftools convert -O z -o {output.vcf} {output.bcf} &>> {log}
		tabix {output.vcf}
		"""

# one might need to reindex all bam files before running delly 
rule delly_somatic_joint:
	input: 
		bam = ['%s/%s/aligned/raw.mdups.recal.bam' % (in_data, t) for t in tumors],
		control_bam = '%s/%s/aligned/raw.mdups.recal.bam' % (in_data, "WT_Ctrl")
	output:
		bcf = '%s/sv_calls/RPE_TP53/delly/raw_calls.bcf' % (derived),
		svprops = '%s/sv_calls/RPE_TP53/delly/raw_calls.tab' % (derived),
		sampleprops = '%s/sv_calls/RPE_TP53/delly/raw_calls.sampleprops.tsv' % (derived)
	threads: config['threads']
	conda:
		'../envs/trans.yaml'
	params:
		g = genome["fasta"],
		excl = '%s/delly.excl.tsv' % (metadata), 
		sv_root = config["svprops"]
	log:
		'%s/joint/RPE_TP53/02_delly_joint.log' % (logs)
	shell:
		"""
		export OMP_NUM_THREADS=7
		{params.run_delly} call -x {params.excl} -o {output.bcf} -g {params.g} \
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
		sam_info = '%s/delly_samples_RPE.tsv' % (metadata),
		sv_root = config["svprops"]
	log:
		'%s/joint/RPE_TP53/03_delly_filter.log' % (logs)
	shell:
		"""
		touch {input.bcf_raw}.csi
		# no support in the matched normal and an overall confident VCF filter equal to PASS.
		# Not in place. just keeping code > -a 0.01 # equire a minimum variant allele frequency of 1%
		
		delly filter -a 0.01 -p -f somatic -o {output.bcf} \
			-s {params.sam_info} {input.bcf_raw} &> {log}
		{params.sv_root}/svprops {output.bcf} > {output.svprops}
		{params.sv_root}/sampleprops {output.bcf} > {output.sampleprops}

		bcftools convert -O v -o {output.vcf} {output.bcf} &>> {log}
		"""


# 		"""
# 		# outputs uncompressed vcf file
# 		bcftools convert -O v -o {output} {input}
# 		"""


