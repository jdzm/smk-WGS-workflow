rule freebayes_call:
	input: 
		'%s/{s}/aligned/filt.bam' % (in_data)
	output:
		'%s/{s}/freebayes/variants.vcf' % (derived)
	threads: config["threads"] // 4
	conda:
		'../envs/bwa-gatk.yaml'
	params:
		g = config["ref_genome"]["bwa_idx"]
	log:
		'%s/joint/RPE_TP53/03_freebayes_joint.log' % (logs)
	shell:
		"""
		freebayes -f {params.g} -g 200 -C 4 {input} > {output}
		"""
		
rule freebayes_joint_call:
	input: 
		['%s/%s/aligned/filt.bam' % (in_data, t) for t in tumors]
	output:
		vcf = '%s/snv_calls/RPE_TP53/freebayes/variants.vcf' % (derived)
	threads: config["threads"]
	conda:
		'../envs/bwa-gatk.yaml'
	params:
		g = config["ref_genome"]["bwa_idx"]
	log:
		'%s/joint/RPE_TP53/03_freebayes_joint.log' % (logs)
	shell:
		"""
		freebayes-parallel <(fasta_generate_regions.py {params.g}.fai 100000) {threads} \
			-f {params.g} -g 200 -C 4 {input} > {output.vcf}
		"""

# discards alignments on positions with cov>200 reads to reduce comp time
# Requires at least 4 supporting observations to consider a variant
# -p flag can be used to force different ploidy

# use an input VCF (bgzipped + tabix indexed) to force calls at particular alleles
# freebayes -f ref.fa -@ in.vcf.gz aln.bam >var.vcf

# Run freebayes in parallel on 100000bp chunks of the ref (fasta_generate_regions.py is also
# located in the scripts/ directory in the freebayes distribution).  Use 36 threads.
# freebayes-parallel <(fasta_generate_regions.py ref.fa.fai 100000) 36 -f ref.fa aln.bam >out.vcf

