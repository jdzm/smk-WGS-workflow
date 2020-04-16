rule igv_files:
	input: 
		bam = '%s/{s}/aligned/filt.bam' % (in_data),
	output:
		tdf = '%s/{s}.filt.bam.tdf' % (config["igv_vis"]),
	threads: 1
	conda: '../envs/trans.yaml'
	params:
		chromsizes = config["chr_sizes"]["main"]
	log:
		'%s/{s}/04_igvtools.log' % (logs)
	shell: 
		"""
		igvtools count 2> {log} \
			-f min,max,mean {input.bam} {output.tdf} {params.chromsizes}
		"""

rule vcf_to_bedpe:
	input: 
		vcf1 = '%s/{s}/delly/filtered_calls.vcf' % (derived),
	output:
		bedpe = '%s/{s}/delly/filtered_calls.bedpe' % (derived),
	conda: '../envs/trans.yaml'
	shell: 
		"""
		SURVIVOR vcftobed {input.vcf} 0 -1 {output.bedpe}
		"""

### SURVIVOR has many other options for vcf comparison
# -- Comparison/filtering
# 	merge	Compare or merge VCF files to generate a consensus or multi sample vcf files.
# 	filter	Filter a vcf file based on size and/or regions to ignore
# 	stats	Report multipe stats over a VCF file
# 	compMUMMer	Annotates a VCF file with the breakpoints found with MUMMer (Show-diff).