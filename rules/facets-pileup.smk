# Facets (MSKCC) calling of allele-specific copy number variants
# Use conda to download bioconda r-facets

#1. Count preparation. Generate SNP pileup either with the dbsnp file or
#	our own variant calls. Code from HTSlib
# snp-pileup -g <dbsnp.vcf.gz> <outputfile.csv.gz> <sequence files...bam>
# -g option produces gzip output
# -d maximum depth is set to 4000 
# -P, --pseudo-snps=MULTIPLE Every MULTIPLE positions, if there is no SNP,
                             # insert a blank record with the total count at the
                             # position.
# -r, --min-read-counts=READS   Comma separated list of minimum read counts for
#                              a position to be output. Default is 0.
# Run on two sumples just to test

rule snp_pileup:
	input: 
		['%s/%s/aligned/filt.bam' % (in_data, t) for t in tumors]
	output:
		'%s/snv_calls/RPE_TP53/pileup/snp_pileup_all.csv.gz' % (derived)
	threads: config["threads"]
	conda:
		'../envs/freec-useR.yaml'
	params:
		dbsnp = config["ref_genome"]["dbsnp_vcf"],
		mincounts = 5,
		con_bam = '%s/WT_Ctrl/aligned/filt.bam' % (in_data)
	log:
		'%s/joint/RPE_TP53-bypass/01_mutect_joint.log' % (logs)
	shell:
		"""
		# -P, --pseudo-snps=MULTIPLE Every M positions, if there is no SNP, insert a blank record with the total count at the position.
		snp-pileup -g -P100 -d2000 -Q20 -r{params.mincounts} \
			{params.vcf} {output} {params.con_bam} {input}
		"""

rule run_facets:
	input: 
		pileup = '%s/snv_calls/RPE_TP53/pileup/snp_pileup_all.csv.gz' % (derived)
	output:
		plot = '%s/summary_facets_ploidy.png' % (derived),
		diagplot = '%s/diag_spider.png' % (derived),
		'%s/snv_calls/RPE_TP53/pileup/snp_pileup_all.csv.gz' % (derived)
	threads: config["threads"]
	conda:
		'../envs/freec-useR.yaml'
	params:
		dbsnp = config["ref_genome"]["dbsnp_vcf"],
		mincounts = 5,
		con_bam = '%s/WT_Ctrl/aligned/filt.bam' % (in_data)
	log:
		'%s/joint/RPE_TP53-bypass/01_mutect_joint.log' % (logs)
	script:
		'../scripts/facets.R'







