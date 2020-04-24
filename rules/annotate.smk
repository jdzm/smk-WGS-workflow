#### FUNCOTATOR
rule funcotate_vcf: 
	input: 
		'%s/{s}/mutect_calls/somatic.vcf.gz' % (derived)
	ouput:
		passfilt = '%s/{s}/mutect_calls/filtered_PASS.vcf' % (derived),
		vcf = '%s/{s}/annotate/somatic_func.vcf' % (derived),
		maf = '%s/{s}/annotate/somatic_func.maf' % (derived),
	conda: '../envs/bwa-gatk.yaml'
	params:
		g = config["ref_genome"]["bwa_idx"],
		func_path = config["annotate_vars"]["func_dir"]
	log: '%s/{s}/08_funcotate.log' % (logs)
	shell:
		"""
		bcftools view -f PASS {input} > {output.passfilt}
		
		gatk Funcotator 2>{log} --variant {input} --reference {params.g} --ref-version hg19 \
			--data-sources-path {params.func_path} --output {output.vcf} --output-file-format VCF \
			--remove-filtered-variants true
		
		gatk Funcotator 2>>{log} --variant {input} --reference {params.g} --ref-version hg19 \
			--data-sources-path {params.func_path} --output {output.maf} --output-file-format MAF \
			--remove-filtered-variants true
		"""

#### VEP
# rule vcf2maf:
# 	input: 
# 		'%s/{s}/mutect_calls/somatic.vcf.gz' % (derived)
# 	ouput:
# 		unzip_vcf = '%s/{s}/mutect_calls/somatic.vcf.gz' % (derived)
# 		passfilt = '%s/{s}/mutect_calls/filtered_PASS.vcf' % (derived),
# 		vcf = '%s/{s}/annotate/somatic_func.vcf' % (derived),
# 		maf = '%s/{s}/annotate/somatic_func.maf' % (derived),
# 	conda: '../envs/vep.yaml'
# 	params:
# 		g = config["ref_genome"]["bwa_idx"],
# 		func_path = config["annotate_vars"]["func_dir"],
# 		dbsnp = config["ref_genome"]["dbsnp_vcf"],
# 		gnomadaf = config["ref_genome"]["gnomad_af_only"],
# 		sam_id = '{s}', 
# 		con_id = samples[wc.s]['control']
# 	log: '%s/{s}/09_vep.log' % (logs)
# 	shell:
# 		"""
# 		### Hardcoded bit to run vcf2maf and vep
# 		export PERL_BASE="/home/diazmiya/miniconda2/envs/vep"
# 		export PERL5LIB="$PERL_BASE/lib/site_perl/5.26.2/:$PERL_BASE/lib/site_perl/5.26.2/x86_64-linux-thread-multi:$PERL5LIB" # export PERL5LIB="$PERL_BASE/perl-5.22.2/lib/perl5:$PERL_BASE/perl-5.22.2/lib:$PERL_BASE/perl-5.22.2/lib/perl5/x86_64-linux:$PERL5LIB"
# 		export PATH="$PERL_BASE/bin/:$PATH"
# 		export VEP_PATH="/home/diazmiya/miniconda2/envs/vep/share/ensembl-vep-99.2-0"
# 		export VEP_DATA="/mnt/data3/Juan/repository/ensembleVEP/hg19"
# 		export PERL5LIB=$VEP_PATH:$PERL5LIB
# 		export PATH=$VEP_PATH/htslib:$PATH
		
# 		vcf_file={{input}%.gz}
# 		echo $vcf_file
# 		# gunzip -c ${vcf_file}".gz" > ${vcf_file}
# 		# vcf2maf.pl --input-vcf ${vcf_file} --output-maf ${vcf_file%.vcf}.maf \
# 		# 	--tumor-id TP53_Ctrl --normal-id WT_Ctrl \
# 		#  	--vcf-tumor-id TP53_Ctrl --vcf-normal-id WT_Ctrl \
# 		#  	--ref-fasta {params.g} --vep-data $VEP_DATA --vep-path $VEP_PATH \
# 		#  	--filter-vcf {params.dbsnp} --maf-center ISREC 

		"""