#### FUNCOTATOR
rule funcotate_vcf: 
	input: 
		'%s/{s}/mutect_calls/somatic.vcf.gz' % (derived)
	output:
		vcf = '%s/{s}/annotate/somatic_func.vcf' % (derived),
		maf = '%s/{s}/annotate/somatic_func.maf' % (derived),
	conda: '../envs/bwa-gatk.yaml'
	params:
		g = genome["fasta"],
		func_path = config["annotate_vars"]["func_dir"]
	log: '%s/{s}/09_funcotate.log' % (logs)
	shell:
		"""		
		gatk Funcotator 2>{log} --variant {input} --reference {params.g} --ref-version hg19 \
			--data-sources-path {params.func_path} --output {output.vcf} --output-file-format VCF \
			--remove-filtered-variants true
		
		gatk Funcotator 2>>{log} --variant {input} --reference {params.g} --ref-version hg19 \
			--data-sources-path {params.func_path} --output {output.maf} --output-file-format MAF \
			--remove-filtered-variants true
		"""

#### VEP & vcf2maf
# This one is tricky... You need to download and install vep on the computer first. 
## conda env create --file ../envs/vep.yaml 
# Then execute the following:
## vep_install -a cf -s homo_sapiens -y GRCh37 -c $VEP_DATA -d $VEP_PATH --NO_HTSLIB -v $VER
## vep_convert_cache --species homo_sapiens --version 99_GRCh37 --dir $VEP_DATA
# This will install vep and download the reference libraries (~15G)
# Now, create the same conda env with the yaml file inside the smk rule 
# In the shell part, I export the perl paths and the vep paths so that vep is able to find where it is
# but I use the paths from the local vep environmen, the smk one will be generated separately and
# it has random names. The paths are HARDCODED below
## If it crashes and you need ot re-run, remember to delete the vep vcf file created in destination
rule vcf2maf:
	input: 
		'%s/{s}/mutect_calls/somatic.vcf.gz' % (derived)
	output:
		unzip_vcf = temp('%s/{s}/annotate/somatic.vcf' % (derived)),
		maf = '%s/{s}/annotate/somatic.maf' % (derived),
	threads: config["threads"] // 4
	params:
		g = genome['fasta'],
		gnomadaf = genome["gnomad_af_only"],
		sam_id = '{s}', 
		con_id = lambda wc: '%s' % (samples[wc.s]['control']),
		center ='ChezDiaz',
		vcf2maf = '/home/diazmiya/miniconda2/envs/vep/bin/vcf2maf.pl'
	log: '%s/{s}/08_vep.log' % (logs)
	shell:
		"""
		### Hardcoded bit to run vcf2maf and vep
		export PERL_BASE="/home/diazmiya/miniconda2/envs/vep"
		export PERL5LIB="$PERL_BASE/lib/site_perl/5.26.2/:$PERL_BASE/lib/site_perl/5.26.2/x86_64-linux-thread-multi:$" 
		# export PERL5LIB="$PERL_BASE/lib/site_perl/5.26.2/:$PERL_BASE/lib/site_perl/5.26.2/x86_64-linux-thread-multi:$PERL5LIB" 
		export PATH="$PERL_BASE/bin/:$PATH"
		export VEP_PATH="/home/diazmiya/miniconda2/envs/vep/share/ensembl-vep-99.2-0"
		export VEP_DATA="/mnt/data3/Juan/repository/ensembleVEP/hg19"
		export PERL5LIB=$VEP_PATH:$PERL5LIB
		export PATH=$VEP_PATH/htslib:$PATH
		
		gunzip -c {input} > {output.unzip_vcf}
		{params.vcf2maf} 2> {log} --input-vcf {output.unzip_vcf} --output-maf {output.maf} \
			--tumor-id {params.sam_id} --normal-id {params.con_id} \
		 	--vcf-tumor-id {params.sam_id} --vcf-normal-id {params.con_id} \
		 	--ref-fasta {params.g} --vep-data $VEP_DATA --vep-path $VEP_PATH \
		 	--filter-vcf {params.gnomadaf} --maf-center {params.center} 
		# --tmp-dir = input folder by default
		# --retain-info --retain-fmt # optional
		"""

rule concat_maf: 
	input: 
		maf = ['%s/%s/annotate/somatic.maf' % (derived,t) for t in tumors],
	output:
		 joint_maf = '%s/snv_calls/{m}/annotated/somatic_all.maf' % (derived),
		 joint_raw = '%s/snv_calls/{m}/annotated/somatic_raw.maf' % (derived),
	conda: '../envs/bwa-gatk.yaml'
	params:
		filt= 'PASS'
	log: '%s/joint/{m}/04_join_mafannot.log' % (logs)
	script:
		'../scripts/mergeMaf.R'










