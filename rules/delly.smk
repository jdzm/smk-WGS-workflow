## Delly can be used with a file with regions to exclude. 
## svprops provide a set of scripts that are quite useful to quickly transform vcf output
## Lastly, Delly requires a tsv file when run for several samples at once. 

rule delly_joint_call:
	input: 
		bam = ['%s/%s/aligned/raw.mdups.recal.bam' % (in_data, t) for t in tumors],
		control_bam = '%s/%s/aligned/raw.mdups.recal.bam' % (in_data, "WT_Ctrl")
	output:
		bcf = '%s/sv_calls/{m}/delly/raw_{sv}.bcf' % (derived),
		csi = '%s/sv_calls/{m}/delly/raw_{sv}.bcf.csi' % (derived),
		vcf = '%s/sv_calls/{m}/delly/raw_{sv}.vcf.gz' % (derived),
		tbi = '%s/sv_calls/{m}/delly/raw_{sv}.vcf.gz.tbi' % (derived),
		svprops = '%s/sv_calls/{m}/delly/raw_{sv}.tab' % (derived)
	threads: 7
	priority: 10
	conda:
		'../envs/trans.yaml'
	params:
		g = genome["fasta"],
		excl = '%s/delly.excl.tsv' % (metadata), 
		intervals = genome['mutect_exclude'],
		run_delly = config["delly_precomp"],
		svtype='{sv}'
	log:
		'%s/joint/{m}/02_delly_joint_{sv}.log' % (logs)
	shell:
		"""
		export OMP_NUM_THREADS=7
		{params.run_delly} call -x {params.intervals} -t {params.svtype} -o {output.bcf} -g {params.g} \
			{input.bam} {input.control_bam} &> {log}
		
		bcftools convert -O z -o {output.vcf} {output.bcf} &>> {log}
		tabix {output.vcf}
		"""

rule delly_joint_merge:
	input: 
		raw = ['%s/sv_calls/{m}/delly/raw_%s.vcf.gz' % (derived, sv) for sv in svtypes]
	output:
		bcf_raw='%s/sv_calls/{m}/delly/raw_svs.bcf' % (derived),
		vcf_raw='%s/sv_calls/{m}/delly/raw_svs.vcf.gz' % (derived),
		bcf_filt=temp ('%s/sv_calls/{m}/delly/somatic_svs.bcf' % (derived)),
		vcf_filt='%s/sv_calls/{m}/delly/somatic_svs.vcf.gz' % (derived),
		svprops_raw = '%s/sv_calls/{m}/delly/raw_svs.tab' % (derived),
		svprops_filt = '%s/sv_calls/{m}/delly/somatic_svs.tab' % (derived),
	threads: 2
	conda:
		'../envs/trans.yaml'
	params:
		sam_info = '%s/delly_samples_RPE.tsv' % (metadata),
		sv_root = config["svprops"]
	log:
		'%s/joint/{m}/02_delly_joint_merge.log'  % (logs)
	shell:
		"""		
		# merge and index raw calls
		delly merge {input.raw} -o {output.bcf_raw} &> {log}
		bcftools convert -O z -o {output.vcf_raw} {output.bcf_raw} &>> {log}
		tabix {output.vcf_raw}
		{params.sv_root}/svprops {output.bcf_raw} > {output.svprops_raw}

		### now filter
		# no support in the matched normal and an overall confident VCF filter equal to PASS.
		# Not in place. just keeping code > -a 0.01 # require a minimum variant allele frequency of 25%
		
		delly filter -a 0.01 -p -f somatic -o {output.bcf_filt} \
			-s {params.sam_info} {output.bcf_raw} &> {log}
		{params.sv_root}/svprops {output.bcf_filt} > {output.svprops_filt}
		
		bcftools convert -O z -o {output.vcf_filt} {output.bcf_filt} &>> {log}
		tabix {output.vcf_raw}
		tabix {output.vcf_filt}
		rm {output.bcf_filt}.csi 
		"""

rule delly_svcall:
	input: 
		bam = '%s/{s}/aligned/raw.mdups.recal.bam' % (in_data), 
		control_bam = lambda wc: '%s/%s/aligned/raw.mdups.recal.bam' % (in_data, samples[wc.s]['control'])
	output:
		bcf = '%s/{s}/delly/raw_{sv}.bcf' % (derived),
		csi = '%s/{s}/delly/raw_{sv}.bcf.csi' % (derived),
		vcf = '%s/{s}/delly/raw_{sv}.vcf.gz' % (derived),
		tbi = '%s/{s}/delly/raw_{sv}.vcf.gz.tbi' % (derived),
		svprops = '%s/{s}/delly/raw_{sv}.tab' % (derived)
	threads: 2
	conda:
		'../envs/trans.yaml'
	params:
		g = genome["fasta"],
		excl = '%s/delly.excl.tsv' % (metadata), 
		intervals = genome['mutect_exclude'],
		sv_root = config["svprops"],
		run_delly = config["delly_precomp"],
		svtype='{sv}'
	log:
		'%s/{s}/06_delly_{sv}.log' % (logs)
	shell:
		"""
		export OMP_NUM_THREADS=2
		{params.run_delly} call -x {params.intervals} -t {params.svtype} -o {output.bcf} -g {params.g} \
			{input.bam} {input.control_bam} &> {log}
		
		bcftools convert -O z -o {output.vcf} {output.bcf} &>> {log}
		tabix {output.vcf}

		{params.sv_root}/svprops {output.bcf} > {output.svprops}
		"""

rule delly_svfilter:
	input: 
		bcf_raw = '%s/{s}/delly/raw_{sv}.bcf' % (derived),
		csi_raw = '%s/{s}/delly/raw_{sv}.bcf.csi' % (derived)
	output:
		bcf = temp('%s/{s}/delly/somatic_{sv}.bcf' % (derived)),
		csi = temp('%s/{s}/delly/somatic_{sv}.bcf.csi' % (derived)),
		vcf = '%s/{s}/delly/somatic_{sv}.vcf.gz' % (derived),
		tbi = '%s/{s}/delly/somatic_{sv}.vcf.gz.tbi' % (derived),
		svprops = '%s/{s}/delly/somatic_{sv}.tab' % (derived),
	threads: 2
	conda:
		'../envs/trans.yaml'
	params:
		sam_info = '%s/delly_samples_RPE.tsv' % (metadata),
		sv_root = config["svprops"]
	log:
		'%s/{s}/06_delly_filter_{sv}.log' % (logs)
	shell:
		"""
		# no support in the matched normal and an overall confident VCF filter equal to PASS.
		# Not in place. just keeping code > -a 0.01 # require a minimum variant allele frequency of 1%
		
		delly filter -a 0.01 -p -f somatic -o {output.bcf} \
			-s {params.sam_info} {input.bcf_raw} &> {log}
		{params.sv_root}/svprops {output.bcf} > {output.svprops}

		bcftools convert -O z -o {output.vcf} {output.bcf} &>> {log}
		tabix {output.vcf}
		"""

#svtypes= config['svtypes']

rule delly_svmerge:
	input: 
		raw = ['%s/{s}/delly/raw_%s.vcf.gz' % (derived, sv) for sv in svtypes],
		filt = ['%s/{s}/delly/somatic_%s.vcf.gz' % (derived, sv) for sv in svtypes],
	output:
		bcf_raw=temp ('%s/{s}/delly/raw_svs.bcf' % (derived)),
		bcf_filt=temp ('%s/{s}/delly/somatic_svs.bcf' % (derived)),
		vcf_raw='%s/{s}/delly/raw_svs.vcf.gz' % (derived),
		vcf_filt='%s/{s}/delly/somatic_svs.vcf.gz' % (derived),
	threads: 2
	conda:
		'../envs/trans.yaml'
	log:
		'%s/{s}/06_delly_merge_svs.log' % (logs)
	shell:
		"""		
		delly merge {input.filt} -o {output.bcf_filt} &> {log}
		delly merge {input.raw} -o {output.bcf_raw} &> {log}

		bcftools convert -O z -o {output.vcf_filt} {output.bcf_filt} &>> {log}
		bcftools convert -O z -o {output.vcf_raw} {output.bcf_raw} &>> {log}
		tabix {output.vcf_filt}
		tabix {output.vcf_raw}
		rm {output.bcf_raw}.csi {output.bcf_filt}.csi 
		"""


