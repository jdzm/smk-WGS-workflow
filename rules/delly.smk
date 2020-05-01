## Delly can be used with a file with regions to exclude. 
## svprops provide a set of scripts that are quite useful to quickly transform vcf output
## Lastly, Delly requires a tsv file when run for several samples at once. 

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
		bedbreaks = '%s/{s}/delly/somatic_svs.bp.bedpe' % (derived),
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
		
		tail -n +2 {output.svprops} | awk '{{print $1"\t"($2-500)"\t"($2+500)"\t"$5"L";}}' > {output.bedbreaks}
		tail -n +2 {output.svprops} | awk '{{print $3"\t"($4-500)"\t"($4+500)"\t"$5"R";}}' >> {output.bedbreaks}
		sort -k4,4 {output.bedbreaks} > {output.bedbreaks}

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

rule delly_joint_BND:
	input: 
		bam = ['%s/%s/aligned/raw.mdups.recal.bam' % (in_data, t) for t in tumors],
		control_bam = '%s/%s/aligned/raw.mdups.recal.bam' % (in_data, "WT_Ctrl")
	output:
		bcf = '%s/sv_calls/RPE_TP53/delly/raw_BND.bcf' % (derived),
		csi = temp('%s/sv_calls/RPE_TP53/delly/raw_BND.bcf.csi' % (derived)),
		vcf = '%s/sv_calls/RPE_TP53/delly/raw_BND.vcf.gz' % (derived),
		tbi = '%s/sv_calls/RPE_TP53/delly/raw_BND.vcf.gz.tbi' % (derived),
		svprops = '%s/sv_calls/RPE_TP53/delly/raw_calls.tab' % (derived),
		sampleprops = '%s/sv_calls/RPE_TP53/delly/raw_calls.sampleprops.tsv' % (derived)
	threads: 7
	priority: 10
	conda:
		'../envs/trans.yaml'
	params:
		g = genome["fasta"],
		excl = '%s/delly.excl.tsv' % (metadata), 
		intervals = genome['mutect_exclude'],
		sv_root = config["svprops"],
		run_delly = config["delly_precomp"]
	log:
		'%s/joint/RPE_TP53/02_delly_joint_BND.log' % (logs)
	shell:
		"""
		export OMP_NUM_THREADS=7
		{params.run_delly} call -x {params.intervals} -t BND -o {output.bcf} -g {params.g} \
			{input.bam} {input.control_bam} &> {log}
		
		bcftools convert -O z -o {output.vcf} {output.bcf} &>> {log}
		tabix {output.vcf}

		{params.sv_root}/svprops {output.bcf} > {output.svprops}
		{params.sv_root}/sampleprops {output.bcf} > {output.sampleprops}
		"""

rule delly_joint_BND_filter:
	input: 
		bcf_raw = '%s/sv_calls/RPE_TP53/delly/raw_BND.bcf' % (derived),
		csi_raw = '%s/sv_calls/RPE_TP53/delly/raw_BND.bcf.csi' % (derived)
	output:
		bcf = '%s/sv_calls/RPE_TP53/delly/somatic_BND.bcf' % (derived),
		csi = '%s/sv_calls/RPE_TP53/delly/somatic_BND.bcf.csi' % (derived),
		vcf = '%s/sv_calls/RPE_TP53/delly/somatic_BND.vcf.gz' % (derived),
		tbi = '%s/sv_calls/RPE_TP53/delly/somatic_BND.vcf.gz.tbi' % (derived),
		svprops = '%s/sv_calls/RPE_TP53/delly/somatic_BND.tab' % (derived),
		sampleprops = '%s/sv_calls/RPE_TP53/delly/somatic_BND.sampleprops.tsv' % (derived),
		bedbreaks = '%s/sv_calls/RPE_TP53/delly/somatic_BND.bp.bedpe' % (derived),
	threads: 2
	conda:
		'../envs/trans.yaml'
	params:
		sam_info = '%s/delly_samples_RPE.tsv' % (metadata),
		sv_root = config["svprops"]
	log:
		'%s/joint/RPE_TP53/02_delly_joint_BND_filt.log'  % (logs)
	shell:
		"""
		# no support in the matched normal and an overall confident VCF filter equal to PASS.
		# Not in place. just keeping code > -a 0.01 # require a minimum variant allele frequency of 25%
		
		delly filter -a 0.01 -p -f somatic -o {output.bcf} \
			-s {params.sam_info} {input.bcf_raw} &> {log}
		{params.sv_root}/svprops {output.bcf} > {output.svprops}
		{params.sv_root}/sampleprops {output.bcf} > {output.sampleprops}
		
		tail -n +2 {output.svprops} | awk '{{print $1"\t"($2-500)"\t"($2+500)"\t"$5"L";}}' > {output.bedbreaks}
		tail -n +2 {output.svprops} | awk '{{print $3"\t"($4-500)"\t"($4+500)"\t"$5"R";}}' >> {output.bedbreaks}
		sort -k4,4 {output.bedbreaks} > {output.bedbreaks}

		bcftools convert -O z -o {output.vcf} {output.bcf} &>> {log}
		tabix {output.vcf}
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
		bedbreaks = '%s/{s}/delly/somatic_{sv}.bp.bedpe' % (derived),
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
		
		tail -n +2 {output.svprops} | awk '{{print $1"\t"($2-500)"\t"($2+500)"\t"$5"L";}}' > {output.bedbreaks}.tmp
		tail -n +2 {output.svprops} | awk '{{print $3"\t"($4-500)"\t"($4+500)"\t"$5"R";}}' >> {output.bedbreaks}.tmp
		sort -k4,4 {output.bedbreaks}.tmp > {output.bedbreaks}
		rm {output.bedbreaks}.tmp

		bcftools convert -O z -o {output.vcf} {output.bcf} &>> {log}
		tabix {output.vcf}
		"""
