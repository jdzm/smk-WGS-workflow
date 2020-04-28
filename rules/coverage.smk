rule coverage:
	input:
		bam = '%s/{s}/aligned/raw.mdups.recal.bam' % (in_data)
	output: 
		hist = '%s/{s}/QC_plots/cov.hist' % (derived), 
		plot_cov1 = '%s/{s}/QC_plots/cov.gen.pdf' % (derived), 
		plot_cov2 = '%s/{s}/QC_plots/cov.chr.pdf' % (derived) 
	threads: 1
	params: 
		sam_id = '{s}',
		chromsizes = genome['chrom_sizes']['full']
	conda: 
		'../envs/freec.yaml'
	shell: 
		"""		
		bedtools genomecov -ibam {input.bam} -max 40 > {output.hist}

		Rscript --vanilla scripts/plotCov.R {params.sam_id} {params.chromsizes} {output.hist} 
        """

rule alfred_qc:
	input:
		bam = '%s/{s}/aligned/raw.mdups.recal.bam' % (in_data)
	output:
		qc = '%s/{s}/QC_plots/qc.tsv.gz' % (derived),
		cov = '%s/{s}/QC_plots/cov.gz' % (derived)
	threads: config['threads'] // 2
	conda: '../envs/trans.yaml'
	params: 
		g = genome["fasta"]
	log: '%s/{s}/04_alfred-qc.log' % (logs)
	shell: 
		"""
		alfred qc -r {params.g} -o {output.qc} {input.bam} &> {log}
		alfred count_dna -o {output.cov} {input.bam} &>> {log}
		"""

rule plot_qc: 
	input:
		qc = '%s/{s}/QC_plots/qc.tsv.gz' % (derived)
	output:
		qc_report =  '%s/{s}/QC_plots/qc_report.pdf' % (derived)
	threads: 1
	conda: '../envs/trans.yaml'
	log: '%s/{s}/04_alfred-qc.log' % (logs)
	shell: 
		"""
		Rscript scripts/alfred_stats.R {input.qc} {output.qc_report} &>> {log}
		# zgrep ^ME qc.tsv.gz | cut -f 2- | datamash transpose | column -t^C 
		# # code to extract the stats table in a readable format. 
		# Rscript scripts/alfredrd input.cov output.cov_pdf &>> log # not necessary
		"""



