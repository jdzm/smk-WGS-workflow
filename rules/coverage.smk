# Generates and plots histograms of coverage per chr, genome and contig
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

# Generates qc report using alfred qc from geargenomics
rule alfred_qc:
	input:
		bam = '%s/{s}/aligned/raw.mdups.recal.bam' % (in_data)
	output:
		qc = '%s/{s}/QC_plots/qc.tsv.gz' % (derived),
		qc_report =  '%s/{s}/QC_plots/qc_report.pdf' % (derived)
	threads: 1
	conda: '../envs/trans.yaml'
	params: 
		g = genome["fasta"]
	log: '%s/{s}/04_alfred-qc.log' % (logs)
	shell: 
		"""
		alfred qc -su -r {params.g} -o {output.qc} {input.bam} &> {log} 
		## -su flag serves to add correct information about hard clip rates. Not essential
		
		Rscript scripts/alfred_stats.R {output.qc} {output.qc_report} &>> {log}
		## code to extract the stats table in a readable format. 
		## zgrep ^ME qc.tsv.gz | cut -f 2- | datamash transpose | column -t 
		"""

# Plots genome coverage in 10kb bins
rule plot_cov:
	input:
		bam = '%s/{s}/aligned/raw.mdups.recal.bam' % (in_data)
	output:
		cov = '%s/{s}/QC_plots/cov.gz' % (derived),
		cov_pdf = '%s/{s}/QC_plots/cov.gz.pdf' % (derived),
	threads: 1
	conda: '../envs/trans.yaml'
	params: 
		g = genome["fasta"]
	log: '%s/{s}/04_alfred-qc.log' % (logs)
	shell: 
		"""
		alfred count_dna -o {output.cov} {input.bam} &>> {log} 
		Rscript scripts/alfred_cov.R {output.cov} {output.cov_pdf} &>> log 
		"""

# Gathers alfred qc stats from all experiments in the run to one table
rule gather_QC: 
	input: 
		qc_reps = ['%s/%s/QC_plots/qc.tsv.gz' % (derived, s) for s in samples]
	output:
		info_table = '%s/QC_summary/qc_info.tsv' % (derived)
	conda: '../envs/useR.yaml'
	params: 
		outdir = '%s/QC_summary/' % (derived)
	script:	
		'../scripts/qc_gather.R'


