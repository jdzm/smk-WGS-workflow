###### Rules for Alignment preprocessing ######
## 1. prep_fastq: Starts by symlinking the fasta files to the input folder
## 2. bwa_sort: Aligns fasta to hg19 using bwa, sorts by name for fixing mate information, and sorts again by coordinate. Results in raw.bam
## 3. mark_duplicates: Marks and removes duplicates. Temporary file, leaves mdups_metrics.txt as stats file
## 4. recal_quals: Recalibrates base quality scores. Required step to improve variant calling later. Temp, leaves recal.table as stats
## 5. filter_low_mapq: Removes all reads below the mapping quality threshold specified in the config. Returns filt.bam ready to work. 

## All rules have a priority of 50-30 so that they are ran the first thing if there are new samples. 
## All rules require the same conda environment, conda: '../envs/bwa-gatk.yaml'
import os

rule prep_fastq:
    input:
        r1 = lambda wc: '%s/%s' % (samples[wc.s]['fastq'], samples[wc.s]['read1']),
        r2 = lambda wc: '%s/%s' % (samples[wc.s]['fastq'], samples[wc.s]['read2'])
    output:
        r1 = '%s/{s}/fastq/1.fastq.gz' % (in_data),
        r2 = '%s/{s}/fastq/2.fastq.gz' % (in_data)
    priority: 50
    run:
        os.symlink(input.r1, output.r1)
        os.symlink(input.r2, output.r2)
#RG = '@RG\\tID:A{params.replicate}\\tSM:{s}\\tLB:lib1\\tPU:run1\\tPL:ILLUMINA'
rule bwa_sort:
    input:
        r1 = '%s/{s}/fastq/1.fastq.gz' % (in_data),
        r2 = '%s/{s}/fastq/2.fastq.gz' % (in_data)
    output:
        prefix = temp ('%s/{s}/aligned/prefix.bam' % (in_data)),
        bam = '%s/{s}/aligned/raw.bam' % (in_data), 
    params:
        bwa_idx = config["ref_genome"]["bwa_idx"],
        proc_mem = config['sam_mem'], 
        replicate = lambda wc: '%s' % (samples[wc.s]['replicate']),
        label = '{s}'
    threads: config['threads']
    conda: '../envs/bwa-gatk.yaml'
    priority: 50
    log:
        bwa =    '%s/{s}/01_bwa.log' % (logs),
        to_bam = '%s/{s}/02_to_bam.log' % (logs)
    shell:
        """
        bwa mem -R "@RG\\tID:A{params.replicate}\\tSM:{params.label}\\tLB:lib1\\tPU:run1\\tPL:ILLUMINA" \
            -t {threads} {params.bwa_idx} {input.r1} {input.r2} 2> {log.bwa}\
        | samtools view -Shb -@ {threads} - \
        | samtools sort -n -@ {threads} -m {params.proc_mem} - > {output.prefix} 2> {log.to_bam}
        
        samtools fixmate {output.prefix} - 2>> {log.to_bam} \
        | samtools sort -@ {threads} -m {params.proc_mem} > {output.bam} 2>> {log.to_bam}

        """

rule mark_duplicates:
	input:
		bam = '%s/{s}/aligned/raw.bam' % (in_data),
	output:
		mdups = temp('%s/{s}/aligned/mdups.bam' % (in_data)),
		bai = temp('%s/{s}/aligned/mdups.bam.bai' % (in_data)),
		metrics = '%s/{s}/aligned/mdups_metrics.txt' % (in_data)
	conda: '../envs/bwa-gatk.yaml'
	threads: config['threads'] // 2
	priority: 40
	log:
		'%s/{s}/03_mark_dups.log' % (logs)
	shell:
		"""
		gatk MarkDuplicates \
		 -I {input.bam} \
		 -O {output.mdups} \
		 --REMOVE_DUPLICATES true \
		 -M {output.metrics} \
		 2> {log}
		samtools index {output.mdups}
		"""

rule recal_quals:
    input:
        '%s/{s}/aligned/mdups.bam' % (in_data)
    output: 
        recal_table = '%s/{s}/aligned/recal.table' % (in_data), 
        bam = temp('%s/{s}/aligned/recal.bam' % (in_data)),
        bai = temp('%s/{s}/aligned/recal.bam.bai' % (in_data))
    conda: '../envs/bwa-gatk.yaml'
    threads: config['threads'] // 2
    priority: 40
    params:
        g = config["ref_genome"]["bwa_idx"],
        dbsnp = config["ref_genome"]["dbsnp_vcf"]
    log:
        '%s/{s}/04_base_recal.log' % (logs)
    shell:
        """
        gatk BaseRecalibrator 2> {log} \
            -R {params.g} \
            -I {input} \
            --known-sites {params.dbsnp} \
            -O {output.recal_table}
        
        gatk ApplyBQSR 2>> {log} \
            -R {params.g} \
            -I {input} \
            --bqsr {output.recal_table} \
            -O {output.bam} 

        samtools index {output.bam}
        """

rule filter_low_mapq:
	input:
		'%s/{s}/aligned/recal.bam' % (in_data)
	output: 
		filt = '%s/{s}/aligned/filt.bam' % (in_data),
		bai = '%s/{s}/aligned/filt.bam.bai' % (in_data),
		flag = '%s/{s}/aligned/filt.flagstat' % (in_data)
	params: 
		map_qual = config['min_map_qual'],
		proc_mem = config['sam_mem']
	threads: config['threads'] // 2
	priority: 30
	conda: '../envs/bwa-gatk.yaml'
	log:
		'%s/{s}/04_mapq_filter.log' % (logs)
	shell: 
		"""
		samtools view -@ {threads} -h -b -q {params.map_qual} {input} > {output.filt}
		samtools index {output.filt} > {output.bai}
		samtools flagstat {output.filt} > {output.flag}
        touch {output.bai}
        """

# rule exclude_blacklist: 
# 	input: 
# 		filtered_bam
# 	output:
# 		nobl_bam
# 	threads: config["threads"]
# 	params:
# 		coordinates
# 	shell:
# 		"do things"
