###### Rules for Alignment preprocessing ######
## 1. prep_fastq: Starts by symlinking the fasta files to the input folder
## 2. bwa_sort: Aligns fasta to hg19 using bwa, sorts by name for fixing mate information. Results in raw.bam which will be temp
## 3. mark_duplicates: Marks duplicates. Produces mdups_metrics.txt as stats file
## 4. recal_quals: Recalibrates base quality scores. Required step to improve variant calling later. Final 'raw' bamfile
## # 5. rmdups: snippet to remove duplicates. Not used as it is not needed for variant calling 
## bamfiles output by the last rule are protected and snakemake will not overwrite them once produced

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
        bam = temp('%s/{s}/aligned/raw.nsort.bam' % (in_data))
    params:
        g = genome["fasta"],
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
            -t {threads} -M {params.g} {input.r1} {input.r2} 2> {log.bwa}\
        | samtools view -Shb -@ {threads} - \
        | samtools sort -n -@ {threads} -m {params.proc_mem} - > {output.prefix} 2> {log.to_bam}
        
        samtools fixmate -@ {threads} {output.prefix} {output.bam} 2>> {log.to_bam} 

        """

# Running an instance of MarkDuplicatesSpark. Will not write a PG line on the header
# Duplicates are not removed in this step. If you use MarkDuplicates, might have problems with sort
# sorted by coordinate automatically
rule mark_duplicates:
    input:
        bam = '%s/{s}/aligned/raw.nsort.bam' % (in_data)
    output:
        mdups = temp('%s/{s}/aligned/mdups.bam' % (in_data)),
        bai = temp('%s/{s}/aligned/mdups.bam.bai' % (in_data)),
        sbi = temp('%s/{s}/aligned/mdups.bam.sbi' % (in_data)),
        metrics = '%s/{s}/aligned/mdups_metrics.txt' % (in_data)
    conda: '../envs/bwa-gatk.yaml'
    threads: config['threads']
    priority: 40
    params: tmpdir = config["scratch"]
    log: '%s/{s}/03_mark_dups.log' % (logs)
    shell:
        """
        gatk MarkDuplicatesSpark 2> {log} \
         -I {input.bam} \
         -O {output.mdups} \
         -M {output.metrics} \
         --QUIET=true \
         --conf 'spark.executor.cores={threads}' \
         --conf 'spark.local.dir={params.tmpdir}'
        
        """
# samtools index {output.mdups}
# bai = temp('%s/{s}/aligned/mdups.bam.bai' % (in_data)),
# grep "## MET" -A 2 {output.metrics} >> {output.recordstats}
#--REMOVE_DUPLICATES true \

rule recal_quals:
    input:
        mdups = '%s/{s}/aligned/mdups.bam' % (in_data)
    output: 
        recal_table = '%s/{s}/aligned/recal.table' % (in_data), 
        bam = protected('%s/{s}/aligned/raw.mdups.recal.bam' % (in_data)),
        bai = '%s/{s}/aligned/raw.mdups.recal.bai' % (in_data)
    conda: '../envs/bwa-gatk.yaml'
    threads: config['threads'] 
    priority: 40
    params:
        g = genome["fasta"],
        dbsnp = genome["dbsnp_vcf"]
    log:
        '%s/{s}/04_base_recal.log' % (logs)
    shell:
        """
        gatk BaseRecalibrator 2> {log} \
            -R {params.g} \
            -I {input.mdups} \
            --known-sites {params.dbsnp} \
            -O {output.recal_table}
        
        gatk ApplyBQSR 2>> {log} \
            -R {params.g} \
            -I {input.mdups} \
            --bqsr {output.recal_table} \
            -O {output.bam} 

        touch {output.bai} 
        """

# rule rmdups:
#     input:
#         bam = '%s/{s}/aligned/raw.mdups.recal.bam' % (in_data),
#         recordstats = '%s/{s}/aligned/recordstats.txt' % (in_data)
#     output: 
#         nodups = protected('%s/{s}/aligned/nodups.bam' % (in_data)),
#         dups = '%s/{s}/aligned/dups.bam' % (in_data),
#         bai = '%s/{s}/aligned/nodups.bam.bai' % (in_data)
#     threads: config['threads'] // 2
#     priority: 30
#     conda: '../envs/bwa-gatk.yaml'
#     log:
#         '%s/{s}/04_mapq_filter.log' % (logs)
#     shell: 
#         """
#         samtools view -@ {threads} -h -b -F 0X400 {input.bam} > {output.nodups}
#         samtools view -@ {threads} -h -b -f 0X400 {input.bam} > {output.dups}

#         samtools index {output.nodups} > {output.bai}

#         line5=$(samtools view -c {output.nodups})
#         echo $line5 "rmdups" >> {input.recordstats}
#         """
# samtools view -@ {threads} -h -b -q {params.map_qual} {input} > {output.filt}