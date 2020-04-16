# merge bams if necessary
## define set of samples to merge
## 

merge_bams:
	input: expand ('%s/{m}/aligned.filt.bam' % (input), m=[wc.to_merge])
	output: 
		bam = '%s/{to_merge}/aligned/filt.bam' % (input), 
		bai = '%s/{to_merge}/aligned/filt.bam.bai' % (input)
	conda: '../envs/bwa-gatk.yaml'
	threads: config['threads'] // 4
	shell:
		"""
		samtools merge {output.bam} {input}
		samtools index {output.bam}
		# gatk addOrReplaceReadGroups
		## reag groups have old names as samples, need merged names
		"""