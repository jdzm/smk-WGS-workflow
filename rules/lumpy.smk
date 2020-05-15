rule lumpy_prep: 
	input:
		bam = '%s/{s}/aligned/raw.mdups.recal.bam' % (in_data), 
	output:
		disc = '%s/{s}/aligned/discordants.bam' % (in_data),
		split = '%s/{s}/aligned/splitters.bam' % (in_data), 
	threads: config['threads'] 
	conda:
	log:
	shell:
		"""
		# Extract the discordant paired-end alignments.
		samtools view -b -F 1294 {input.bam} \
			| samtools sort -@ {threads} - > {output.disc}

		# Extract the split-read alignments
		samtools view -h {input} \
		    | extractSplitReads_BwaMem -i stdin -d False \
		    | samtools view -Sb - \
		    | samtools sort -@ {threads} - \
		    > {output.split}
		"""



fasta=/mnt/data3/Juan/repository/human_g1k_hs37d5/human_g1k_hs37d5.fasta
excl=/mnt/data3/Juan/tools/delly/excludeTemplates/human.hg19.excl.tsv

smoove call -x --genotype --name $name --outdir . \
	-f $fasta --processes 12 --exclude $bed *.bam


smoove call --outdir results-smoove/ --exclude $bed --name $sample \
	--fasta $fasta -p 1 --genotype /path/to/$sample.bam