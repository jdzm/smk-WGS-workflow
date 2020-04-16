#!/bin/bash
# conda activate align

projectDir=/mnt/data3/Juan/20200103-WGS
cd $projectDir
files=$(find -name *_filtered.bam | tr '\n' ' ')
labs='TP53_48Rec1 TP53_Aneu2 TP53_Aneu1 TP53_48Rec2 TP53_Ctrl WT_Ctrl'
blacklist=metadata/hg19/hg19-blacklist.v2.bed
threads=16

multiBamSummary bins --bamfiles ${files} \
	--outFileName derived_data/multiBamResults.npz \
	--smartLabels \
	--binSize 10000 \
	--numberOfProcessors ${threads} \
	--blackListFileName $blacklist \
	--outRawCounts derived_data/multiBamResults_rawCounts.tsv \
	2> derived_data/multiBamResults.log


# plot spearman correlation
plotCorrelation -in derived_data/multiBamResults.npz \
	--corMethod spearman --whatToPlot heatmap \
	--plotNumbers \
	--colorMap bwr \
	--plotFile derived_data/coverage_sp_heat.pdf \
	--outFileCorMatrix derived_data/coverage_corrMat_sp.tsv

# plot pearson correlation 
plotCorrelation -in derived_data/multiBamResults.npz \
	--corMethod pearson --whatToPlot heatmap \
	--plotNumbers \
	--colorMap bwr \
	--plotFile derived_data/coverage_pe_heat.pdf \
	--outFileCorMatrix derived_data/coverage_corrMat_pe.tsv

##### Plot Cover
# the tool samples one milion bp, counts the n of overlapping reads and reports how many bases are covered hoe many times
out=derived_data/coverage
for chr in chr{1..22} chrX
do 
	plotCoverage -b $files -o ${out}/multiCoveragePlot_${chr}.pdf \
		--smartLabels --numberOfSamples 1000000 \
		--region ${chr} \
		--blackListFileName $blacklist \
		--numberOfProcessors 15 \
		--outRawCounts ${out}/multiCoverage_rawcounts${chr}.tsv \
		2> ${out}/multiCoveragePlot_{chr}.log
done