# code to make a reduced merged_nodups.txt hic reads file
awk '{ if (($2 == "chr22") && ($6 =="chr22")) {print}}' merged_nodups.txt \
	| awk '{print $15, $2, $3, $1, $6, $7, $5, $9, $12}' > dummy.txt

# code to loop over CNV calls, subtract and count lines
control=WT_Ctrl/freec_single/freec_CNVs
for sample in `ls -1d TP53*` WT_Ctrl 
	do 
		bedfile=${sample}/freec_single/freec_CNVs 
		wc -l $bedfile
		bedtools subtract -a $control -b $bedfile | wc -l
done


## DELLY Merged calls and filtering
exclude=/mnt/data3/Juan/tools/delly/excludeTemplates/human.hg19.excl.tsv
genome=/mnt/data3/Juan/repository/hg19/hg19_chr1_22XYM.fa
############ K562 ############
sampleinfo=/mnt/data3/Juan/20200103-WGS/metadata/delly_samples_K562.tsv
logfile=delly.log
delly call -x $exclude -o svs_delly.bcf -g $genome \
	../../../input/K562_{Treated1,Treated2,Control}/aligned/filt.bam \
	&> $logfile

#bcftools view svs_delly.bcf > svs_delly.vcf
$svprops svs_delly.bcf > svs_delly.tab
$sampleprops svs_delly.bcf > svs_delly.sampleprops.tsv
#There are many parameters available to tune the somatic SV filtering. 
# Below we require 
# no support in the matched normal and an overall confident 
# SV site prediction with the VCF filter field being equal to PASS.
# -a 0.25 # equire a minimum variant allele frequency of 25%
delly filter -p -f somatic -o somatic_delly.bcf \
	-s $sampleinfo svs_delly.bcf &>> $logfile
$sampleprops somatic_delly.bcf > somatic_delly.sampleprops.tsv
$svprops somatic_delly.bcf > somatic_delly.tab

############ RPE ############
sampleinfo=/mnt/data3/Juan/20200103-WGS/metadata/delly_samples_RPE.tsv
logfile=delly.log
delly call -x $exclude -o svs_delly.bcf -g $genome \
	../../../input/TP53_{Ctrl,48Rec1,48Rec2,Aneu1,Aneu2}/aligned/filt.bam \
	../../../input/WT_Ctrl/aligned/filt.bam &> $logfile

$svprops svs_delly.bcf > svs_delly.tab
$sampleprops svs_delly.bcf > svs_delly.sampleprops.tsv

# no support in the matched normal and an overall confident VCF filter equal to PASS.
# NO> -a 0.25 # equire a minimum variant allele frequency of 25%
delly filter -p -f somatic -o somatic_delly.bcf \
	-s $sampleinfo svs_delly.bcf &>> $logfile
$sampleprops somatic_delly.bcf > somatic_delly.sampleprops.tsv
$svprops somatic_delly.bcf > somatic_delly.tab


#########

cd derived_data

for sample in `ls WT_* TP53_* K562_* -d`
for sample in `ls TP53_20w* -d`
do
	echo $sample
	zcat $sample/alfred_qc/qc.tsv.gz | grep ^ME | datamash transpose | column -t | grep Median
	zcat $sample/alfred_qc/qc.tsv.gz | grep ^ME | datamash transpose | column -t | grep "#MappedPairs"
done
#zcat $sample/alfred_qc/qc.tsv.gz | grep ^ME | datamash transpose | column -t | grep MedianInsertSize


# Code to extract qc table. Need to optimize. 
cd ../input

for sample in `ls WT_* TP53_* RPE_* K562_* -d`
for sample in `ls TP53_20w* -d`
do
	raw=`awk 'NR==1{print $1}' $sample/aligned/raw.flagstat`
	filt=`awk 'NR==1{print $1}' $sample/aligned/filt.flagstat`
	echo "sample" $sample "TotalRawReads" $raw "TotalFiltReads" $filt 
done

for sample in `ls WT_* TP53_* RPE_* K562_* -d`
for sample in `ls TP53_20w* -d`
do
	if [ ! -f $sample/aligned/raw.flagstat ]; then
		echo "process" $sample
	    samtools flagstat -@ 15 $sample/aligned/raw.bam > $sample/aligned/raw.flagstat
	fi
done


########
$svprops somatic_svs.bcf | tail -n +2 | awk '{print $1"\t"($2-500)"\t"($2+500)"\t"$5"L";}' > somatic_svs.bp.bed
$svprops somatic_svs.bcf | tail -n +2 | awk '{print $3"\t"($4-500)"\t"($4+500)"\t"$5"R";}' >>  somatic_svs.bp.bed
sort -k4,4 somatic_svs.bp.bed

########
# Replace read groups 
samtools view -H $sample/aligned/filt.bam | grep '@RG' | awk 'FNR == 1 {print $3}' | cut -d':' -f 2
# need to replace SM with the new name of the sample
# RG = '@RG\\tID:A1\\tSM:{s}\\tLB:lib1\\tPU:run1\\tPL:ILLUMINA'
oldname=$(samtools view -H $sample/aligned/filt.bam | grep '@RG' | awk 'FNR == 1 {print $3}' | cut -d':' -f 2)


oldname=test
bamfile=filt.bam
newname=sprinkles
newfile=filt_RG.bam


gatk AddOrReplaceReadGroups -I $bamfile -O $newfile \
	-ID A1 -LB lib1 -PU run1 -PL ILLUMINA \
	-SM $newname
samtools index $newfile


### gatk way, it makes the file 10Gb bigger
# gatk AddOrReplaceReadGroups &>> ${logfile} -I ${sample}/aligned/filt.bam -O ${sample}/aligned/filt_RG.bam \
# 	-ID A1 -SM ${sample} -LB lib1 -PU run1 -PL ILLUMINA 


for sample in TP53_{48h2,6w1,6w2} TP53_6wTumo{1,2,3}
do
	oldname=$(samtools view -H $sample/aligned/filt.bam | grep '@RG' | awk 'FNR == 1 {print $3}' | cut -d':' -f 2) 
	if [ ${sample} != ${oldname} ] 
	then
		echo ${sample} ${oldname} 
		#samtools view -H $sample/aligned/filt.bam | sed -z "s/${oldname}/${sample}/" | samtools reheader - $sample/aligned/filt.bam > $sample/aligned/filt_RG.bam 
		mv $sample/aligned/filt_RG.bam $sample/aligned/filt.bam
		samtools index $sample/aligned/filt.bam > $sample/aligned/filt.bam.bai
	fi
done

# do some touching
for s in TP53_{48h1,48h2,6w1,6w2} TP53_6wTumo{1,2,3}
do
	touch ../derived_data/${s}/
	touch ../derived_data/${s}/coverage/cov.hist
	touch ../derived_data/${s}/coverage/cov.gen.pdf
	touch ../derived_data/${s}/coverage/cov.chr.pdf
	# 
	touch ../derived_data/${s}/alfred_qc/qc.tsv.gz
	touch ../derived_data/${s}/alfred_qc/cov.gz
	touch ../derived_data/${s}/alfred_qc/qc_report.pdf
	
	for m in control single
	do
		touch ../derived_data/${s}/freec_${m}/config_${m}.txt
		touch ../derived_data/${s}/freec_${m}/freec_CNVs
		touch ../derived_data/${s}/freec_${m}/freec_ratio.txt
		touch ../derived_data/${s}/freec_${m}/freec_ratio.png
		touch ../derived_data/${s}/freec_${m}/freec_CNVs.p.value.txt
	done
done

touch ../derived_data/snv_calls/RPE_TP53/mutect2/somatic.vcf.gz
touch ../derived_data/snv_calls/RPE_TP53/mutect2/somatic.vcf.gz.tbi

gzip -dc somatic.vcf.gz | sed "s/_Aneu/_6w/g" | sed "s/_48Rec/_48h/g" | sed "s/RPE_Tumo1_6w/TP53_6wTumo1/g" | sed "s/RPE_Tumo2_6w/TP53_6wTumo2/g" | sed "s/RPE_Tumo3_6w/TP53_6wTumo3/g" > mod.vcf

#####

genome=/mnt/data3/Juan/repository/hg19/hg19_chr1_22XYM.fa
chromsizes=/mnt/data3/Juan/repository/hg19/hg19.main.chrom.sizes
cd /mnt/data3/Juan/igv-files/WGS

for sample in TP53_Ctrl TP53_6w1
do
	#igvtools toTDF $WGS/input/$sample/aligned/filt.bam ${sample}.filt.bam.tdf ${chromsizes}
	igvtools count -f min,max,mean $WGS/input/$sample/aligned/filt.bam ${sample}.filt.bam.tdf ${chromsizes}
done


SURVIVOR vcftobed somatic_svs.vcf 0 -1 somatic_svs.bedpe

#### FUNCOTATOR
vars=/mnt/data3/Juan/20200103-WGS/derived_data/snv_calls/RPE_TP53/mutect2/somatic.vcf.gz
genome=/mnt/data3/Juan/repository/hg19/hg19_chr1_22XYM.fa
func_dir=/mnt/data3/Juan/repository/funcotator_v1.6
funvars=/mnt/data3/Juan/20200103-WGS/derived_data/snv_calls/RPE_TP53/funcotator/variants.funcotated.maf

gatk Funcotator 2>local-maf.log --variant $vars --reference $genome --ref-version hg19 \
	--data-sources-path $func_dir --output $funvars --output-file-format MAF \
	--remove-filtered-variants true 



gatk VariantsToTable \
     -V filtered_PASS.vcf \
     -F CHROM -F POS -F TYPE -GF AD \
     -O output.table

gatk VariantsToTable \
     -V filtered_PASS.vcf \
     -F HOM-REF -F HET -F HOM-VAR -F NO-CALL \
     -O output2.table

gatk VariantsToTable \
	 --QUIET \
     -V filtered_PASS.vcf \
     -F CHROM -F POS -F TYPE -GF AF \
     -O output3.table


## bcftools
bcftools view -Ou -s WT_Ctrl,TP53_Ctrl filtered_PASS.vcf | bcftools query -f %INFO/AC\t%INFO/AN\n

bcftools view -f PASS somatic.vcf.gz > filtered_PASS.vcf

bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\n" somatic_svs.vcf

bcftools stats -s - filtered_PASS.vcf > filtered_PASS.vcf.stats



### MEtrics mdups

# TP53_{48h1,48h2,6w1,6w2} TP53_6wTumo{1,2,3}
grep "LIB" $mdups | sed "s/$/\tSAMPLE/" > dups.txt

for sam in WT_Ctrl TP53_Ctrl TP53_{48h1,48h2,6w1,6w2} TP53_6wTumo{1,2,3} TP53_20w TP53_20wTumo{1,2,3}
do
	grep "lib" input/$sam/aligned/mdups_metrics.txt | sed "s/$/\t$sam/" >> dups.txt
done


# recordstats
line1=$(samtools view -c {output.prefix})
line2=$(samtools view -c {output.bam})
echo $line1 "raw entries" > {output.recordstats}
echo $line2 "fixmate entries" >> {output.recordstats}
line3=$(samtools view -c {output.mdups})
echo $line3 "mdups entries" >> {input.recordstats}
line4=$(samtools view -c {output.bam})
echo $line4 "recal entries" >> {input.recordstats}
line5=$(samtools view -c -F 0X400 {output.bam})
echo $line5 "rmdups" >> {input.recordstats}

##
conda activate lumpy
bamfile=raw.mdups.recal.bam
bamfile=raw.bam
threads=15

samtools view -b -F 1294 $bamfile \
			| samtools sort -@ $threads - > discordants.bam

samtools view -h $bamfile \
		    | extractSplitReads_BwaMem -i stdin -d False \
		    | samtools view -@ $threads -Sb - \
		    | samtools sort -@ $threads - \
		    > splitters.bam


fasta=/mnt/data3/Juan/repository/human_g1k_hs37d5/human_g1k_hs37d5.fasta
excl=/mnt/data3/Juan/tools/delly/excludeTemplates/human.human_g1k_hs37d5.excl.tsv
threads=1
smoove call -x --genotype --name $name --outdir . \
	-f $fasta --processes 12 --exclude $bed *.bam


smoove call --outdir results-smoove/ --exclude $excl --name test_tumor \
	--fasta $fasta -p $threads --genotype raw.mdups.recal.bam

lumpyexpress \
	-B ../$bamfile \
	-S test_tumor.split.bam \
	-D test_tumor.disc.bam


for sample in `ls -1 input/`
do
	qcfile=derived_data/$sample/QC_plots/qc.tsv.gz
	qcinfo=derived_data/$sample/QC_plots/qcinfo.tsv
	if [ -f "$qcfile" ]; then
		zgrep ^ME $qcfile | cut -f 2- | datamash transpose | awk 'OFS="\t" {print $1,$2}' > $qcinfo
	else
		echo $qcfile "does not exist"
	fi
done

cd $WGS
for sample in `ls -1 input/`
do
	qcfile=derived_data/$sample/alfred_qc/qc.tsv.gz
	qcinfo=derived_data/$sample/alfred_qc/qcinfo.tsv
	if [ -f "$qcfile" ]; then
		zgrep ^ME $qcfile | cut -f 2- | datamash transpose | awk 'OFS="\t" {print $1,$2}' > $qcinfo
	else
		echo $qcfile "does not exist"
	fi
done


qcfiles=$(find ./ -name "qcinfo.tsv")

# Alfredd to evaluate suppl and secondary

fasta=/mnt/data3/Juan/repository/human_g1k_hs37d5/human_g1k_hs37d5.fasta
bam1
qcrep1

alfred qc -su -r $fasta $bam1 -o $qcrep1 


 write for loop for samples

exclude=/mnt/data3/Juan/repository/hg19/hg19.exclude.bed
project=20200103-WGS
sample=TP53_6w1
vcf_raw=/mnt/data3/Juan/$project/derived_data/$sample/mutect_calls/somatic.vcf
vcf_PASS=somatic_PASS.vcf
bcftools view -f PASS $vcf_raw > $vcf_PASS

lines_raw=$(zgrep -v ^# $vcf_raw | wc -l )
lines_pass=$(grep -v ^#  $vcf_PASS| wc -l )

#bedtools subtract -a $vcf_PASS -b $exclude > test.bs.pass.excl.vcf
bedtools intersect -header -v -a $vcf_PASS -b $exclude > test.is.pass.excl.vcf
#bedtools intersect -v -a $vcf_PASS -b $exclude > test.is.pass.excl.bed # same, results, use intersect -header if want to keep header

wc -l test.is.pass.excl.bed
grep -v "#" test.bs.pass.excl.vcf | wc -l


bedtools intersect -v -a test.bs.pass.excl.vcf -b test.is.pass.excl.vcf

cmp --silent test.bs.pass.excl.vcf test.is.pass.excl.vcf && echo '### SUCCESS: Files Are Identical! ###' || echo '### WARNING: Files Are Different! ###'
cmp --silent test.bs.pass.excl.bed test.is.pass.excl.bed && echo '### SUCCESS: Files Are Identical! ###' || echo '### WARNING: Files Are Different! ###'

### Loop to generate files
project=20200103-WGS
#project=20200424-Variants
statsfile=/mnt/data3/Juan/${project}/derived_data/snv_calls/RPE_TP53/indiv_calls.txt
statsfile=/mnt/data3/Juan/${project}/indiv_calls.txt
touch $statsfile

# TP53_Ctrl TP53_{48h1,6w1,6wTumo1,20w,20wTumo1}
for sample in TP53_6wTumo1
do
	# vcf_raw=/mnt/data3/Juan/$project/derived_data/$sample/mutect_calls/somatic.vcf.gz
	# vcf_PASS=/mnt/data3/Juan/$project/derived_data/$sample/mutect_calls/somatic_PASS.vcf
	# vcf_filt=/mnt/data3/Juan/$project/derived_data/$sample/mutect_calls/somatic_PASS_excl.vcf
	vcf_raw=/mnt/data3/Juan/$project/derived_data/$sample/mutect_backup_noexcl/somatic.vcf.gz
	vcf_PASS=/mnt/data3/Juan/$project/derived_data/$sample/mutect_backup_noexcl/somatic_PASS.vcf
	vcf_filt=/mnt/data3/Juan/$project/derived_data/$sample/mutect_backup_noexcl/somatic_PASS_excl.vcf
	bcftools view -f PASS $vcf_raw > $vcf_PASS
	bedtools intersect -header -v -a $vcf_PASS -b $exclude > $vcf_filt
	echo $sample "raw" $(zgrep -v ^# $vcf_raw | wc -l ) >> $statsfile
	echo $sample "PASS" $(grep -v ^#  $vcf_PASS| wc -l ) >> $statsfile
	echo $sample "PASS_excl" $(grep -v ^#  $vcf_filt| wc -l ) >> $statsfile
done

####
file1=./derived_data/TP53_Ctrl/mutect_calls/somatic_PASS_excl.vcf
file2=./derived_data/TP53_6w1/mutect_calls/somatic_PASS_excl.vcf
file2=./derived_data/TP53_20w/mutect_calls/somatic_PASS_excl.vcf

project=20200103-WGS
#filtered_vcfs=$(find -name 'somatic_PASS_excl.vcf')
#sample_names=$(ls -1 $filtered_vcfs | cut -d/ -f3)

intersect_vcfs=/mnt/data3/Juan/${project}/derived_data/snv_calls/RPE_TP53/vcf_intersect_all.bed

## !!! Watch out order!!! 
sample_names=$(echo TP53_{Ctrl,48h1,6w1,6wTumo1,20w,20wTumo1})

filtered_vcfs="./derived_data/TP53_Ctrl/mutect_calls/somatic_PASS_excl.vcf ./derived_data/TP53_48h1/mutect_calls/somatic_PASS_excl.vcf \
./derived_data/TP53_6w1/mutect_calls/somatic_PASS_excl.vcf ./derived_data/TP53_6wTumo1/mutect_calls/somatic_PASS_excl.vcf \
./derived_data/TP53_20w/mutect_calls/somatic_PASS_excl.vcf ./derived_data/TP53_20wTumo1/mutect_calls/somatic_PASS_excl.vcf"

multiIntersectBed -i $filtered_vcfs -header -names $sample_names > $intersect_vcfs
CMD="multiIntersectBed -i $filtered_vcfs -header -names $sample_names > $intesect_vcfs"




bedtools jaccard -a $file1 -b $file2 > jacc.test

bedfiles=$(find -name "somatic_PASS_excl.vcf")

bedtools multiIntersectBed -i $filtered_vcfs




