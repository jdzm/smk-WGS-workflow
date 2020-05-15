# code to make a reduced merged_nodups.txt hic reads file
awk '{ if (($2 == "chr22") && ($6 =="chr22")) {print}}' merged_nodups.txt \
	| awk '{print $15, $2, $3, $1, $6, $7, $5, $9, $12}' > dummy.txt


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

gzip -dc somatic.vcf.gz | sed "s/_Aneu/_6w/g" | sed "s/_48Rec/_48h/g" | sed "s/RPE_Tumo1_6w/TP53_6wTumo1/g" | sed "s/RPE_Tumo2_6w/TP53_6wTumo2/g" | sed "s/RPE_Tumo3_6w/TP53_6wTumo3/g" > mod.vcf


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

### grab Metrics mdups
# TP53_{48h1,48h2,6w1,6w2} TP53_6wTumo{1,2,3}
grep "LIB" $mdups | sed "s/$/\tSAMPLE/" > dups.txt

for sam in WT_Ctrl TP53_Ctrl TP53_{48h1,48h2,6w1,6w2} TP53_6wTumo{1,2,3} TP53_20w TP53_20wTumo{1,2,3}
do
	grep "lib" input/$sam/aligned/mdups_metrics.txt | sed "s/$/\t$sam/" >> dups.txt
done

## Lumpy manual try
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


## Concordance analysis code 
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

####
### Fix issue with delly ###

BCF/VCF file has no sample genotypes!

# merge and index raw calls
delly merge ../derived_data/sv_calls/RPE_TP53/delly_joint/raw_BND.vcf.gz ../derived_data/sv_calls/RPE_TP53/delly_joint/raw_DEL.vcf.gz ../derived_data/sv_calls/RPE_TP53/delly_joint/raw_DUP.vcf.gz ../derived_data/sv_calls/RPE_TP53/delly_joint/raw_INV.vcf.gz -o ../derived_data/sv_calls/RPE_TP53/delly_joint/merged_raw.bcf 
bcftools convert -O z -o ../derived_data/sv_calls/RPE_TP53/delly_joint/merged_raw.vcf.gz ../derived_data/sv_calls/RPE_TP53/delly_joint/merged_raw.bcf
tabix ../derived_data/sv_calls/RPE_TP53/delly_joint/merged_raw.vcf.gz
/mnt/data3/Juan/tools/svprops/src/svprops ../derived_data/sv_calls/RPE_TP53/delly_joint/merged_raw.bcf > ../derived_data/sv_calls/RPE_TP53/delly_joint/merged_raw.tab

### now filter
# no support in the matched normal and an overall confident VCF filter equal to PASS.
# Not in place. just keeping code > -a 0.01 # require a minimum variant allele frequency of 25%

delly filter -a 0.01 -p -f somatic -o ../derived_data/sv_calls/RPE_TP53/delly_joint/merged_somatic.bcf                        -s ../metadata/delly_samples_RPE.tsv ../derived_data/sv_calls/RPE_TP53/delly_joint/merged_raw.bcf 
/mnt/data3/Juan/tools/svprops/src/svprops ../derived_data/sv_calls/RPE_TP53/delly_joint/merged_somatic.bcf > ../derived_data/sv_calls/RPE_TP53/delly_joint/merged_somatic.tab

bcftools convert -O z -o ../derived_data/sv_calls/RPE_TP53/delly_joint/merged_somatic.vcf.gz ../derived_data/sv_calls/RPE_TP53/delly_joint/merged_somatic.bcf 
tabix ../derived_data/sv_calls/RPE_TP53/delly_joint/merged_raw.vcf.gz
tabix ../derived_data/sv_calls/RPE_TP53/delly_joint/merged_somatic.vcf.gz
# rm ../derived_data/sv_calls/RPE_TP53/delly_joint/merged_somatic.bcf.csi


destdir=/Users/diazmiyar/Desktop/copy_number_profiles_all
source_path=/Users/diazmiyar/mnt/20200103-WGS/derived_data
cd $source_path
for sample in `ls -d K562* TP53* WT_Ctrl`
do
	for type in single control
	do
		file_name=${source_path}/${sample}/freec_${type}/freec_ratio.png
		if [ -e ${file_name} ] 
		then
			cp ${file_name} ${destdir}/${type}_run/${sample}_freec_ratio.png
		fi
	done
done





