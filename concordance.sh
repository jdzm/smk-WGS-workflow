# concordance analysis
## take a set of mutect calls called with the different versions of the pipeline
## Need to check how much results overlap and what is more reliable 

workdir=/mnt/data3/Juan/common-scripts/tests/pipeline-version-concordance
cd $workdir
####
file1=./mutect_6wTumo1_v1_excl/somatic_PASS.vcf # version one (hg19) of pipeline with blacklist excluded

file2=./mutect_6wTumo1_v1_noexcl/somatic_PASS_excl.vcf # version one (hg19) of pipeline with blacklist excluded in post 
file3=./mutect_6wTumo1_v1_noexcl/somatic_PASS.vcf # version one (hg19) of pipeline 

file4=./mutect_6wTumo1_v2/somatic_PASS.vcf # version two (hs37d5) of the pipeline with bl excluded 
##

### Consider liftover
### USeless. go to R

## code to check if two files are equal
cmp --silent $file1 $file2 && echo '### SUCCESS: Files Are Identical! ###' || echo '### WARNING: Files Are Different! ###'
##

### INTERSECTION analysis on pipeline V1
#project=20200103-WGS
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
