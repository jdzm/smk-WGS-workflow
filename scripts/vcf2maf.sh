## Adapted code from Daniele
conda activate vep
# homemade version
# perl is here: ./lib/site_perl/5.26.2/
#VepRefFasta="/mnt/data3/Juan/repository/ensembleVEP/hg19/homo_sapiens/99_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz"
export PERL_BASE="/home/diazmiya/miniconda2/envs/vep"
export PERL5LIB="$PERL_BASE/lib/site_perl/5.26.2/:$PERL_BASE/lib/site_perl/5.26.2/x86_64-linux-thread-multi:$PERL5LIB"
export PATH="$PERL_BASE/bin/:$PATH"
export VEP_PATH="/home/diazmiya/miniconda2/envs/vep/share/ensembl-vep-99.2-0"
export VEP_DATA="/mnt/data3/Juan/repository/ensembleVEP/hg19"
export PERL5LIB=$VEP_PATH:$PERL5LIB
export PATH=$VEP_PATH/htslib:$PATH
#export PATH=$VEP_PATH/samtools/bin:$PATH

dbsnp=/mnt/data3/Juan/repository/hg19/dbsnp/dbsnp_chr.vcf.gz
gnomad=/mnt/data3/Juan/repository/hg19/gnomAD/af-only-gnomad.raw.sites.hg19.vcf.gz
genome=/mnt/data3/Juan/repository/hg19/hg19_chr1_22XYM.fa
VER=99
# vep_install -a cf -s homo_sapiens -y GRCh37 -c $VEP_DATA -d $VEP_PATH --NO_HTSLIB -v $VER
# vep_convert_cache --species homo_sapiens --version 99_GRCh37 --dir $VEP_DATA

derived=/mnt/data3/Juan/20200103-WGS/derived_data/
logs=/mnt/data3/Juan/20200103-WGS/logs/
cd $derived
# filter_vcf=${dbsnp}
filter_vcf=${gnomad}
control=WT_Ctrl
for sam in TP53_Ctrl
do
	workdir=./$sam/annotate/
	mkdir -p $workdir
	vcf=./$sam/mutect_calls/somatic.vcf.gz
	maf=./$sam/annotate/somatic_vep_gnomad.maf
	log=$logs/$sam/09_annotate.log

	vcf_file=${vcf%.gz}
	gunzip -c ${vcf_file}".gz" > ${vcf_file}
	vcf2maf.pl 2> $log --input-vcf ${vcf_file} --output-maf ${maf} \
		--tumor-id $sam --normal-id ${control} \
	 	--vcf-tumor-id $sam --vcf-normal-id ${control} \
	 	--ref-fasta $genome --vep-data $VEP_DATA --vep-path $VEP_PATH \
	 	--filter-vcf $filter_vcf --maf-center ISREC --tmp-dir ${workdir} \
	 	--retain-info --retain-fmt
done




vcf_file=${MutectDir}"/"${normalid}"_vs_"${tumorid}"_FilterContaminFfpeOxo_PON.vcf" # R1700054-1-A_vs_R1700054-1-B_FilterContaminFfpeOxo.vcf.gz       --ref-fasta ${VepRefFasta} 
gunzip -c ${vcf_file}".gz" > ${vcf_file}
maf_file=${Mafsdir}"/"${normalid}"_vs_"${tumorid}"_FilterContaminFfpeOxo_PON.maf"
perl ${ToolsDir}"vcf2maf-1.6.16/"vcf2maf.pl --tmp-dir ${Mafsdir} --input-vcf ${vcf_file} --output-maf ${maf_file} --ncbi-build "GRCh38" --ref-fasta ${VepRefFasta} --filter-vcf ${VepFilterVcf} --vep-path $VEP_PATH --vep-data $VEP_DATA --retain-info --custom-enst ${ResourceDir}"/isoform_overrides_at_mskcc.txt" \
											  --tumor-id ${tumorid} --normal-id ${normalid} --vcf-tumor-id ${tumorid} --vcf-normal-id ${normalid}