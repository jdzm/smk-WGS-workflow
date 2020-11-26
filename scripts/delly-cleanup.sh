# script to cleanup delly calls and output txt tables
# How to get AF
# Sure, using bcftools you first need to subset the VCF to the tumor sample then filter for precise or imprecise SVs and output the
# REF and ALT support (DR and DV for imprecise SVs and RR and RV for precise SVs).

# Precise SVs:
# bcftools view -s <tumor_vcf_sample_name> delly.sv.bcf | bcftools filter -i 'PRECISE==1' - | bcftools query -f "%CHROM %POS %INFO/CHR2 %INFO/END [ %RV %RR]\n" - 

# Imprecise SVs:
# bcftools view -s <tumor_vcf_sample_name> delly.sv.bcf | bcftools filter -i 'IMPRECISE==1' - | bcftools query -f "%CHROM %POS %INFO/CHR2 %INFO/END [ %DV %DR]\n" -


# DR and DV for imprecise SVs and RR and RV for precise SVs

# also TP53_48h and TP53_6w
for callset in TP53_6w{1,2,3} TP53_6wTumo{11,12,13,21,22,23,31,32,33}; do
	cd /mnt/data3/Juan/20200610-WGD-E02/derived_data/${callset}/delly/split_svs
	querystr='"%CHROM %POS %ID %REF %ALT %QUAL %INFO/CHR2 %INFO/END %INFO/PRECISE %INFO/IMPRECISE %INFO/SVLEN %INFO/INSLEN %INFO/HOMLEN %INFO/PE %INFO/SR [%RV %RR %DR %DV %CN %GT %GQ %GL %FT]\n"'
	for sample in ${callset} WT_Ctrl; do
		for svtype in DUP DEL INV BND; do
			in_BCF=raw_${svtype}.bcf

			CMD="bcftools view -f PASS -s $sample $in_BCF | bcftools query -f ${querystr}  - "

			eval $CMD > passing_${svtype}_${sample}.txt
		done
		mkdir -p ../merged_svs
		cat passing_{DUP,DEL,INV,BND}_${sample}.txt > ../merged_svs/merged_svs_${sample}.txt
		rm passing_*.txt
	done
done

# %CHROM %POS %ID %REF %ALT %INFO/CHR2 %INFO/CHR2 %INFO/PRECISE %INFO/IMPRECISE %INFO/INSLEN %INFO/HOMLEN
# %INFO/PE %INFO/SR %INFO/QUAL %INFO/SOMATIC
# %
#  [ %RV %RR %DR %DV %CN %GT %GQ %GL]\n

#### Extract for joint call
cd /mnt/data3/Juan/20200424-Variants/derived_data/sv_calls/RPE_TP53/delly_joint

querystr='"%CHROM %POS %ID %REF %ALT %QUAL %INFO/CHR2 %INFO/END %INFO/PRECISE %INFO/IMPRECISE \
%INFO/SVLEN %INFO/INSLEN %INFO/HOMLEN %INFO/PE %INFO/SR \
[ %RV %RR %DR %DV %CN %GT %GQ %GL %FT]\n"'
for sample in WT_Ctrl TP53_{Ctrl,48h1,48h2,6w1,6w2,20w} TP53_6wTumo{1,2,3} TP53_20wTumo{1,2,3}; do
	for svtype in DUP DEL INV BND; do
		in_BCF=raw_${svtype}.bcf

		CMD="bcftools view -f PASS -s $sample $in_BCF | bcftools query -f ${querystr}  - "

		eval $CMD > passing_${svtype}_${sample}.txt
	done
	cat passing_{DUP,DEL,INV,BND}_${sample}.txt > merged_svs_${sample}.txt
done
rm passing_*.txt 
### Problems with passing in INV DUP DEL for TP53_48h2 and TP53_6w2





