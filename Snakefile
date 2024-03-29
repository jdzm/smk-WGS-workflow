configfile: "config.yaml"

##### Helper variables #####
# Directories
in_data = config["in_data"]
derived = config["derived_data"]
metadata = config["metadata"]
logs = config["logs"]

# Sample definitions, Genome, SVs
samples = config["samples"]
# merged_samples = config["merged"]
controls = ["WT_Ctrl", "test_control", "NA", "CPA_Spike_Ctrl", "K562_Ctrl"]
tumors = list (set (samples) - set (controls))
models = list(set([samples[sample]['model'] for sample in tumors]))

genome = config["ref_genome"]["human_g1k_hs37d5"]
CHRS = genome["chromnames"]
svtypes= config["svtypes"]

##### Rules
include: "rules/align.smk" 
include: "rules/coverage.smk"
include: "rules/freec-control.smk"
include: "rules/mutect2-bychr.smk"
include: "rules/annotate.smk"
# include: "rules/delly.smk"


##### Define target files
qcplots=['cov.gen','qc_report','cov.gz']

		## Alignments
		# ['%s/%s/aligned/prefix.bam' % (in_data, s) for s in samples],
		## QC files
		# expand(['%s/%s/QC_plots/{qcplot}.pdf' % (derived, s) for s in samples], qcplot=qcplots),
		# ['%s/QC_summary/qc_info.tsv' % (derived)],

		## Delly files
		# expand(['%s/%s/delly/split_svs/raw_{sv}.vcf.gz' % (derived, t) for t in tumors], sv=svtypes),
		# ['%s/%s/delly/merged_somatic.vcf.gz' % (derived, t) for t in tumors],
		# ['%s/sv_calls/%s/delly_joint/raw_BND.vcf.gz' % (derived,m) for m in models],

		## Mutect results
		# expand(['%s/%s/mutect_calls/chrs/chr{chr}_somatic.vcf.gz' % (derived, t) for t in tumors], chr=CHRS),
		# ['%s/snv_calls/%s/mutect2-joint/somatic.vcf.gz' % (derived, m) for m in models],
		# ['%s/%s/annotate/somatic_func.maf' % (derived, t) for t in tumors],
		# ['%s/%s/mutect_calls/somatic.vcf.gz' % (derived, t) for t in tumors],
		# ['%s/%s/annotate/somatic.maf' % (derived, t) for t in tumors],
		# ['%s/snv_calls/%s/annotated/somatic_all.maf' % (derived, m) for m in models],

		## CNVs control-FREEC
		# ['%s/%s/freec_single/freec_CNVs.p.value.txt' % (derived,s) for s in samples],
		# expand(['%s/%s/freec_{runmode}_{resol}/freec_CNVs.p.value.txt' % (derived,t) for t in tumors], runmode = ['control'], resol= ['10k','50k']),
		# expand(['%s/cnv_calls/%s/control-freec/cnvs_{runmode}.tsv' % (derived,m) for m in models], runmode=['single', 'control']),
		# expand(['%s/%s/freec_{runmode}_{resol}/freec_CNVs.p.value.txt' % (derived,t) for t in tumors], runmode = ['control'], resol= ['10k','50k']),
		# expand(['%s/cnv_calls/%s/control-freec/cnvs_{runmode}_{resol}.tsv' % (derived,m) for m in models], runmode=['single', 'control'], resol= ['10k','50k']),		

rule all:
	input:
		expand(['%s/%s/QC_plots/{qcplot}.pdf' % (derived, s) for s in samples], qcplot=qcplots),
		# ['%s/QC_summary/qc_info.tsv' % (derived)],
		expand(['%s/%s/freec_{runmode}_{resol}/freec_CNVs.p.value.txt' % (derived,s) for s in samples], runmode = ['single'], resol= ['10k','50k']),
		expand(['%s/%s/freec_{runmode}_{resol}/freec_CNVs.p.value.txt' % (derived,t) for t in tumors], runmode = ['control'], resol= ['10k','50k']),
		# ['%s/%s/mutect_calls/somatic.vcf.gz' % (derived, t) for t in tumors],
		# ['%s/%s/annotate/somatic.maf' % (derived, t) for t in tumors],


