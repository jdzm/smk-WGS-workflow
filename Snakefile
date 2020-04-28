configfile: "config.yaml"

#### Helper variables

samples = config["samples"]
controls = ["WT_Ctrl", "test_control"]
tumors = list (set (samples) - set (controls))
in_data = config["in_data"]
derived = config["derived_data"]
metadata = config["metadata"]
genome = config["ref_genome"]["human_g1k_hs37d5"]
logs = config["logs"]

models = list(set([samples[sample]['model'] for sample in tumors]))
#print (models)

#### Rules
include: "rules/align.smk"
include: "rules/coverage.smk"
# include: "rules/freec-control.smk"
# include: "rules/merge-cnvs.smk"
include: "rules/mutect2.smk"
include: "rules/delly.smk"
# include: "rules/data_vis.smk"

# possible targets for rule all. 
		# ['%s/%s/aligned/filt.bam' % (in_data, s) for s in samples],
		# ['%s/%s/alfred_qc/qc_report.pdf' % (derived, s) for s in samples],
		# ['%s/%s/coverage/cov.gen.pdf' % (derived, s) for s in samples],
		# expand (['%s/%s/freec_{c}/freec_CNVs.p.value.txt' % (derived, s) for s in tumors], c=['control','single']),
		# expand (['%s/%s/freec_{c}/freec_CNVs.p.value.txt' % (derived, s) for s in controls], c=['single'])
		# ['%s/snv_calls/%s/mutect2-gnomAD-joint/somatic.vcf.gz' % (derived, m) for m in models]
		# ['%s/snv_calls/RPE_TP53-bypass/mutect2/somatic.vcf.gz' % (derived)],
		# ['%s/snv_calls/%s/mutect2-gnomAD-joint/somatic.vcf.gz' % (derived, m) for m in models],
		# expand(['%s/cnv_calls/%s/control-freec/cnvs_{c}.tsv' % (derived, m) for m in models], c=['control','single'])
		# ['%s/%s/delly/filtered_calls.bcf' % (derived, s) for s in tumors],

rule all:
	input:
		['%s/%s/aligned/raw.mdups.recal.bam' % (in_data, s) for s in samples],
		expand(['%s/%s/QC_plots/{qcplot}.pdf' % (derived, s) for s in samples], qcplot=['cov.gen','qc_report']),
		# ['%s/%s/delly/somatic_svs.vcf.gz' % (derived, s) for s in tumors],
		# ['%s/%s/mutect_calls/somatic.vcf.gz' % (derived, s) for s in tumors],

