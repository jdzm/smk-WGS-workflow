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


#### Define target files

rule all:
	input:
		# ['%s/%s/aligned/raw.mdups.recal.bam' % (in_data, s) for s in samples],
		expand(['%s/%s/QC_plots/{qcplot}.pdf' % (derived, s) for s in samples], qcplot=['cov.gen','qc_report', 'cov.gz']),
		# ['%s/QC_summary/qc_info.tsv' % (derived)],
		# ['%s/%s/delly/somatic_svs.vcf.gz' % (derived, s) for s in tumors],
		expand(['%s/%s/delly/somatic_{sv}.vcf.gz' % (derived, s) for s in tumors], sv=['BND','DEL']),
		['%s/%s/mutect_calls/somatic.vcf.gz' % (derived, s) for s in tumors],
		# ['%s/sv_calls/RPE_TP53/delly/somatic_BND.vcf.gz' % (derived)],

