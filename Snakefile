configfile: "config.yaml"

#### Helper variables

samples = config["samples"]
controls = ["WT_Ctrl"]
tumors = list (set (samples) - set (controls))
in_data = config["in_data"]
derived = config["derived_data"]
metadata = config["metadata"]
genome = config["ref_genome"]["name"]
logs = config["logs"]

models = list(set([samples[sample]['model'] for sample in tumors]))
#print (models)

#### Rules
#include: "rules/align.smk"
include: "rules/coverage.smk"
include: "rules/freec-control.smk"
include: "rules/mutect2.smk"
include: "rules/delly.smk"
include: "rules/data_vis.smk"

# possible targets for rule all. 
		# ['%s/%s/aligned/filt.bam' % (in_data, s) for s in samples],
		# ['%s/%s/alfred_qc/qc_report.pdf' % (derived, s) for s in samples],
		# ['%s/%s/coverage/cov.gen.pdf' % (derived, s) for s in samples],
		# expand (['%s/%s/freec_{m}/freec_CNVs.p.value.txt' % (derived, s) for s in tumors], m=['control','single']),
		# expand (['%s/%s/freec_{m}/freec_CNVs.p.value.txt' % (derived, s) for s in controls], m=['single'])
		# ['%s/snv_calls/%s/mutect2/somatic.vcf.gz' % (derived, c) for c in models]
		# ['%s/snv_calls/RPE_TP53-bypass/mutect2/somatic.vcf.gz' % (derived)],

rule all:
	input:
		# ['%s/%s/freec_BAF/freec_CNVs.p.value.txt' % (derived, s) for s in tumors],
		# ['%s/sv_calls/%s/delly/somatic_svs.bcf' % (derived, c) for c in models],
		expand (['%s/%s/freec_{m}/freec_CNVs.p.value.txt' % (derived, s) for s in tumors], m=['control','single']),
		['%s/snv_calls/%s/mutect2/somatic.vcf.gz' % (derived, c) for c in models],
		['%s/%s/delly/filtered_calls.bcf' % (derived, s) for s in tumors],
		#['%s/%s.filt.bam.tdf' % (config["igv_vis"], s) for s in samples],
		['%s/%s/mutect_calls/somatic.vcf.gz' % (derived, s) for s in tumors],

