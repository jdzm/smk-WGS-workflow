# smk-WGS-workflow
From FASTA to SNV, SV and CNV calls


The usual structure to follow with this workflow is the following: 
```
.
├── derived_data
├── input
├── logs
├── metadata
└── workflow
```
The `input` folder will contain the symlinks of the fasta files and the alignment files. All the output of the analysis with the only exception of the bamfiles will be stored in `derived_data` in subfolders corresponding to the individual samples. In `metadata` is where you can symlink the reference genome or other files that are necessary for individual steps like the breakdown of your samples into tumor and normal and the template config for Control-FREEC.

In general, this pipeline follows the GATK4 workflow to preprocess the alignment files and to call Single Nucleotide Variants with Mutect2. It also implements vcf2maf for SNV annotation, Copy Number Variant calling with control-FREEC and Structural Variant calling with Delly.

## 01. Preprocessing

Samples are aligned individually to human_g1k_hs37d5 (human genome equivalent to hg19 with decoy sequences and uptated bases from the 1000 genomes project) genome with decoy sequences. Following GATK4 best practices, datasets are aligned with `bwa mem`, followed a `samtools fixmate` step. Next, duplicates are marked with `MarkDuplicatesSpark` (picard) and base scores recalibrated. This results in a ready-to-analyse raw bamfile containing all original reads (unfiltered): `input/sample_name/aligned/raw.mdups.recal.bam`. 

### 01-1. Quality Control and Coverage plots

A coverage histogram is extracted using `bedtools genomecov` and plotted genome-wide and per-chomosome. Coverage is also extracted in 10k bins by `alfred count_dna` and a QC report is produced using [Alfred](https://www.gear-genomics.com/docs/alfred/).  

All QC is stored in `derived_data/sample_name/plots`. The rule `gather_QC` will also produce a summry file gathering all `qc.tsv.gz` files into `derived_data/QC_summary/qc_info.tsv`. 

## 02. Single-Nucleotide Variant Calling (Mutect2)

Somatic variants are called using GATK4 Mutect2 comparing each sample to its control (wild-type sample in most cases). The process is done Chromosome-wise to parallelize and reduce time. For filtering, a subsequent step of `LearnReadOrientationModel` and `FilterMutectCalls` are applied. Then, all variants from each sample are merged and stored in `derived_data/sample_name/mutect_calls/somatic.vcf.gz`. 

Somatic variant calling can be done in all samples together using the rules in `rules/mutect2.smk`. Given that it took almost a week to run it once you increase the number of samples, we decided to call individually for each tumor-control pair and then merge samples followed by manual checks of variant quality. 

Variant annotation is subsequently performed using Ensembl Variant Effect Predictot [VEP](http://www.ensembl.org/info/docs/tools/vep/index.html) and Funcotator (which was not used later on). The process of installing VEP is not very straight-forward and is documented in `rules/annotate.smk`. 

After the annotation step, all maf files belonging to one model cell line are concatenated and stored in `derived_data/snv_calls/model/annotated/`. Quality control of the variants and filtering is not included in this pipeline.

## 03. Copy-Number Variant Calling (Control-FREEC)

Control-FREEC is run on different modes: per-sample (`freec_single`), and compared to a control(`freec_control`). It is also run with two different resolutions, 10kb and 50kb bins. In each run, the algorithm tests which ploidy gives a better fit out of ploidies 2, 3 and 4; choosing the best fit (specified in `derived_data/sample_name/freec_control_50k/freec_info.txt`. Sample purity is assumed to be 1 since we work with cell lines. 

The output is stored in the derived folder under each sample name for single and control runs. A per-sample copy number profile plot is generated from the `freec_ratio.txt` raw calls. A segmentation file for all CNVs with significance values can be found in `freec_CNVs.p.value.txt`.

The script `scripts/mergeFreec.R` combines all output per model and per call type and stores it in `derived_data/cnv_calls/model/control-freec/`. 

An alternative Snakefile and config file can be used to change the sample used as control (code in `alt.snakefile` and `config-alt.yaml`).

## 04. Structural Variant Calling (Delly)

As with SNVs, somatic SVs are called per sample-control pair or in a merged fashion. In order to multi-thread the process, it is usefult to use the precompiled version of Delly (check `config.yaml` for specific version). To parallelize, Delly is run 4 times per sample pair, one for each SV type (BND-translocations, DUP, DEL and INV) and then merged into `derived_data/sample_name/delly/merged_somatic.vcf.gz`. 

Due to repeated problems with the merging of the variants, it was done in a more 'rudimentary way' out of the smk workflow. The script `scripts/delly-cleanup.sh` contains info on how to merge variants and extract the coverage and support per sample to a txt file. Additionally, the script `scripts/delly-cleanup.R` contains code to plot variant support and to produce bedfiles for IGV visualization. Bed files can be found in `derived_data/sv_calls/model/delly_sv_annot_juicer`. 


#### Extra Files

Some rules in `rules/zz_archive` and scripts are no longer used and correspond to some tools that were tested but not implemented like `freebayes`, `facets` adnd `lumpy`. 

Conda environment specifications with program versions can be found in `envs`.

 
   
