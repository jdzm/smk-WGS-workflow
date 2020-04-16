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


