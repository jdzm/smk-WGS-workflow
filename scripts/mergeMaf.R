#!/usr/bin/env Rscript
## Concatenate MAF files and filter them with filter == PASS
suppressPackageStartupMessages(require (readr))
suppressPackageStartupMessages(require (tidyverse))

maf_files=unlist(snakemake@input['maf'])

maf.all = tibble()
for (filename in maf_files) {
  # load file
  maf_df = read_tsv(filename, col_types = cols(.default = "c"), comment = "#") 
  maf.all = bind_rows(maf.all, maf_df)# %>% as.tbl()
}
write.table(maf.all %>% filter (FILTER == 'PASS'), 
	file = as.character(snakemake@output["joint_maf"]), quote = F, row.names = F, sep = '\t')
write.table(maf.all, 
	file = as.character(snakemake@output["joint_raw"]), quote = F, row.names = F, sep = '\t')