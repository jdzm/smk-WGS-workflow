# merge 
set.seed(1225)
suppressPackageStartupMessages(library (GenomicRanges)) # some of the functions are masked by dplyr
suppressPackageStartupMessages(library (tidyverse))
suppressPackageStartupMessages(library (readr))

cnv_files = unlist(snakemake@input['cnvs'])
cnv.all = tibble()
newcolnames = c("chr","start","end","copyNumber","status","Wilcoxpval","KSpval", "sample")
for (filename in cnv_files) {
  # load file
  sam_name = unlist(strsplit(filename, split = "/"))[3]
  cnvs = read_tsv(filename, col_types = cols(.default = "c")) %>% mutate (sample = sam_name)
  names (cnvs) = newcolnames
  cnv.all = bind_rows(cnv.all, cnvs)# %>% as.tbl()
}
write.table(cnv.all, file = as.character(snakemake@output["cnv_sam"]), quote = F, row.names = F, sep = '\t')

ratio_files = unlist(snakemake@input['ratios'])
ratio.all = tibble()
basenames = c("chr","start","ratio","medRatio","copyNumber","subclone_CN","subclone_popul")
for (filename in ratio_files) {
  # load file
  sam_name = unlist(strsplit(filename, split = "/"))[3]
  ratios = read_tsv(filename, col_types = cols(.default = "c"), progress = FALSE) %>% mutate (Start = as.integer(Start)-1)
  
  names (ratios) = c(basenames[1:2], paste0(basenames[3:7], '.',sam_name))
  if (which (filename == ratio_files) == 1){
    ratio.all = ratios
  } else {
    ratio.all = ratio.all %>% left_join(ratios, by = basenames[1:2])   
  }
}

blacklist_file = as.character(snakemake@params["blacklist"])
blacklist = read.table(blacklist_file, sep = '\t', stringsAsFactors = F) %>% separate (V1, into = c("V1", "chr"), sep = "chr") %>% select (-V1) %>% 
  mutate (chr = factor (chr, levels = c(1:22, "X", "Y"))) %>% arrange(chr) %>% as.tbl() 
dellyExcl_file = as.character(snakemake@params["excl"])
dellyExcl = read.table(dellyExcl_file, sep = '\t',nrows = 140, stringsAsFactors = F) %>% filter(!grepl("chr", V1)) %>% rename (chr=V1) %>%
  mutate (chr = factor (chr, levels = c(1:22, "X", "Y"))) %>% arrange(chr) %>% as.tbl() 

blacklist.gr = GRanges(blacklist$chr,IRanges (blacklist$V2,blacklist$V3), "Type" = blacklist$V4)
dellyExcl.gr = GRanges(dellyExcl$chr,IRanges (dellyExcl$V2,dellyExcl$V3), "Type" = dellyExcl$V4)
rat.gr = GRanges(ratio.all$chr,IRanges (ratio.all$start,ratio.all$start+49000))
o_bl = findOverlaps(rat.gr, blacklist.gr)
o_delly = findOverlaps(rat.gr, dellyExcl.gr)

dellyTypes = o_delly %>% as.data.frame() %>% as.tbl() %>% rowwise() %>%
  mutate(dellyType = as.character(dellyExcl[subjectHits,"V4"])) %>%
  mutate(chr = as.character(ratio.all[queryHits,1]), start = as.integer(ratio.all[queryHits,2])) %>% select (-queryHits, -subjectHits)
blTypes = o_bl %>% as.data.frame() %>% as.tbl() %>% rowwise() %>%
  mutate(blType = as.character(blacklist[subjectHits,"V4"])) %>%
  mutate(chr = as.character(ratio.all[queryHits,1]), start = as.integer(ratio.all[queryHits,2])) %>% select (-queryHits, -subjectHits)

ratio.all_bl = ratio.all %>% left_join(blTypes, by = c("chr", "start")) %>% left_join(dellyTypes, by = c("chr", "start")) 
write.table(ratio.all_bl, file = as.character(snakemake@output["ratio_mat"]), quote = F, row.names = F, sep = '\t')


info_files = unlist(snakemake@input['info'])
info.all = tibble()
for (filename in info_files) {
  # load file
  sam_name = unlist(strsplit(filename, split = "/"))[3]
  infos = read.table(filename, header = F, sep = '\t', stringsAsFactors = F) %>% 
    mutate (sample = sam_name, V1 = replace (V1, V1=="Sample_Name", "Bam_File"))

  info.all = info.all %>% bind_rows(infos)   
}

write.table(info.all %>% pivot_wider(names_from = sample, values_from = V2), 
            file = as.character(snakemake@output["freec_info"]), quote = F, row.names = F, sep = '\t')



