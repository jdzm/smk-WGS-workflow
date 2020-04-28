#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# plots histogram of genome coverage
set.seed(1225)
suppressPackageStartupMessages(library(tidyverse))

sample = args[1]
#sample = snakemake@params[["sam_id"]]
#read_tsv(snakemake@input[["hist"]]
chromsizes = read.table (args[2], stringsAsFactors = F) %>% mutate (V1 = factor (V1, levels = V1))

covhist = read_tsv(args[3], col_names = c('chr', 'cov','bases','total','ratio'),
  col_types = cols(chr = col_character(),cov = col_double(),
    bases = col_double(),total = col_double(),ratio = col_double())) %>%
  mutate (chr = factor (chr, levels = c(levels (chromsizes$V1), 'genome')))

contig_other = c("MT", as.character(chromsizes$V1[which (str_length(chromsizes$V1) > 3)]))

gencov = covhist %>% filter (chr == 'genome', cov != 0) 

histplot = gencov %>% filter (cov != 0) %>% 
  ggplot (aes (x = cov, y = ratio)) +
  geom_col (fill = 'steelblue') + ylab ('genome fraction') + xlab ('coverage')+
  theme_bw(base_size=14)+ggtitle (paste (sample, "Coverage"))

chrplot = covhist %>% filter (cov != 0, !(chr %in% contig_other)) %>%
  ggplot (aes (x = cov, y = ratio)) +
  geom_col (fill = 'steelblue') + ylab ('genome fraction') + xlab ('coverage')+
  theme_bw(base_size=14)+ggtitle (paste (sample, "Coverage"))+
  facet_wrap(~chr)

contigplot = covhist %>% filter (cov != 0, chr %in% contig_other) %>%
  ggplot (aes (x = cov, y = ratio)) +
  geom_col (fill = 'steelblue') + ylab ('genome fraction') + xlab ('coverage')+
  theme_bw(base_size=14)+ggtitle (paste (sample, "Coverage"))+
  facet_wrap(~chr)

outdir = dirname (args[3])
ggsave (filename = paste0(outdir, "/cov.gen.pdf"), 
  plot = histplot, device = "pdf", width=6, height=5)
ggsave (filename = paste0(outdir, "/cov.chr.pdf"), 
  plot = chrplot, device = "pdf", width=7, height=7)
ggsave (filename = paste0(outdir, "/cov.contigs.pdf"), 
        plot = contigplot, device = "pdf", width=7, height=7)
# ggsave (filename = snakemake@output[['plot_cov1']], 
#   plot = histplot, device = "pdf", width=6, height=5)
# ggsave (filename = snakemake@output[['plot_cov2']], 
#   plot = chrplot, device = "pdf", width=7, height=7)