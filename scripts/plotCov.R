#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# plots histogram of genome coverage
set.seed(1225)
suppressPackageStartupMessages(library(tidyverse))

sample = args[1]
#sample = snakemake@params[["sam_id"]]
#read_tsv(snakemake@input[["hist"]]

covhist = read_tsv(args[2], col_names = c('chr', 'cov','bases','total','ratio'), 
	col_types = cols(chr = col_character(),cov = col_double(),
		bases = col_double(),total = col_double(),ratio = col_double())) %>% 
  filter (chr != "chrM") %>% mutate (chr = factor (chr, levels = c(paste0('chr',c(1:22, 'X','Y')), 'genome')))

gencov = covhist %>% filter (chr == 'genome') 

histplot = gencov %>%
  ggplot (aes (x = cov, y = ratio)) +
  geom_col (fill = 'steelblue') + ylab ('genome fraction') + xlab ('coverage')+
  theme_bw(base_size=14)+ggtitle (paste (sample, "Coverage"))

chrplot = covhist %>% 
  ggplot (aes (x = cov, y = ratio)) +
  geom_col (fill = 'steelblue') + ylab ('genome fraction') + xlab ('coverage')+
  theme_bw(base_size=14)+ggtitle (paste (sample, "Coverage"))+
  facet_wrap(~chr)

outdir = dirname (args[2])
ggsave (filename = paste0(outdir, "/cov.gen.pdf"), 
  plot = histplot, device = "pdf", width=6, height=5)
ggsave (filename = paste0(outdir, "/cov.chr.pdf"), 
  plot = chrplot, device = "pdf", width=7, height=7)
# ggsave (filename = snakemake@output[['plot_cov1']], 
# 	plot = histplot, device = "pdf", width=6, height=5)
# ggsave (filename = snakemake@output[['plot_cov2']], 
# 	plot = chrplot, device = "pdf", width=7, height=7)