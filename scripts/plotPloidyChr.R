#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# looking at WGS processed with freec_control
set.seed(1225)
suppressPackageStartupMessages(library (tidyverse))

# ## Runs on single samples
# sample = snakemake@params[["sam_id"]]
# cnvs = read.table(snakemake@input[["cnp"]], header = T) %>% as.tbl()
# #maxPloidyLevel = snakemake@params[["maxploidy"]]
# ploidy = snakemake@params[["ploidy"]] # expected ploidy
# ylimits = c(0, maxPloidyLevel)
# control = snakemake@params[["withcontrol"]]

## Runs on single samples
sample = args[1]
ratios = read.table(args[2], header = T) %>% as.tbl()
ploidy = args[3] # expected ploidy
if (ploidy == '2,3,4'){ploidy = 2}
outdir = args[4]
control = ifelse (grepl ("single", outdir), "single", "control") # takes also BAF as control
#maxPloidyLevel = snakemake@params[["maxploidy"]]
maxPloidyLevel = 6
ylimits = c(0, maxPloidyLevel)


message ("Running with params: ", sample, " ploidy of ",ploidy, ", maxPloidy of ",maxPloidyLevel)

chromsizes = ratios %>% 
  group_by (Chromosome) %>% mutate (maxChr = max (Start)) %>% 
  select (Chromosome, maxChr) %>% distinct() %>% ungroup() %>% 
  mutate (Chromosome = factor (Chromosome, levels = c(1:22, "X"))) %>% arrange(Chromosome) %>% 
  mutate (newEnd = cumsum (as.numeric(maxChr)), newStart =cumsum(as.numeric(maxChr))- maxChr+1)
#for (chr in 1:23){chromsizes[chr,"newStart"] = ifelse (chr == 1, 1, sum (chromsizes [1:(chr-1),"maxChr"]))}

newCoords = ratios %>% mutate (Chromosome = factor (Chromosome, levels = c(1:22, "X"))) %>% arrange(Chromosome) %>%
  left_join(chromsizes, by = "Chromosome") %>% mutate (Start = Start + newStart) %>%
  mutate (PloidyCol = ifelse(CopyNumber==ploidy, "2", ifelse (CopyNumber<ploidy, "1", "3")))

# plt_ploidy = newCoords %>% 
#   filter (Ratio > -1) %>% mutate (CopyNumber = replace (CopyNumber, CopyNumber >= maxPloidyLevel, maxPloidyLevel)) %>%
#   ggplot (aes (x = Start, y = CopyNumber, color = PloidyCol))+
#   geom_vline(aes (xintercept = newStart), linetype = 'dashed', color = '#7b8280')+
#   geom_point(size = 1)+theme_bw(base_size = 16)+
#   theme (legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())+
#   scale_color_manual(values = c('#428bca', '#5cb85c', '#d9534f'))+ 
#   ggtitle (ifelse (control =="single", paste(sample), paste (sample, 'with', control)))+
#   ylab ("predicted ploidy") +xlab ("Chromosome") +coord_cartesian (ylim = ylimits)

plt_ratio = newCoords %>% 
  filter (Ratio > -1) %>%
  mutate (Ratio = replace (Ratio, Ratio >= maxPloidyLevel/2, maxPloidyLevel/2)) %>%
  ggplot (aes (x = Start, y = Ratio*ploidy, color = PloidyCol))+
  geom_hline(yintercept = c(1,2,3,4,6), color = '#cacdcc', alpha = 0.6)+
  geom_vline(aes (xintercept = newStart), linetype = 'dashed', color = '#7b8280', alpha = 0.7)+
  geom_point(size = .8) + theme_bw(base_size = 16)+
  theme (legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  scale_color_manual(values = c('#428bca', '#5cb85c', '#d9534f'))+
  ggtitle (ifelse (control =="single", paste(sample), paste (sample, 'with', control)))+
  ylab ("normalized copy number profile") +xlab ("Chromosome") + coord_cartesian (ylim = ylimits)

# png(filename = snakemake@output[["plo"]], width = 1000, height = 400)
# plot (plt_ploidy)
# dev.off()

png(filename = paste0(outdir, "freec_ratio.png"), width = 1000, height = 400)
plot (plt_ratio)
dev.off()


plt_chrom = ratios %>% mutate (Chromosome = factor (Chromosome, levels = c(1:22, "X"))) %>% arrange(Chromosome) %>%
  mutate (PloidyCol = ifelse(CopyNumber==ploidy, "2", ifelse (CopyNumber<ploidy, "1", "3"))) %>%
  filter (Ratio > -1) %>%
  mutate (Ratio = replace (Ratio, Ratio >= maxPloidyLevel/2, maxPloidyLevel/2)) %>%
  ggplot (aes (x = Start/10^6, y = Ratio*ploidy, color = PloidyCol))+
  geom_hline(yintercept = c(1,2,3,4,6), color = '#cacdcc', alpha = 0.6)+
  #geom_vline(aes (xintercept = newStart), linetype = 'dashed', color = '#7b8280', alpha = 0.7)+
  geom_point(size = .8) + theme_bw(base_size = 12)+
  theme (legend.position = "none")+ facet_wrap(~Chromosome, scales="free_x", ncol = 5)+
  scale_color_manual(values = c('#428bca', '#5cb85c', '#d9534f'))+
  ggtitle (ifelse (control =="single", paste(sample), paste (sample, 'with', control)))+
  ylab ("normalized copy number profile") +xlab ("Chromosome position (Mb)")+
  coord_cartesian (ylim = ylimits)

png(filename = paste0(outdir, "freec_bychrom.png"), width = 1000, height = 1000)
plot (plt_chrom)
dev.off()

