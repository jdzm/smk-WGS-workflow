# looking at WGS processed with freec_control
set.seed(1225)
suppressPackageStartupMessages(library (GenomicRanges))
library (tidyverse); library (RColorBrewer)
home = T
basepath = ifelse (home, "/Users/diazmiyar/mnt/", "/mnt/data3/Juan/")
setwd (basepath)
## Runs on single samples
projectpath = paste0 (basepath, "20200103-WGS/")
blacklist = read.table(paste0(basepath, "repository/hg19/hg19-blacklist.v2.bed"), sep = '\t') %>% 
  separate (V1, into = c("V1", "Chromosome"), sep = "chr") %>% select (-V1) %>% 
  mutate (Chromosome = factor (Chromosome, levels = c(1:22, "X"))) %>% arrange(Chromosome) %>% as.tbl() %>%
  mutate (V2 = V2+1, V3 = V3+1)
control = "/freec_single/"
ploidy = 2
maxPloidyLevel = 6
ylimits = c(0, maxPloidyLevel)


all_samples = c('WT_Ctrl', 'TP53_Ctrl', 'TP53_48Rec1', 'TP53_48Rec2','TP53_Aneu1','TP53_Aneu2')
cnp.all = data.frame()
for (sam in all_samples){
  path_freec = paste0(projectpath, "derived_data/", sam, control)
  cnvs = read.table(paste0(path_freec, "freec_ratio.txt"), header = T) %>% as.tbl() %>% mutate (sample = sam)
  cnp.all = bind_rows(cnp.all, cnvs) %>% as.tbl()
}
chromsizes = read.table("~/mnt/repository/hg19/hg19.main.chrom.sizes", header = F, stringsAsFactors = F) %>% 
  mutate (Chromosome = row_number(), Chromosome = ifelse (V1=='chrX', 'X', ifelse(V1=='chrY','Y',Chromosome))) %>% 
  mutate (Chromosome = factor (Chromosome, levels = c(1:22, "X", 'Y'))) %>% arrange(Chromosome) %>% 
  mutate (newEnd = cumsum (as.numeric(V2)), newStart =cumsum(as.numeric(V2))- V2+1) %>% as.tbl()

newCoords = cnp.all %>% mutate (Chromosome = factor (Chromosome, levels = c(1:22, "X", 'Y'))) %>% arrange(Chromosome) %>%
  left_join(chromsizes, by = "Chromosome") %>% mutate (Start = Start + newStart) %>%
  mutate (PloidyCol = ifelse(CopyNumber==2, "2", ifelse (CopyNumber<2, "1", "3")))%>% 
  mutate (sample = factor (sample, levels = all_samples))

#median ratio for the whole fragment resulted from segmentation
newCoords %>% 
  filter (Ratio > -1) %>%
  mutate (Ratio = replace (Ratio, Ratio >= maxPloidyLevel/2, maxPloidyLevel/2), 
          CopyNumber = replace (CopyNumber, CopyNumber >= maxPloidyLevel, maxPloidyLevel)) %>%
  ggplot (aes (x = Start, y = Ratio*ploidy, color = PloidyCol))+
  geom_hline(yintercept = c(1,2,3,4,6,8), color = '#cacdcc', alpha = 0.6)+
  geom_vline(aes (xintercept = newStart), linetype = 'dashed', color = '#7b8280', alpha = 0.7)+
  geom_point(size = .8) + theme_bw(base_size = 16)+
  theme (legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  scale_color_manual(values = c('#428bca', '#5cb85c', '#d9534f'))+
  ylab ("normalized copy number profile") +xlab ("Chromosome") + coord_cartesian (ylim = ylimits)+
  facet_wrap(~sample, ncol = 1)

newCoords %>% 
  filter (MedianRatio > -1) %>%
  mutate (MedianRatio = replace (MedianRatio, MedianRatio >= maxPloidyLevel/2, maxPloidyLevel/2), 
          CopyNumber = replace (CopyNumber, CopyNumber >= maxPloidyLevel, maxPloidyLevel)) %>%
  ggplot (aes (x = Start, y = MedianRatio*ploidy, color = PloidyCol))+
  geom_hline(yintercept = c(1,2,3,4,6,8), color = '#cacdcc', alpha = 0.6)+
  geom_vline(aes (xintercept = newStart), linetype = 'dashed', color = '#7b8280', alpha = 0.7)+
  geom_point(size = .8) + theme_bw(base_size = 16)+
  theme (legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  scale_color_manual(values = c('#428bca', '#5cb85c', '#d9534f'))+
  ylab ("normalized copy number profile") +xlab ("Chromosome") + coord_cartesian (ylim = ylimits)+
  facet_wrap(~sample, ncol = 1)


cnp.all %>% group_by(MedianRatio, sample) %>% mutate (id = row_number())
cnp_ids = data.frame() 
for (sample in unique(cnp.all$sample)) {
  setsub = cnp.all %>% filter (sample == sample) %>% mutate (id=1)
  id=1
  for (rowid in 1:nrow(setsub)){
    if (rowid ==1){
      # skip
    } else {
      test = setsub[rowid,"MedianRatio"] == setsub[(rowid-1),"MedianRatio"] & setsub[rowid,"Chromosome"] == setsub[(rowid-1),"Chromosome"]
      if (!test) {id = id + 1}
      setsub[rowid,"id"] = id
    }
  }
  cnp_ids =bind_rows(cnp_ids, setsub)
}

binsize = 50000
setsub_reloaded = setsub %>% filter (!is.na(id)) %>%
  group_by (sample, id, sample, Chromosome) %>% 
  summarise (start = min (Start), end = (max(Start) + binsize -1), MedianRatio = unique(MedianRatio), CopyNumber = unique(CopyNumber), 
             CopyNumber = factor (CopyNumber, levels = c(1:maxPloidyLevel)))

setsub_reloaded %>% #filter (Chromosome==1) %>% 
  mutate (MedianRatio = replace (MedianRatio, MedianRatio < 0, 0), MedianRatio = replace (MedianRatio, MedianRatio*2>6, 6)) %>%
  ggplot (aes(color = CopyNumber))+
  geom_segment(aes (x = start, xend = end, y = MedianRatio*2, yend = MedianRatio*2), 
               size = 1) +theme_bw() +facet_wrap(~Chromosome, scales = "free_x")+ coord_cartesian(ylim = c(0,6))



# tranform this ploidy segments to a bedfile and mask blacklisted regions, maybe also low complexity regions

blacklist = read.table(paste0(basepath, "repository/hg19/hg19-blacklist.v2.bed"), sep = '\t') %>% 
  mutate (V1 = factor (V1, levels = paste0 ('chr', c(1:22, "X", 'Y')))) %>% arrange(V1) %>% as.tbl() %>%
  mutate (V2 = V2+1, V3 = V3+1, V1 = as.character(V1)) %>% dplyr::rename (chr=V1, start = V2, end=V3, type=V4) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T, ignore.strand = T)

cnvs = read.table(paste0(path_freec, "freec_ratio.txt"), header = T, stringsAsFactors = F) %>% as.tbl() %>% 
  mutate (sample = sam) %>% dplyr::rename (chr=Chromosome, start=Start, score = MedianRatio) %>% 
  mutate (chr = paste0('chr', chr), end = start+49999)%>% 
  select (-Ratio, -CopyNumber) %>% makeGRangesFromDataFrame(keep.extra.columns = T, ignore.strand = T)



# loop to collapse positions
subsetByOverlaps(cnvs, blacklist, invert = T) %>% as.data.frame() %>% 
  
  
findOverlaps (cnvs, blacklist)
countOverlaps(cnvs, blacklist) %>% sum()



cnvs


cnv.col = data.frame()
cnp.all %>% filter (sample == 'WT_Ctrl') %>% group_by ()

####
all_samples = c('WT_Ctrl', 'TP53_Ctrl', 'TP53_48Rec1','TP53_Aneu1', 
                'K562_Control', 'K562_Treated1')
cnp.all = data.frame()
for (sam in all_samples){
  path_freec = paste0(projectpath, "derived_data/", sam, "/freec_single/")
  cnvs = read.table(paste0(path_freec, "freec_ratio.txt"), header = T) %>% as.tbl() %>% mutate (sample = sam)
  cnp.all = bind_rows(cnp.all, cnvs) %>% as.tbl()
}

bin.order = c(1:22, "X", "Y")
cnp.all %>% filter (Ratio >0) %>% group_by(sample) %>% mutate (sample_Median = median (Ratio)) %>%
  group_by(sample, Chromosome) %>% mutate (chrom_Median = median (Ratio)) %>% 
  select (sample, Chromosome, sample_Median, chrom_Median) %>% distinct %>% ungroup %>%
  mutate (Chromosome = factor(Chromosome, levels = bin.order)) %>% arrange (sample, Chromosome) %>% 
  pivot_longer(cols = c(3,4), values_to = "value") %>%
  ggplot (aes (x = Chromosome, y = value*2, color = sample))+
  geom_hline(yintercept = c(2,3,4), color = '#cacdcc', alpha = 0.6)+
  geom_point () + facet_wrap(~name, ncol = 1)+
  theme_bw(base_size = 14)+ggsci::scale_color_futurama() +coord_cartesian(ylim = c(0,4))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  

sum_cnp = cnp.all %>% filter (Ratio >0) %>% group_by(sample) %>% summarise (sample_Median = median (Ratio), sample_Mean = mean (Ratio)) %>% as.data.frame
sum_cnp$ploidy = c(3,3,2,2,2,2)
sum_cnp %>% mutate (sample_Median = sample_Median*ploidy,sample_Mean = sample_Mean*ploidy )

#### Look at CNVs directly

cnvs_single
cnvs_control 

all.samples = c('WT_Ctrl', 'TP53_Ctrl', 'TP53_48Rec1','TP53_48Rec2','TP53_Aneu1', 'TP53_Aneu2', 
                'K562_Control', 'K562_Treated1', 'K562_Treated2')
cnvs.all = data.frame()
for (sam in all_samples){
  for (control in c("/freec_single/", "/freec_control/")){
    if (!((sam == "K562_Control" | sam == "WT_Ctrl") & (control == "/freec_control/"))){
      path_freec = paste0(projectpath, "derived_data/", sam, control)
      #message("loading ", path_freec)
      cnvs = read.table(paste0(path_freec, "freec_CNVs.p.value.txt"), header = T, sep = '\t', stringsAsFactors = F) %>% 
        as.tbl() %>% mutate (sample = sam, control = ifelse (grepl('control', control), T, F), chr = as.character(chr))
      cnvs.all = bind_rows(cnvs.all, cnvs) %>% as.tbl()
    }
  }
}
cnvs.all$sample = factor (cnvs.all$sample, levels = all.samples)

cnvs.all %>% filter (control == F)

cnvs.all %>% filter (control == F, !grepl("K562", sample)) %>% 
  ggplot ()+
  geom_bar (aes (x = sample, fill = WilcoxonRankSumTestPvalue< 0.05), position = 'dodge')+
  theme_classic()+ ggsci::scale_fill_jco(name = "p<0.05")+
  theme (axis.text.x = element_text (angle = 45, vjust = 0.6)) + facet_wrap(~status)


cnv_filt_control = cnvs.all %>% filter (control == T, !grepl("K562", sample), 
                     WilcoxonRankSumTestPvalue<0.05) 
bin.order = c(1:22, "X","Y")
custom.order = bin.order[which(bin.order %in% unique (cnv_filt_control$chr))]

cnv_filt_control %>% mutate (chr = factor (chr, levels = custom.order)) %>%
  ggplot ()+
  geom_bar (aes (x = sample, fill = status), position = position_dodge(preserve='single') )+
  theme_classic(base_size = 14)+ ggsci::scale_fill_jco(name = "CNA_status")+
  theme (axis.text.x = element_text (angle = 45, vjust = 0.6)) 

cnv_filt_control %>% mutate (chr = factor (chr, levels = rev(custom.order)), 
                             sample = factor (sample, levels = all.samples[c(3:6,2)])) %>%
  ggplot ()+
  geom_bar (aes (x = chr, fill = status), position = position_dodge(preserve='single') )+
  theme_classic(base_size = 14)+ ggsci::scale_fill_jco(name = "CNA_status")+
  #theme (axis.text.x = element_text (angle = 45, vjust = 0.6)) +
  facet_wrap(~sample,ncol = 2)+ coord_flip()
  

cnvs.bed = cnv_filt_control %>% select (-contains("Pvalue"), -control) %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
cnvs.TP = cnvs.bed[cnvs.bed$sample =="TP53_Ctrl"]

cnvs.altset = IRanges::subsetByOverlaps(cnvs.bed, cnvs.TP, invert = T) %>% 
  as.data.frame %>% as.tbl %>% select (-strand) %>% rename (chr = seqnames)

custom.order = bin.order[which(bin.order %in% unique (cnvs.altset$chr))]
cnvs.altset %>% mutate (chr = factor (chr, levels = rev(custom.order))) %>%
  group_by (chr, status, sample) %>% count() %>%
  ggplot ()+
  geom_col (aes (x = chr, y = factor (n), fill = status), position = position_dodge(preserve='single'))+
  theme_classic(base_size = 14)+ ggsci::scale_fill_jco(name = "CNA_status")+
  #theme (axis.text.x = element_text (angle = 45, vjust = 0.6)) +
  facet_wrap(~sample,ncol = 1)+ ylab ("count")

cnvs.bed[cnvs.bed$sample =="TP53_Ctrl"]%>% 
  as.data.frame %>% as.tbl %>% select (-strand) %>% rename (chr = seqnames, CN = copy.number) %>% select (-width, -control)


# We can call them differently if we assume tetraploidy... compare! 
ploidy.cnvs = data.frame()
for (sam in c('TP53_48Rec1','TP53_48Rec2','TP53_Aneu1', 'TP53_Aneu2')){
  for (control in c("/freec_single/", "/freec_control/")){
    path_freec = paste0(projectpath, "derived_data/", sam, control)
    cnvs = read.table(paste0(path_freec, "freec_CNVs.p.value.txt"), header = T, sep = '\t', stringsAsFactors = F) %>% 
      as.tbl() %>% mutate (ploidy = "2n",sample = sam, control = ifelse (grepl('control', control), T, F), chr = as.character(chr))
    ploidy.cnvs = bind_rows(ploidy.cnvs, cnvs) %>% as.tbl()
  }
  for (control in c("/freec_single_4n/", "/freec_control_4n/")){
    path_freec = paste0(projectpath, "derived_data/", sam, control)
    cnvs = read.table(paste0(path_freec, "freec_CNVs.p.value.txt"), header = T, sep = '\t', stringsAsFactors = F) %>% 
      as.tbl() %>% mutate (ploidy = "4n",sample = sam, control = ifelse (grepl('control', control), T, F), chr = as.character(chr))
    ploidy.cnvs = bind_rows(ploidy.cnvs, cnvs) %>% as.tbl()
  }
}
ploidy.cnvs$sample = factor (ploidy.cnvs$sample, levels = all.samples)
ploidy_filt = ploidy.cnvs %>% filter (control == T, WilcoxonRankSumTestPvalue<0.05) 

cnvs.ploidybed = ploidy_filt %>% select (-contains("Pvalue"), -control) %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
cnvs.TP = cnvs.bed[cnvs.bed$sample =="TP53_Ctrl"]

cnvs.altset = IRanges::subsetByOverlaps(cnvs.ploidybed, cnvs.TP, invert = T) %>% 
  as.data.frame %>% as.tbl %>% select (-strand) %>% rename (chr = seqnames)

custom.order = bin.order[which(bin.order %in% unique (cnvs.altset$chr))]
cnvs.altset %>% mutate (chr = factor (chr, levels = custom.order)) %>%
  group_by (chr, status, sample,ploidy) %>% count() %>%
  ggplot ()+
  geom_col (aes (x = chr, y = factor (n), fill = status), position = position_dodge(preserve='single'))+
  theme_classic(base_size = 14)+ ggsci::scale_fill_jco(name = "CNA_status")+
  #theme (axis.text.x = element_text (angle = 45, vjust = 0.6)) +
  facet_grid(ploidy~sample)+ ylab ("count")


######





