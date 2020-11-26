my_packages = c ("tidyverse", "maftools", "vcfR",'RColorBrewer', 'ggsci', 'patchwork', 'ggridges', 'scales')
msg = paste0(suppressPackageStartupMessages(lapply (my_packages, require, character.only = T)))

all.samples = c('TP53_48h','TP53_6w','TP53_6w1', 'TP53_6w2', 'TP53_6w3',
                 'TP53_6wTumo11','TP53_6wTumo12','TP53_6wTumo13',
                 'TP53_6wTumo21','TP53_6wTumo22','TP53_6wTumo23',
                 'TP53_6wTumo31','TP53_6wTumo32','TP53_6wTumo33')

# sample_df = cbind.data.frame(samples = c(all.samples1, all.samples2), batch = c(rep('wgd1', length(all.samples1)), rep ('wgd2', length(all.samples2)))) %>% 
#   mutate (path = ifelse (batch == 'wgd1', 'mnt/20200424-Variants/derived_data/', 'mnt/20200610-WGD-E02/derived_data/'))
# 
# sample_info = read_tsv('Desktop/sample-rename.txt', col_types = cols())
# all.samples = sample_info$newname

outdir = 'Desktop/WGS_PlotGallery/20200730-dellysv/'
dir.create(outdir)

chrs = c(1:22,'X')
mypal = ggsci::pal_nejm()(4)
names (mypal) = c('BND','DUP','DEL','INV')
mypalRGB = col2rgb(mypal) %>% t() %>% as.data.frame() %>% rownames_to_column(var = 'svtype') %>% 
  mutate (color = paste (red, green, blue, sep = ',')) %>% select (svtype, color)


#####
dellycolnames = c('chr', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'CHR2', 'END','PRECISE', 'IMPRECISE', 'SVLEN',
                  'INSLEN', 'HOMLEN', 'PE', 'SR', 'RV', 'RR', 'DR', 'DV', 'CN', 'GT', 'GQ', 'GL', 'FT')
sv_list = list()
for (sam1 in all.samples){
  sv_sample = read_delim(paste0('mnt/20200610-WGD-E02/derived_data/',sam1,'/delly/merged_svs/merged_svs_',sam1,'.txt'), 
                         col_names = dellycolnames, delim = ' ',
                         col_types = cols(chr = col_character(), CHR2 = col_character())) %>%
    mutate (RV = as.integer(RV))
  
  sv_control = read_delim(paste0('mnt/20200610-WGD-E02/derived_data/',sam1,'/delly/merged_svs/merged_svs_WT_Ctrl.txt'), 
                          col_names = c(dellycolnames[1:15], paste0(dellycolnames[16:24], '_N')), delim = ' ', 
                          col_types = cols(chr = col_character(), CHR2 = col_character()))%>%
    mutate (RV_N = as.integer(RV_N))
  
  # precise, PE and SR are common to both 
  bnd = sv_sample %>% filter (grepl ('BND', ID), PRECISE ==1, GQ > 0, chr != 'Y', CHR2 != 'Y') %>% select (chr, POS, END, ID, ALT, CHR2, RR, RV, CN:GQ, FT) %>%
    left_join(sv_control %>% select (chr, POS, ID, ALT, CHR2, RR_N, RV_N, CN_N:GQ_N, FT_N), by = c("chr", "POS", "ID", "ALT", "CHR2")) %>%
    filter (RV_N <= 2, GQ > 10) %>% mutate (svtype = substr(ID, 1, 3)) %>% left_join(mypalRGB, by ='svtype')
  
  other_sv = sv_sample %>% filter (!grepl ('BND', ID), PRECISE ==1, GQ > 0, chr != 'Y', CHR2 != 'Y') %>% select (chr, POS, END, ID, ALT, CHR2, RR, RV, CN:GQ, FT) %>%
    left_join(sv_control %>% select (chr, POS, ID, ALT, CHR2, RR_N, RV_N, CN_N:GQ_N, FT_N), by = c("chr", "POS", "ID", "ALT", "CHR2")) %>%
    filter (RV_N == 0, GQ > 10) %>% mutate (svtype = substr(ID, 1, 3)) %>% left_join(mypalRGB, by ='svtype')
  
  filtered_svs = bind_rows(bnd, other_sv) %>% mutate_at (vars (POS, END), as.integer)

  bedpe_svs = bnd %>% select (chr, POS, END, CHR2, ID, ALT, color) %>%
    rename (chr1 = chr, x1 = POS, x2 = END, chr2=CHR2) %>%
    separate(ALT, into = c('a','b','y1','d'), remove = F) %>% select (-a,-b,-d) %>%
    mutate (y1 = as.integer (y1), y2=y1+1) %>%
    select(chr1:chr2, y1,y2, color, ID, ALT) %>%
    bind_rows(other_sv %>% select (chr, POS, END, ID, ALT, color) %>%
                rename (chr1 = chr, x1 = POS, x2 = END) %>%
                mutate (chr2= chr1, y1 = x1, y2 = x2) %>%
                select(chr1:x2, chr2 ,y1, y2, color, ID, ALT) ) %>%
    mutate (score='.',strand1='.',strand2='.') %>%
    select (chr1:y2, ID,score:strand2, ALT, color) %>% mutate_at (vars (x1, x2, y1, y2), as.integer)

  bed_svs = bnd %>% bind_rows(other_sv) %>%
    select (chr:ALT, GQ, svtype, color) %>%
    rename (`#chrom` = chr, start = POS, end = END, name = ID, score = GQ, itemRgb = color) %>%
    mutate (strand='.', thickStart='.', thickEnd='.', blockCount='.', blockSizes='.', blockStarts='.') %>%
    select (`#chrom`, start, end, name, score,strand ,thickStart ,thickEnd, itemRgb, blockCount, blockSizes, blockStarts) %>%
    mutate_at (vars (start, end), as.integer)

  # chrom start end name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts
  filename = paste0('mnt/20200610-WGD-E02/derived_data/sv_calls/RPE_TP53/delly_sv_annot_juicer/', gsub('TP53_', '',sam1), '_precise_svs')
  cat(paste0('track name="',sam1,'" description="precise SVs" itemRgb="On" \n'),file=paste0(filename, ".bed"))
  write.table(bed_svs, paste0(filename, ".bed"), sep = '\t', row.names = F,quote = F, append = TRUE,col.names = F)
  cat(paste0('track name="',sam1,'" description="precise SVs" itemRgb="On" \n'),file=paste0(filename, ".bedpe"))
  write.table(bedpe_svs, paste0(filename, ".bedpe"), sep = '\t', row.names = F,quote = F, append = TRUE,col.names = F)
  
  write_tsv(filtered_svs %>% mutate (Sample = sam1), paste0('mnt/20200610-WGD-E02/derived_data/',sam1,'/delly/filt_precise_svs.tsv') ,col_names = T)
}

sv_list = list()
for (sam1 in all.samples){
  filtered_svs = read_tsv (paste0('mnt/20200610-WGD-E02/derived_data/',sam1,'/delly/filt_precise_svs.tsv'), col_types = cols(chr = 'c'))
  sv_list[sam1] = list (filtered_svs)
}
 
## copy and rename to projects/WGD-WholeGenome/svs/delly-calls

all_delly = tibble()
for (i in seq_along(sv_list)){
  all_delly = bind_rows(all_delly, sv_list[[i]])
}

sample_info = read_tsv('Desktop/sample-rename.txt', col_types = cols())

to_export = all_delly %>% left_join(sample_info %>% rename (Sample = oldname) %>% select (Sample, newname), by = 'Sample')
write_tsv(to_export, "mnt/projects/WGD-WholeGenome/svs/delly-calls/delly_calls_precise_filtered_wgd2.tsv")

####### IMPRECISE
dellycolnames = c('chr', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'CHR2', 'END','PRECISE', 'IMPRECISE', 'SVLEN',
                  'INSLEN', 'HOMLEN', 'PE', 'SR', 'RV', 'RR', 'DR', 'DV', 'CN', 'GT', 'GQ', 'GL', 'FT')
sv_list = list()
for (sam1 in all.samples){
  sv_sample = read_delim(paste0('mnt/20200610-WGD-E02/derived_data/',sam1,'/delly/merged_svs/merged_svs_',sam1,'.txt'), 
                         col_names = dellycolnames, delim = ' ',
                         col_types = cols(chr = col_character(), CHR2 = col_character())) %>%
    mutate (RV = as.integer(RV))
  
  sv_control = read_delim(paste0('mnt/20200610-WGD-E02/derived_data/',sam1,'/delly/merged_svs/merged_svs_WT_Ctrl.txt'), 
                          col_names = c(dellycolnames[1:15], paste0(dellycolnames[16:24], '_N')), delim = ' ', 
                          col_types = cols(chr = col_character(), CHR2 = col_character()))%>%
    mutate (RV_N = as.integer(RV_N))
  
  # imprecise, PE and SR are common to both 
  bnd = sv_sample %>% filter (grepl ('BND', ID), IMPRECISE ==1, GQ > 0, chr != 'Y', CHR2 != 'Y', FT == 'PASS') %>% select (chr, POS, END, ID, ALT, CHR2, DR, DV, CN:GQ, FT) %>%
    left_join(sv_control %>% select (chr, POS, ID, ALT, CHR2, DR_N, DV_N, CN_N:GQ_N, FT_N), by = c("chr", "POS", "ID", "ALT", "CHR2")) %>%
    filter (DV_N <= 2, GQ > 10) %>% mutate (svtype = substr(ID, 1, 3)) %>% left_join(mypalRGB, by ='svtype')
  
  other_sv = sv_sample %>% filter (!grepl ('BND', ID), PRECISE ==1, GQ > 0, chr != 'Y', CHR2 != 'Y') %>% select (chr, POS, END, ID, ALT, CHR2, DR, DV, CN:GQ, FT) %>%
    left_join(sv_control %>% select (chr, POS, ID, ALT, CHR2, DR_N, DV_N, CN_N:GQ_N, FT_N), by = c("chr", "POS", "ID", "ALT", "CHR2")) %>%
    filter (DV_N == 0, GQ > 10) %>% mutate (svtype = substr(ID, 1, 3)) %>% left_join(mypalRGB, by ='svtype')
  
  filtered_svs = bind_rows(bnd, other_sv) %>% mutate_at (vars (POS, END), as.integer)
  # 
  bedpe_svs = bnd %>% select (chr, POS, END, CHR2, ID, ALT, color) %>%
    rename (chr1 = chr, x1 = POS, x2 = END, chr2=CHR2) %>%
    separate(ALT, into = c('a','b','y1','d'), remove = F) %>% select (-a,-b,-d) %>%
    mutate (y1 = as.integer (y1), y2=y1+1) %>%
    select(chr1:chr2, y1,y2, color, ID, ALT) %>%
    bind_rows(other_sv %>% select (chr, POS, END, ID, ALT, color) %>%
                rename (chr1 = chr, x1 = POS, x2 = END) %>%
                mutate (chr2= chr1, y1 = x1, y2 = x2) %>%
                select(chr1:x2, chr2 ,y1, y2, color, ID, ALT) ) %>%
    mutate (score='.',strand1='.',strand2='.') %>%
    select (chr1:y2, ID,score:strand2, ALT, color) %>% mutate_at (vars (x1, x2, y1, y2), as.integer)
  
  bed_svs = bnd %>% bind_rows(other_sv) %>%
    select (chr:ALT, GQ, svtype, color) %>%
    rename (`#chrom` = chr, start = POS, end = END, name = ID, score = GQ, itemRgb = color) %>%
    mutate (strand='.', thickStart='.', thickEnd='.', blockCount='.', blockSizes='.', blockStarts='.') %>%
    select (`#chrom`, start, end, name, score,strand ,thickStart ,thickEnd, itemRgb, blockCount, blockSizes, blockStarts) %>%
    mutate_at (vars (start, end), as.integer)
  
  # chrom start end name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts
  filename = paste0('mnt/20200610-WGD-E02/derived_data/sv_calls/RPE_TP53/delly_sv_annot_juicer/', gsub('TP53_', '',sam1), '_imprecise_svs')
  cat(paste0('track name="',sam1,'" description="precise SVs" itemRgb="On" \n'),file=paste0(filename, ".bed"))
  write.table(bed_svs, paste0(filename, ".bed"), sep = '\t', row.names = F,quote = F, append = TRUE,col.names = F)
  cat(paste0('track name="',sam1,'" description="precise SVs" itemRgb="On" \n'),file=paste0(filename, ".bedpe"))
  write.table(bedpe_svs, paste0(filename, ".bedpe"), sep = '\t', row.names = F,quote = F, append = TRUE,col.names = F)
  
  write_tsv(filtered_svs %>% mutate (Sample = sam1), paste0('mnt/20200610-WGD-E02/derived_data/',sam1,'/delly/filt_imprecise_svs.tsv') ,col_names = T)
}

sv_list = list()
for (sam1 in all.samples){
  filtered_svs = read_tsv (paste0('mnt/20200610-WGD-E02/derived_data/',sam1,'/delly/filt_imprecise_svs.tsv'), col_types = cols(chr = 'c'))
  sv_list[sam1] = list (filtered_svs)
}

## copy and rename to projects/WGD-WholeGenome/svs/delly-calls

all_delly = tibble()
for (i in seq_along(sv_list)){
  all_delly = bind_rows(all_delly, sv_list[[i]])
}

sample_info = read_tsv('Desktop/sample-rename.txt', col_types = cols())

to_export = all_delly %>% left_join(sample_info %>% rename (Sample = oldname) %>% select (Sample, newname), by = 'Sample')
write_tsv(to_export, "mnt/projects/WGD-WholeGenome/svs/delly-calls/delly_calls_imprecise_filtered_wgd2.tsv")