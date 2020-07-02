# qc_gather.R
# input qc_reps
# outdirs params.outdir
set.seed(1225)
suppressPackageStartupMessages(library (tidyverse))

qc_files = unlist(snakemake@input['qc_reps'])
# need to add bit to extract 

qc_project=data.frame()
for (filename in qc_files){
  if (file.exists(filename)){
    if (nrow(qc_project) ==0 ){
      cmd = paste("zgrep ^ME", filename,"| cut -f 2- | datamash transpose | awk '{print $1,$2}'") 
      qc_project = bind_rows (qc_project, read.table(pipe (cmd), stringsAsFactors = F, header = F, comment.char = "") )
    } else {
      cmd = paste("zgrep ^ME", filename,"| cut -f 2- | datamash transpose | awk '{print $1,$2}'") 
      qc_project = left_join (qc_project, read.table(pipe (cmd), stringsAsFactors = F, header = F, comment.char = "") , by = 'V1')
    }
  }
}

qc_mat = t(as.matrix (qc_project %>%rowwise() %>% mutate (V1 = replace(V1, grepl("#",V1), gsub("#", "N", V1)))))
colnames (qc_mat) = qc_mat[1,]
to_plot = as.data.frame(qc_mat[-1,], stringsAsFactors = F) %>% as_tibble()

to_export = to_plot %>% pivot_longer(cols = 2:ncol(to_plot), names_to = "metric", values_to = "value") %>% pivot_wider(names_from = Sample)
# outdir = as.character(snakemake@params["outdir"])
outdir = dirname (as.character(snakemake@output["info_table"]))
dir.create(outdir, showWarnings = F)
# write.table(to_export, paste0(outdir, "qc_info.tsv" ), quote = F, sep = "\t", row.names = F)
write.table(to_export, as.character(snakemake@output["info_table"]), quote = F, sep = "\t", row.names = F)
