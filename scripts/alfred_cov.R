library(ggplot2)
library(scales)

#chrs = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX")
chrs = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X")

args = commandArgs(trailingOnly=TRUE)
if (length(args) > 1) { pdffile = args[2]; } else { pdffile = paste0(args[1], ".pdf"); }
x = read.table(args[1], header=T)
x = x[x$chr %in% chrs,]
x$chr = factor(x$chr, levels=chrs)

# Iterate samples
for (i in 5:ncol(x)) {
    sample = colnames(x)[i]
    # print(sample)
    df = data.frame(chr=x$chr, start=x$start + (x$end - x$start) / 2, rd=log(x[,i] / median(x[,i]))/log(2))
    p1 = ggplot(data=df, aes(x=start, y=rd))
    p1 = p1 + geom_point(pch=21, size=0.5)
    p1 = p1 + xlab("Chromosome")
    p1 = p1 + ylab("Log2 median normalized read depth")
    p1 = p1 + scale_x_continuous(labels=comma)
    p1 = p1 + facet_grid(. ~ chr, scales="free_x", space="free_x")
    p1 = p1 + theme_bw()+theme(axis.text.x = element_text(angle=45, hjust=1))
    ggsave(pdffile, width=28, height=6)
    print(warnings())
}
