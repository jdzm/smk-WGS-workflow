
library (facets)

set.seed(1223)
#snakemake@input[["pileup"]]
datafile = 'TP53_6wTumo1_mskcc_pileup_dbsnp.txt.gz'
rcmat = readSnpMatrix(datafile)
rcmat[1:10,]
xx = preProcSample(rcmat)

# A bivariate genome segmentation is performed on logR and logOR by extending the 
# CBS algotithm (Olshen et al., 2004; Venkatraman and Olshen, 2007) to the 
# bivariate scenario using a T2 statistic for identifiying change points.

# First some preprocessing to correct by several params

# Read depth ratio between tumor and normal gives information on total copy number

oo=procSample(xx,cval=150) # lower critical value cval leads to higher sensitivity for small changes

# logR changes are proportional to diploid. logR=0 is 2n. 
oo$dipLogR 

# Call allele-specific copy number and associated cellular fraction, estimate tumor purity and ploidy.
fit=emcncf(oo)

# once tha logR of diploid state is obtained, we calc obs CN for each cluster
# as exp (logRc - logR0) where logRc is the logR summary for the cluster

# once we have obs total number we obtain allele specific CNs m and p and 
# the cell fraction phi using logOR. For clonal CNAs, phi equals tumor puity
# For subconal events phi will be lower than overall sample purity

# Now Expectation-maximization EM algo is used for estimation refinement
# over all SNP loci and segment clusters. Bayes theorem to calculare posteriors
# iin the E-step of CN state g given parameters at the interaction. M-step, 
# given the imputed genotype, we update the model parameters by maximizing 
# the complete data likelihood. Iteration until convergence. X and Y adjusted. 

# segmentation output and EM fit
head (fit$cncf)
fit$purity
fit$ploidy

plt = plotSample(x=oo,emfit=fit)
# snakemake@output[["plot"]
png(filename = 'plot1', width = 1000, height=1500)
print (plt)
dev.off()

spiderplt = logRlogORspider(oo$out, oo$dipLogR)
# snakemake@output[["diagplot"]
png(filename = 'diagplot', width = 500, height=500)
print (spiderplt)
dev.off()

# Write code to export tables. 
