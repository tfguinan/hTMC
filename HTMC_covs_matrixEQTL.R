args <- commandArgs(trailingOnly=TRUE)

# This script adds PFs to our covariates (after PEER completes)
# NOTE: We can now include our covariates in PEER
# We now just reduce the number of columns (for 25PFs + 8 others)
# Really only need args 1 and 3

sample.order <- read.csv(args[1], row.names=1, header=TRUE)$x
covs <- read.csv(args[2], row.names=NULL, header=FALSE)
pfs <- read.csv(args[3], row.names=1, header=TRUE)

# These should follow sample order as specified prior to PEER
# pfs <- pfs[2:ncol(pfs)]
rownames(pfs) <- sample.order
# Redundant as first 8 are not PFs
# colnames(pfs) <- sapply(seq(1, ncol(pfs)), function(x) paste0('PF',x))

pfs <- pfs[,1:33]

# See above note
# rownames(covs) <- sample.order
# covs <- merge(covs, pfs[1:25], by=0, all=TRUE)
# covs <- covs[2:ncol(covs)]

# We must transpose such that rows are features columns are sample
write.table(t(pfs), args[4], sep = ',', row.names = TRUE, col.names = TRUE)
