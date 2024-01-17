# export R_LIBS="/data/menzies_projects/onek1k/share/installs/packages"
library(Matrix)
library(Seurat)
library(SeuratDisk)
library(scPred)
library(BEDMatrix)
args <- commandArgs(trailingOnly=TRUE)

# Loads seurat object (created by merging pools)
object <- readRDS(args[1])
# Subset for only singlets
object <- subset(object, subset = DROPLET.TYPE == 'SNG')
Idents(object) <- 'SNG.BEST.GUESS'

# Filter for only inds with >100 cells in both groups, passed genotyping
# This is temporary as we have retained all IDs... will remove upstream at later date
rm.ids <- c('HTMC-21', 'HTMC-22', 'HTMC-27', 'HTMC-33',
    'HTMC-35', 'HTMC-36', 'HTMC-43', 'HTMC-44', 'HTMC-46',
    'HTMC-51', 'HTMC-54', 'HTMC-58', 'HTMC-61', 'HTMC-62',
    'HTMC-65', 'HTMC-66', 'HTMC-71', 'HTMC-76')
rm.ids.ctl <- c('HTMC-21', 'HTMC-22', 'HTMC-27', 'HTMC-33',
    'HTMC-35', 'HTMC-36', 'HTMC-43', 'HTMC-44', 'HTMC-46',
    'HTMC-51', 'HTMC-54', 'HTMC-58', 'HTMC-61', 'HTMC-62',
    'HTMC-65', 'HTMC-66', 'HTMC-76')
subset.ids <- function(x) {tryCatch({subset(x, idents = rm.ids, invert = TRUE)}, 
    error = function(err) {subset(x, idents = rm.ids.ctl, invert = TRUE)})}
object <- subset.ids(object)
message('Subsetted by ID')

ct <- args[10]
object <- subset(object, subset = scpred_prediction == ct)
message('Subsetted by cell type')


# Preparing for PEER

# For PEER rows are sample and columns feature
# matrix <- AggregateExpression(object, assays = 'SCT')
matrix <- AverageExpression(object, assays = 'SCT', return.seurat=FALSE)
matrix <- matrix$SCT

# Save ordering index for later use
sample.order <- colnames(matrix)
names(sample.order) <- NULL
write.csv(sample.order, args[2])

message("Saving matrix")
message(args[9])
# Here we save MatrixEQTL so rows are feature, columns are sample
write.table(matrix, args[9], sep = ',', row.names = TRUE, col.names = TRUE)

# QC steps following QC #11 from: 
# Pitfalls and opportunities for applying latent variables in single-cell eQTL analyses
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-02873-5

# Transpose so that rows are sample, columns are genes for QC and PEER
matrix <- t(matrix)

pi0 <- 0.9
# Filter zero expression
pistat <- colSums(matrix==0, na.rm=TRUE)/sum(!is.na(matrix[,1]))
matrix <- matrix[,which(pistat < pi0)]
# Log normalisation
matrix <- apply(matrix, 2, function(x)log(x+1))
# Scale and centre
matrix <- scale(matrix, scale=TRUE, center=TRUE)

write.table(matrix, args[3], sep = ',', row.names = FALSE, col.names = FALSE)


# Covariate data
covs <- read.csv(args[4])
rownames(covs) <- covs$sampleid
# Save only age and sex in correct order
covs <- covs[,2:3]

# Load genotype pca data
pca.path <- args[5]
gpca <- read.table(pca.path, sep = '\t', header = TRUE, comment.char="")
rownames(gpca) <- gpca[,1]
gpca <- gpca[2:ncol(gpca)]

gpca <- gpca[rownames(gpca) %in% rownames(covs),]
# gpca <- gpca[match(sample.order,rownames(gpca)),]
covs <- merge(covs, gpca[1:6], by=0, all=TRUE)
rownames(covs) <- covs$Row.names
covs <- covs[2:ncol(covs)]

# This now has age, sex, six genotyping PCs ordered to match the averaged matrix
message("Peek at ordering covariates matrix")
message(head(rownames(covs)))
covs <- covs[match(sample.order,rownames(covs)),]
message(head(rownames(covs)))
# Next step is to include PFs

# PEER isn't accepting covariates it seems...
write.table(covs, args[7], sep = ',', row.names = FALSE, col.names = FALSE)



# Preparing for MatrixEQTL

# Mean of donor aggregated SCT adjusted counts

# Slot arg is deprecated using layer instead




# BED file generated from vcf by:
# plink --vcf .../HTMC_GRCh38.dose.vcf.gz --make-bed --out .../HTMC_GRCh38.dose --allow-extra-chr
snps.file <- args[6]
snps <- BEDMatrix(snps.file)
# Have to reformat the IDs from HTMC-X_HTMC-X to just HTMC-X
# Should work for const-fid too i.e 0_HTMC-X
rownames(snps) <- unlist(lapply(strsplit(rownames(snps), '_'), function(x) return(x[[2]])))
# Correct order
message("Peek at ordering snps matrix")
message(head(rownames(snps)))
snps <- snps[match(sample.order,rownames(snps)),]
message(head(rownames(snps)))
# Must transpose such that columns are samples
snps <- t(as.matrix(snps))

write.table(snps, args[8], sep = ',', row.names = TRUE, col.names = TRUE)


