# source activate peer
library(peer)
# library(Matrix)
args <- commandArgs(trailingOnly=TRUE)

# Need to load including row and column names, then set to empty for PEER to accept matrix
mtx <- as.matrix(read.csv(args[1], header=FALSE))
dim(mtx)
dimnames(mtx) <- NULL
# mtx <- as.matrix(read.csv("/data/menzies_projects/hewittlab/HTMC/share/output_data/HTMC_CTL_indv_SCTmatrix_QC11.csv", header=TRUE))

model <- PEER()
PEER_setPhenoMean(model, mtx)

# Covariates included here in correct order to match per individual aggregated expression matrix
covs <- as.matrix(read.csv(args[2], header=FALSE))
dim(covs)
dimnames(covs) <- NULL
# covs <- as.matrix(read.csv("/data/menzies_projects/hewittlab/HTMC/share/output_data/HTMC_CTL_covs.csv", header=TRUE))
# covs[,2] <- as.numeric(gsub("M", 1, gsub("F", 2, covs[,2])))

# Trying some other modifications
# covs<-sampleinfo[,c("Sex","Age","BMI","Smoking","hypertension","hyperlipidemia","Diabetes")]
# covs<-apply(covs,2,as.numeric)
# rownames(covs)<-sampleinfo$ID


# TODO figure out including covariates
# Will try (12/10/24)
PEER_setCovariates(model, as.matrix(covs))

# Update variable iteration
PEER_setNmax_iterations(model, 2000)

# Include mean effect
# PEER_setAdd_mean(model, TRUE)

# Set number of peer factors
num.pf <- 50
PEER_setNk(model, num.pf)

# These are default bound, variance tolerances
# PEER_setTolerance(model, 0.001)
# PEER_setVarTolerance(model, 0.00001)

PEER_update(model)

# The posterior mean of inferred confounders (NxK ... + known covariates)
factors = PEER_getX(model)
write.csv(factors, args[3])
# weights = PEER_getW(model)
# precision = PEER_getAlpha(model)

residuals = PEER_getResiduals(model)
write.csv(factors, args[4])