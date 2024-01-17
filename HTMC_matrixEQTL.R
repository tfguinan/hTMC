library(MatrixEQTL)
args <- commandArgs(trailingOnly=TRUE)

out.file <- args[4]

# Genotyping data
snps.file <- args[1]
snps = SlicedData$new()
snps$fileDelimiter <- ","
snps$fileOmitCharacters <- ""
snps$fileSkipRows <- 1
snps$fileSkipColumns <- 1
snps$fileSliceSize <- 2000
snps$LoadFile(snps.file)

# Expression data
expr.file <- args[2]
expr = SlicedData$new()
expr$fileDelimiter <- ","
expr$fileOmitCharacters <- ""
expr$fileSkipRows <- 1
expr$fileSkipColumns <- 1
expr$fileSliceSize <- 2000
expr$LoadFile(expr.file)

# Covariate data
covs.file <- args[3]
covs = SlicedData$new()
covs$fileDelimiter <- ","
covs$fileOmitCharacters <- ""
covs$fileSkipRows <- 1
covs$fileSkipColumns <- 1
covs$fileSliceSize <- 2000
covs$LoadFile(covs.file)


me = Matrix_eQTL_engine(
    snps = snps,
    gene = expr,
    cvrt = covs,
    output_file_name = out.file,
    pvOutputThreshold = 1e-2,
    useModel = modelLINEAR,
    errorCovariance = numeric(),
    verbose = TRUE,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE)

