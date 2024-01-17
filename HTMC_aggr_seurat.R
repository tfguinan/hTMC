# export R_LIBS="/data/menzies_projects/onek1k/share/installs/packages"
library(Matrix)
library(Seurat)
library(SeuratDisk)
library(scPred)
library(ggrastr)


# Modified h5 loader code from:
# https://github.com/satijalab/seurat/blob/763259d05991d40721dee99c9919ec6d4491d15e/R/preprocessing.R
# Permitting access to only 'matrix' from cellbender h5
seuratmod.read_h5 <- function(filename, use.names = TRUE, unique.features = TRUE, genome.name = NULL) {
  if (!requireNamespace('hdf5r', quietly = TRUE)) {
    stop("Please install hdf5r to read HDF5 files")
  }
  if (!file.exists(filename)) {
    stop("File not found")
  }
  infile <- hdf5r::H5File$new(filename = filename, mode = 'r')


  output <- list()
  if (hdf5r::existsGroup(infile, 'matrix')) {
    # cellranger version 3
    if (use.names) {
      feature_slot <- 'features/name'
    } else {
      feature_slot <- 'features/id'
    }
  } else {
    if (use.names) {
      feature_slot <- 'gene_names'
    } else {
      feature_slot <- 'genes'
    }
  }

  if (genome.name == FALSE) {
    genomes <- names(x = infile)
  } else {
    genomes <- genome.name
  }

  for (genome in genomes) {
    counts <- infile[[paste0(genome, '/data')]]
    indices <- infile[[paste0(genome, '/indices')]]
    indptr <- infile[[paste0(genome, '/indptr')]]
    shp <- infile[[paste0(genome, '/shape')]]
    features <- infile[[paste0(genome, '/', feature_slot)]][]
    barcodes <- infile[[paste0(genome, '/barcodes')]]
    sparse.mat <- sparseMatrix(
      i = indices[] + 1,
      p = indptr[],
      x = as.numeric(x = counts[]),
      dims = shp[],
      repr = "T"
    )
    if (unique.features) {
      features <- make.unique(names = features)
    }
    rownames(x = sparse.mat) <- features
    colnames(x = sparse.mat) <- barcodes[]
    sparse.mat <- as.sparse(x = sparse.mat)
    # Split v3 multimodal
    if (hdf5r::existsGroup(infile, paste0(genome, '/features'))) {
      types <- infile[[paste0(genome, '/features/feature_type')]][]
      types.unique <- unique(x = types)
      if (length(x = types.unique) > 1) {
        message(
          "Genome ",
          genome,
          " has multiple modalities, returning a list of matrices for this genome"
        )
        sparse.mat <- sapply(
          X = types.unique,
          FUN = function(x) {
            return(sparse.mat[which(x = types == x), ])
          },
          simplify = FALSE,
          USE.NAMES = TRUE
        )
      }
    }
    output[[genome]] <- sparse.mat
  }
  infile$close_all()
  if (length(x = output) == 1) {
    return(output[[genome]])
  } else{
    return(output)
  }
}


# Modified scpredict() code from:
# https://github.com/powellgenomicslab/scPred/blob/master/R/scPredict.R
# Permitting unique name for multiple scPred references in a single seurat object
scpredmod.scPredict <- function(new,
                      reference, 
                      threshold = 0.55, 
                      max.iter.harmony = 20,
                      recompute_alignment = TRUE,
                      seed = 66,
                      name = NA){
  
  # Function validations ----------------------------------------------------
  
  # Validate if provided object is an scPred object
  if(!(is(reference, "Seurat") | is(reference, "scPred"))) stop("'object' must be of class 'scPred' or 'Seurat'")
  
  if(is(reference, "Seurat")){
    spmodel <- reference@misc$scPred
  }else{
    spmodel <- reference
  }
  
  if(is.null(spmodel)) stop("No feature space has been determined!")
  if(!length(spmodel@train)) stop("No models have been trained!")
  if(!is(new, "Seurat")) stop("New data must be a Seurat object")
  

  # Project query data ------------------------------------------------------

  new <- project_query(new, 
                reference = spmodel, 
                max.iter.harmony = max.iter.harmony, 
                recompute_alignment = recompute_alignment,
                seed = seed)
  
  new_embeddings_aligned <- Embeddings(new[["scpred"]])
  colnames(new_embeddings_aligned) <- colnames(spmodel@cell_embeddings)
  
  
  
  # Classify cells using all trained models 
  cellTypeModelNames <- names(spmodel@features)
  .predictCellClass <- function(cellType, spmodel, testEmbeddings){
    
    # Extract features for a given cell type
    as.character(spmodel@features[[cellType]]$feature) -> features
    
    # Extract cell type model
    model <- spmodel@train[[cellType]]
    
    # Perform predictions based on informative PCs
    prediction <- predict(model, 
                          newdata = scPred:::subsetMatrix(testEmbeddings, features), 
                          type = "prob")
    
    # Add cell names to results
    rownames(prediction) <- rownames(testEmbeddings)
    
    # Return positive-class probability
    prediction[,1, drop = FALSE]
    
  }
  
  cat(crayon::green(cli::symbol$record, " Classifying cells...\n"))
  res <- sapply(cellTypeModelNames, .predictCellClass, spmodel, new_embeddings_aligned)
  
  # Gather results
  res <- as.data.frame(res)
  colnames(res) <- cellTypeModelNames
  rownames(res) <- colnames(new)
  
  classes <- cellTypeModelNames
  #plot(res$Lymphoid, col = as.factor(test$CellType))
  # If there is only 2 classes, compute complementary probability for negative class
  if(length(cellTypeModelNames) == 1){
    metadata <- get_metadata(spmodel)
    cellClasses <- levels(metadata$pvar)
    res_comp <- 1 - res[,1]
    negClass <- cellClasses[cellClasses != names(res)]
    res[[negClass]] <- res_comp
    
  }
  
  # Extract maximum probability for each class
  max_props <- as.data.frame(t(apply(res, 1, function(x) c(index = which.max(x),  max = x[which.max(x)]))))
  names(max_props) <- c("index", "max")
  
  # Store classification based on maximum probability
  max_props$generic_class <- names(res)[max_props$index]
  res <- cbind(res, max_props)
  
  # Classify cells according to probability threshold
  pred <- ifelse(res$max > threshold, res$generic_class, "unassigned")
  
  names(pred) <- colnames(new)
  
  # Format results
  res$prediction <- pred
  res$index <- NULL
  res$no_rejection <- res$generic_class
  res$generic_class <- NULL
  
  if (is.na(name)){
    names(res) <- .make_names(paste0("scpred_", names(res)))
  } else{
    names(res) <- .make_names(paste0(name, "_scpred_", names(res)))
  }

  # Return results
  new <- AddMetaData(new, res)
  
  cat(crayon::green("DONE!\n"))
  
  new
}
# We also require .makenames from:
# https://github.com/powellgenomicslab/scPred/blob/master/R/utils.R
.make_names <- function(x){
  x <- gsub("\\+", "_plus", x)
  x <- gsub("\\-", "_minus", x)
  x <- make.names(x)
}


args <- commandArgs(trailingOnly=TRUE)
# args <- c("/data/menzies_projects/hewittlab/HTMC/share/output_data/HTMC_count_array_v015/HTMC_pool-?_GRCh38/*_cellbender/*_filtered.h5")

message("Glob for:")
# TODO look at loading h5 as anndata object to use ambient expression?
message(args[1])
# Locate raw 10x matrix directory using default name from command arg
glob.result <- Sys.glob(args[1])
message(glob.result)

# TODO provide code for training of scPred model
model.1.path <- '/data/menzies_projects/hewittlab/HTMC/share/van_zyl_2020_scpred-fixed.rds'
model.2.path <- '/data/menzies_projects/hewittlab/HTMC/share/van_zyl_2022_model.rds'

trained.scpred.1 <- readRDS(model.1.path)
trained.scpred.2 <- readRDS(model.2.path)



pool.vector <- c()
# Load pool counts, create seurat object, and append to vector
# Load demuxlet information
for (pool in glob.result){
    message(pool)
    temp.pool <- seuratmod.read_h5(pool, genome.name = 'matrix')
    temp.pool <- CreateSeuratObject(counts=temp.pool, project=basename(dirname(dirname(pool))))
    message(length(colnames(temp.pool)))
    temp.pool[["percent.mt"]] <- PercentageFeatureSet(temp.pool, pattern = "^MT-")
    temp.pool <- subset(temp.pool, subset = nFeature_RNA > 100 & percent.mt < 5)
    message(length(colnames(temp.pool)))
    
    # Apply scPred model(s)
    temp.pool <- NormalizeData(temp.pool)
    temp.pool <- scpredmod.scPredict(temp.pool, trained.scpred.1, name='1')
    temp.pool <- scpredmod.scPredict(temp.pool, trained.scpred.2, name='2')

    # Add demuxlet data
    best.path <- Sys.glob(paste0(dirname(dirname(pool)), '/*_demuxlet/*_dmx.best'))
    best <- read.csv(best.path, sep = '\t')
    # Droplet type
    best.type <- data.frame(best$DROPLET.TYPE)
    rownames(best.type) <- c(best$BARCODE)
    temp.pool <- AddMetaData(temp.pool, metadata=best.type, col.name='DROPLET.TYPE')
    # Singlet identity prediction
    best <- best[best$DROPLET.TYPE == 'SNG',]
    rownames(best) <- c(best$BARCODE)
    best <- best['SNG.BEST.GUESS']
    temp.pool <- AddMetaData(temp.pool, metadata=best, col.name='SNG.BEST.GUESS')

    # Further normalisation
    temp.pool <- ScaleData(temp.pool)
    temp.pool <- SCTransform(temp.pool, vars.to.regress = "percent.mt", verbose = FALSE)

    pool.vector <- append(pool.vector, temp.pool) 
}

# This function will return 'A' from HTMC_pool-A_GRCh38, 'dexA' from HTMC_pool-dexA_GRCh38
get.name <- function(obj) {
    return(strsplit(strsplit(as.character(obj$orig.ident[1]), '-')[[1]][2], '_')[[1]][1])
}

project.names <- sapply(pool.vector, get.name)
# Merge all objects to one
pool.merging <- merge(pool.vector[[1]], pool.vector[2:length(pool.vector)], add.cell.ids=project.names, merge.data=TRUE)

saveRDS(pool.merging, args[2])

# pdf("/data/menzies_projects/hewittlab/HTMC/share/output_data/interactive_out/vlnplot.pdf")
# VlnPlot(temp.pool, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# dev.off()

# pdf("/data/menzies_projects/hewittlab/HTMC/share/output_data/interactive_out/scatter.pdf", width=20, height=10)
# plot1 <- FeatureScatter(temp.pool, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(temp.pool, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2
# dev.off()

# Pre QC plots

# QC
# Will use cellbender output to determine min RNA cutoff
# Will include demuxlet output to label cells per individual

# HVG selection (can be input to PEER analysis)
# Linear dimensional reduction
# Cluster (louvain)
# Non linear dimensional reduction (UMAP, tSNE)
# DE
# Post QC plots


# TODO investigate retraining a model based on TM tissues only
# temp.pool[[c('X1_scpred_prediction', 'X2_scpred_prediction')]][temp.pool[['X1_scpred_prediction']] == 'Fibroblast']