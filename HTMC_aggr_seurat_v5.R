# export R_LIBS="/data/menzies_projects/onek1k/share/installs/packages"
library(Matrix)
library(Seurat)
library(SeuratDisk)
library(scPred)
library(harmony)


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


# Modified project_query() code from:
# https://github.com/powellgenomicslab/scPred/blob/master/R/project_query.R
# Update call of GetAssayData() to use layer not assay 'data'
scpredmod.project_query <- function(new, 
                          reference, 
                          max.iter.harmony = 20,
                          recompute_alignment = TRUE,
                          seed = 66, 
                          ...){


# Validate if provided object is an scPred object
if(!(is(reference, "Seurat") | is(reference, "scPred"))) stop("'object' must be of class 'scPred' or 'Seurat'")

if(is(reference, "Seurat")){
  spmodel <- reference@misc$scPred
}else{
  spmodel <- reference
}

if(is.null(spmodel)) stop("No feature space has been determined!")
if(!is(new, "Seurat")) stop("New data must be a Seurat object")


# Dataset alignment -------------------------------------------------------

if("scpred" %in% names(new@reductions)){
  if(recompute_alignment){
    alignment <- TRUE
    cat(crayon::yellow(cli::symbol$figure_dash, "Data has already been aligned to a reference.\n"), sep = "")
    cat(crayon::yellow(cli::symbol$sup_plus, "Skip data alignment using `recompute.alignment = FALSE`.\n"),  sep = "")
      
  } 
  else {
    alignment <- FALSE
  }
  
}else{
  alignment <- TRUE
}

if(alignment){
  
  cat(crayon::green(cli::symbol$record, " Matching reference with new dataset...\n"))
  
  # Subset data
  ref_loadings <- spmodel@feature_loadings
  ref_embeddings <- spmodel@cell_embeddings
  new_features <- rownames(new)
  
  # Get genes
  reference_features <- rownames(ref_loadings)
  
  
  # Get intersection between reference and new datasets
  shared_features <- intersect(reference_features, new_features)
  cat(crayon::cyan("\t", cli::symbol$line, paste(length(reference_features), "features present in reference loadings\n")))
  cat(crayon::cyan("\t", cli::symbol$line, paste(length(shared_features), "features shared between reference and new dataset\n")))
  cat(crayon::cyan("\t", cli::symbol$line, paste0(round(length(shared_features)/length(reference_features) * 100, 2), 
                                                  "% of features in the reference are present in new dataset\n")))
  
  
  # Subset shared genes from reference
  ref_loadings <- ref_loadings[shared_features, ]
  
  
  
  new_data <- GetAssayData(new, assay="RNA", layer="data")[shared_features,]
  means <- spmodel@scaling$means
  stdevs  <- spmodel@scaling$stdevs
  new_data <- Matrix::t(new_data)
  names(means) <- rownames(spmodel@scaling) -> names(stdevs)
  
  # Subset means and standard deviations
  means <- means[shared_features]
  stdevs <- stdevs[shared_features]
  #all(colnames(new_data) == names(means))
  
  i <- stdevs == 0
  
  if(any(i)){
    warning(paste0(sum(i), " features have zero variance but are present in the feature loadings. \nDid you subset or integrated this data before?"))
    cat(crayon::yellow("Removing zero-variance genes from projection\n"))
    
    new_data <- new_data[,!i]
    ref_loadings <- ref_loadings[!i,]
    means <- means[!i]
    stdevs <- stdevs[!i]
    #all(colnames(new_data) == rownames(ref_loadings))
  }
  
  
  scaled_data <- scale(new_data, means, stdevs)
  
  
  new_embeddings <- scaled_data %*% ref_loadings
  
  
  
  dataset <- factor(c(rep("reference", nrow(ref_embeddings)), rep("new", nrow(new_embeddings))), 
                    levels = c("reference", "new"))
  
  
  rownames(ref_embeddings) <- paste0("ref_", rownames(ref_embeddings))
  rownames(new_embeddings) <- paste0("new_", rownames(new_embeddings))
  
  
  eigenspace <- as.data.frame(rbind(ref_embeddings, new_embeddings))
  meta_data <- data.frame(rownames(eigenspace), dataset = dataset)
  
  cat(crayon::green(cli::symbol$record, " Aligning new data to reference...\n"))
  
  set.seed(seed)
  harmony_embeddings <- HarmonyMatrix(eigenspace, 
                                      meta_data, 
                                      'dataset', 
                                      do_pca = FALSE, 
                                      reference_values = "reference",
                                      max.iter.harmony = max.iter.harmony, 
                                      ...)
  
  new_embeddings_aligned <- harmony_embeddings[dataset == "new", , drop = FALSE]
  
}else{
  new_embeddings_aligned <- Embeddings(new, reduction = "scpred")
  colnames(new_embeddings_aligned) <- gsub("scpred_", spmodel@reduction_key, colnames(new_embeddings_aligned))
}

  
  rownames(new_embeddings_aligned) <- gsub("^new_", "", rownames(new_embeddings_aligned))
  new@reductions[["scpred"]] <- CreateDimReducObject(embeddings = new_embeddings_aligned, 
                                                     key = "scpred_",
                                                     assay = DefaultAssay(object = new))
  if(recompute_alignment){
    
    rownames(new_embeddings) <- gsub("^new_", "", rownames(new_embeddings))
    new@reductions[["scpred_projection"]] <- CreateDimReducObject(embeddings = new_embeddings, 
                                                                  key = "Projection_",
                                                                  assay = DefaultAssay(object = new))
  }

  new
  
}

# Modified scPredict() code from:
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

  new <- scpredmod.project_query(new, 
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
# glob.result <- Sys.glob(args[1])
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
    temp.pool <- seuratmod.read_h5(pool, genome.name='matrix')
    temp.pool <- CreateSeuratObject(counts=temp.pool, project=basename(dirname(dirname(pool))))
    temp.pool <- subset(temp.pool, subset=nFeature_RNA > 0)
    # Apply scPred model(s)
    temp.pool <- NormalizeData(temp.pool)
    temp.pool <- scpredmod.scPredict(temp.pool, trained.scpred.1, name='1')
    temp.pool <- scpredmod.scPredict(temp.pool, trained.scpred.2, name='2')

    # Add demuxlet data
    best.path <- Sys.glob(paste0(dirname(dirname(pool)), '/*_demuxlet_default_allbc/*_dmx.best'))
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

    pool.vector <- append(pool.vector, temp.pool) 
}

# Removed code
# message(length(colnames(temp.pool)))
# temp.pool[["percent.mt"]] <- PercentageFeatureSet(temp.pool, pattern = "^MT-")


# Further normalisation
# temp.pool <- ScaleData(temp.pool)
# temp.pool <- SCTransform(temp.pool, vars.to.regress = "percent.mt", verbose = FALSE)

# This function will return 'A' from HTMC_pool-A_GRCh38, 'dexA' from HTMC_pool-dexA_GRCh38
get.name <- function(obj) {
    return(strsplit(strsplit(as.character(obj$orig.ident[1]), '-')[[1]][2], '_')[[1]][1])
}

project.names <- sapply(pool.vector, get.name)
# Merge all objects to one
merged.pool <- merge(pool.vector[[1]], pool.vector[2:length(pool.vector)], add.cell.ids=project.names, merge.data=TRUE)


# SCTransform
merged.pool <- SCTransform(merged.pool)
merged.pool <- RunPCA(merged.pool, npcs=50,verbose=FALSE)
options(future.globals.maxSize = 4e+10)
integrated.merged.pool <- IntegrateLayers(object=merged.pool, method=RPCAIntegration, normalization.method="SCT", verbose=FALSE)

# Not sure this DR will be present
integrated.merged.pool <- FindNeighbors(integrated.merged.pool, dims=1:50, reduction='integrated.dr')
integrated.merged.pool <- FindClusters(integrated.merged.pool, resolution=2)

integrated.merged.pool <- RunUMAP(integrated.merged.pool, dims = 1:50, reduction='integrated.dr')

table(integrated.merged.pool$SNG.BEST.GUESS, useNA='always')
table(integrated.merged.pool$DROPLET.TYPE, useNA='always')
table(integrated.merged.pool$X1_scpred_prediction)
table(integrated.merged.pool$X2_scpred_prediction)

saveRDS(merged.pool, args[2])

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