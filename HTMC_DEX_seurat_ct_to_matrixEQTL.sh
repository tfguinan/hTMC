#PBS -N HTMC_DEX_seurat_ct_to_matrixEQTL
#PBS -M thomas.guinan@utas.edu.au
#PBS -m abe
#PBS -l select=1:ncpus=8:mem=256gb
#PBS -l walltime=192:00:00
#PBS -o /data/menzies_projects/hewittlab/HTMC/share/logs/
#PBS -j oe
#PBS -J 1-2:1
export R_LIBS=/data/menzies_projects/onek1k/share/installs/packages
module load rosalind R

DIR=/data/menzies_projects/hewittlab/HTMC/share
OUTDIR=$DIR/output_data

# Make this an array job to match # of lines in Cell_type.csv
# Read a single line to produce per cell type matrices for PEER and mEQTL
CTCSV=/data/menzies_projects/hewittlab/HTMC/share/raw_data/HTMC_cell_types.csv
# Need to ensure LF encoding
CT=$(sed -n ${PBS_ARRAY_INDEX}p $CTCSV)

CTDIR=$OUTDIR/HTMC_DEX_${CT}
mkdir -p $CTDIR

# DEX
echo "Starting DEX file preparation"
# This script will have run prior (only runs to create object with all)
# COUNTDIR=$OUTDIR/HTMC_count_array_v015
# Rscript $DIR/main/v015/HTMC_aggr_seurat.R "${COUNTDIR}/HTMC_pool-dex?_GRCh38/*_cellbender/*.h5" $AGGOBJ
# Aggregate CellRanger Counts output
# Aggregated Seurat object
AGGOBJ=$OUTDIR/HTMC_DEX_merged.RDS


# Experiment defined covariates
COVIN=$DIR/raw_data/HTMC_sampleid_age_sex.csv
# PCA from genotyping
PCAIN=$DIR/genotyping/g1k/PCA.eigenvec
# plink --vcf .../HTMC_GRCh38.dose.vcf.gz --make-bed --out .../HTMC_GRCh38.dose --allow-extra-chr
BEDIN=$OUTDIR/HTMC_GRCh38.dose.bed

# Output file paths (cell type specific)
IDXPATH=$CTDIR/HTMC_DEX_${CT}_index.csv

# The path for PEER matrix (per donor sums of log1p, z-scaled SCT counts)
PMXPATH=$CTDIR/HTMC_DEX_${CT}_QC11drsum_mtx.csv
PCVPATH=$CTDIR/HTMC_DEX_${CT}_peer_covariates.csv
FACPATH=$CTDIR/HTMC_DEX_${CT}_peer_factors.csv
RESPATH=$CTDIR/HTMC_DEX_${CT}_peer_residuals.csv

# The path for MatrixEQTL matrix (per donor means of SCT counts)
MMXPATH=$CTDIR/HTMC_DEX_${CT}_SCTdrmean_mtx.csv
# The path for MatrixEQTL covariates (will get PFs added)
MCVPATH=$CTDIR/HTMC_DEX_${CT}_matrixEQTL_covariates.csv
SNPPATH=$CTDIR/HTMC_DEX_${CT}_SNPs.csv
QTLPATH=$CTDIR/HTMC_DEX_${CT}_eQTLs

# TODO reorder the arguments in the script
# $AGGOBJ $COVIN $PCAIN $BEDIN $IDXPATH $PMXPATH $PCVPATH $MMXPATH $SNPPATH 
Rscript $DIR/main/v015/HTMC_seurat_ct_to_peer.R \
    $AGGOBJ $IDXPATH $PMXPATH $COVIN $PCAIN $BEDIN $PCVPATH $SNPPATH $MMXPATH $CT


# PEER
echo "Starting DEX PEER"

module load rosalind Anaconda3
source activate peer

Rscript $DIR/main/v015/HTMC_peer.R $PMXPATH $PCVPATH $FACPATH $RESPATH

conda deactivate


# MatrixEQTL
echo "Starting DEX MatrixEQTL file preparation"

Rscript $DIR/main/v015/HTMC_covs_matrixEQTL.R $IDXPATH $PCVPATH $FACPATH $MCVPATH

Rscript $DIR/main/v015/HTMC_matrixEQTL.R $SNPPATH $MMXPATH $MCVPATH $QTLPATH
