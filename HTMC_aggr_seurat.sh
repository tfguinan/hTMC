#PBS -N HTMC_merge
#PBS -M tfguinan@utas.edu.au
#PBS -m abe
#PBS -l select=1:ncpus=1:mem=72gb
#PBS -l walltime=192:00:00
#PBS -o /data/menzies_projects/hewittlab/HTMC/share/logs/
#PBS -j oe


module load rosalind R
export R_LIBS=/data/menzies_projects/onek1k/share/installs/packages
AGGR_R=/data/menzies_projects/hewittlab/HTMC/share/main/v015/HTMC_aggr_seurat.R
DIRECTORY=/data/menzies_projects/hewittlab/HTMC/share/output_data/HTMC_count_array_v015
OUTDIR=/data/menzies_projects/hewittlab/HTMC/share/output_data


# CTL
Rscript $AGGR_R "${DIRECTORY}/HTMC_pool-?_GRCh38/*_cellbender/*.h5" "${OUTDIR}/HTMC_CTL_2_merged.RDS"

# DEX
Rscript $AGGR_R "${DIRECTORY}/HTMC_pool-dex?_GRCh38/*_cellbender/*.h5" "${OUTDIR}/HTMC_DEX_2_merged.RDS"