#PBS -N HTMC_count_array_v015
#PBS -M tfguinan@utas.edu.au
#PBS -m abe
#PBS -l select=1:ncpus=32:mem=72gb
#PBS -l walltime=192:00:00
#PBS -o /data/menzies_projects/hewittlab/HTMC/share/logs/HTMC_count_array_v015/
#PBS -j oe
#PBS -J 1-18:1
VERSION=0.1.5

FULL_TIME="$(date +%s)" 

# Per each pool, this script runs: CellRanger Count; CellBender; Scrublet; SCDS; and (on Scrublet ∩ SCDS singlet barcodes) Demuxlet

### !!! FILL HERE !!! ###

THREADS=32
NAME=HTMC_count_array_v015
# This script requires a main project directory (dir), we recommend within this to have raw_data and output_data
DIR=/data/menzies_projects/hewittlab/HTMC/share
MAIN_CSV=${DIR}/raw_data/HTMC_count_array_v015.csv
FASTQ_DIR=${DIR}/raw_data/fastqs
# Reads aligned to 10x provided GRCh38
FASTA=${DIR}/raw_data/10x_GRCh38/refdata-gex-GRCh38-2020-A/fasta/genome.fa
# Our genotyping vcf
GENO_VCF=${DIR}/raw_data/HTMC_R8M05_protein_sorted_GRCh38.dose.vcf.gz
FIELD=GT


### !!! FILL HERE !!! ###


echo "##### HTMC COUNT ARRAY SCRIPT V${VERSION} #####"
module load rosalind gcc-env

# Cellranger install location
CR_BIN=${DIR}/cellranger-7.1.0/bin
# Helper tool install location
HELPER_DIR=${DIR}/popscle_helper_tools
# Demuxafy install location
DEMUXAFY_DIR=${DIR}/demuxafy-2.0.1


# Read in pool name and sample list path
LINE=$(awk "NR==$PBS_ARRAY_INDEX" "$MAIN_CSV")
IFS=',' read -r POOL SAMPLE INDS NCELLS NGEMS MMAC <<< "$LINE"


##### 1 #####
##### Cellranger count #####
echo && echo "##### CELLRANGER COUNT (RUN1+2) #####"
echo && echo $POOL
echo $SAMPLE

# Initial directory setup
cd $DIR
OUT_DIR=${DIR}/output_data/${NAME}/${POOL}
mkdir -p $OUT_DIR
cd $OUT_DIR
CR_OUT_DIR=${OUT_DIR}/${POOL}_count/outs

date
SECONDS=0

if [ -e $CR_OUT_DIR ]; then
    echo "Cellranger count has run with output"
else
    ${CR_BIN}/cellranger count --id=$POOL --fastqs=$FASTQ_DIR --sample=$SAMPLE \
    --transcriptome=$(dirname $(dirname $FASTA)) --disable-ui
fi

duration=$SECONDS
echo "${duration} seconds elapsed running step"
date
qstat -f $PBS_JOBID | grep 'resources_used.mem\|exec_vnode'



##### 2 #####
##### Cellbender #####
echo && echo "##### CELLBENDER #####"

CELLBENDER_OUT_DIR=${OUT_DIR}/${POOL}_cellbender
mkdir -p $CELLBENDER_OUT_DIR

CELLBENDER_IN=${CR_OUT_DIR}/raw_feature_bc_matrix.h5
CELLBENDER_OUT=${CELLBENDER_OUT_DIR}/${POOL}.h5

module load Anaconda3
# Activate conda using profile (workaround conda init bash)
. "$(conda info --base)/etc/profile.d/conda.sh" && conda activate CellBender

date
SECONDS=0

if [ -e $CELLBENDER_OUT ]; then
    echo "Cellbender has run with output"
else
    echo "Cellbender version: " && cellbender --version
    cellbender remove-background --input $CELLBENDER_IN \
        --output $CELLBENDER_OUT --cpu-threads $THREADS \
        --epochs 200 --checkpoint-mins 180
fi

duration=$SECONDS
echo "${duration} seconds elapsed running step"
date
qstat -f $PBS_JOBID | grep 'resources_used.mem\|exec_vnode' # Get memory info

conda deactivate

# We select the cell bender filtered file with posteriorP(droplet=cell) > 50%
CELLBENDER_H5=${CELLBENDER_OUT_DIR}/${POOL}_filtered.h5
CELLBENDER_BARCODES=${CELLBENDER_OUT_DIR}/${POOL}_cell_barcodes.csv



##### 3 #####
##### Scrublet #####
echo && echo "##### SCRUBLET #####"

SCRUBLET_OUT_DIR=${OUT_DIR}/${POOL}_scrublet
mkdir -p $SCRUBLET_OUT_DIR

# We install and run scrublet outside of demuxafy container (illegal instruction error, unable to rollback)
# Demuxafy scrublet script updated to include explicit scanpy import
module load Python
# pip install scrublet && pip install --upgrade annoy==1.16.3
# wait && pip show scrublet
SCRUBLET_BARCODES=${SCRUBLET_OUT_DIR}/scrublet_singlet_threshold_barcodes.tsv

# TODO evidence to support r,v,k,t args
date
SECONDS=0

if [ -e $SCRUBLET_BARCODES ]; then
    echo "Scrublet has run with output"
else
    python ${DIR}/demuxafy-2.0.1/scrublet_dev.py \
       -m $CELLBENDER_H5 -o $SCRUBLET_OUT_DIR -b $CELLBENDER_BARCODES \
       -r 5 -v 25 -k 50 -t 0.300
fi

duration=$SECONDS
echo "${duration} seconds elapsed running step"
date
qstat -f $PBS_JOBID | grep 'resources_used.mem\|exec_vnode' # Get memory info



##### 4 #####
##### Scds #####
echo && echo "##### SCDS #####"

SCDS_OUT_DIR=${OUT_DIR}/${POOL}_scds
mkdir -p $SCDS_OUT_DIR

module load R
# We require a fixed library location # TODO update to include append user dir
export R_LIBS="/data/menzies_projects/onek1k/share/installs/packages"

# Cellbender gives non-10x h5 file, we convert in modified Scds R script (scds_dev.R)
SCDS_BARCODES=${SCDS_OUT_DIR}/scds_singlet_barcodes.tsv

date
SECONDS=0

if [ -e $SCDS_BARCODES ]; then
    echo && echo "Scds has run with output" && echo
else
    Rscript ${DIR}/demuxafy-2.0.1/scds_dev.R -i $CELLBENDER_H5 \
    -b $CELLBENDER_BARCODES -o $SCDS_OUT_DIR -s ${SCDS_OUT_DIR}/${POOL}_sce.rds
fi

duration=$SECONDS
echo "${duration} seconds elapsed running step"
date
qstat -f $PBS_JOBID | grep 'resources_used.mem\|exec_vnode' # Get memory info


##### 5 #####
##### Barcode intersect #####
echo && echo "##### BARCODE INTERSECT #####"

# Use our python script to get intersect barcodes
date
SECONDS=0

python ${DIR}/main/quick_intersect.py $SCRUBLET_BARCODES $SCDS_BARCODES

duration=$SECONDS
echo "${duration} seconds elapsed running step"
date
qstat -f $PBS_JOBID | grep 'resources_used.mem\|exec_vnode' # Get memory info

INTERSECT_BARCODES_CSV=${OUT_DIR}/scds_scrublet_intersect_barcodes.csv
INTERSECT_BARCODES_TSV=${OUT_DIR}/scds_scrublet_intersect_barcodes.tsv
module purge



##### 6 #####
##### Demuxlet #####

echo && echo "##### DEMUXLET #####"

DEMUXLET_OUT_DIR=${OUT_DIR}/${POOL}_demuxlet
mkdir -p $DEMUXLET_OUT_DIR

module load rosalind gcc-env bcftools bedtools samtools singularity
# TODO include Cellbender, Scrublet, and Scds (Non-empty, likely singlet)

# # Filter the 10x bam file based on barcodes that pass Cellbender (Non-empty), Scrublet and Scds (likely singlet)
BAM=${CR_OUT_DIR}/possorted_genome_bam.bam
FILTERED_BAM=${CR_OUT_DIR}/${POOL}_possorted_genome_filtered.bam
${HELPER_DIR}/filter_bam_file_for_popscle_dsc_pileup.sh $BAM $INTERSECT_BARCODES_TSV $GENO_VCF $FILTERED_BAM

# Use default min-mac arg (could specify 5 to match our $GENO_VCF filtering)
# TODO evidence to support geno-error args
date
SECONDS=0

if [ -e ${DEMUXLET_OUT_DIR}/${POOL}_dmx ]; then
    echo && echo "Demuxlet has run with output" && echo
else
    date
    singularity exec --bind ${DIR} $DEMUXAFY_DIR/Demuxafy.sif \
        popscle dsc-pileup --sam $FILTERED_BAM --vcf $GENO_VCF --group-list $INTERSECT_BARCODES_TSV --sm-list $INDS --out ${DEMUXLET_OUT_DIR}/${POOL}_plp
    echo
    date
    qstat -f $PBS_JOBID | grep 'resources_used.mem\|exec_vnode'

    singularity exec --bind ${DIR} $DEMUXAFY_DIR/Demuxafy.sif \
        popscle demuxlet --plp ${DEMUXLET_OUT_DIR}/${POOL}_plp --vcf $GENO_VCF --group-list $INTERSECT_BARCODES_TSV --sm-list $INDS --field $FIELD \
        --geno-error-coeff 1.0 --geno-error-offset 0.05 --out ${DEMUXLET_OUT_DIR}/${POOL}_dmx
    echo
    date
    qstat -f $PBS_JOBID | grep 'resources_used.mem\|exec_vnode'

    singularity exec --bind ${DIR} $DEMUXAFY_DIR/Demuxafy.sif \
        bash Demuxlet_summary.sh ${DEMUXLET_OUT_DIR}/${POOL}_dmx.best > ${DEMUXLET_OUT_DIR}/${POOL}_demuxlet.txt  
fi

duration=$SECONDS
echo "${duration} seconds elapsed running step"
date
qstat -f $PBS_JOBID | grep 'resources_used.mem\|exec_vnode' # Get memory info

FULL_TIME="$(($(date +%s)-T))"

echo "Full elapsed time in seconds: ${FULL_TIME}"

exit
