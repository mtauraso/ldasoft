#!/bin/bash
# All Sbatch stuff needs to be a continuous comment
#
#SBATCH --job-name=global_fit_LDC2a_dev
#SBATCH --mail-type=NONE
#SBATCH --mail-user=mtauraso
#SBATCH --account=gwastro
#SBATCH --partition=compute
#
# Note: new hyak nodes have 40 cores.
# Tyson's config has 8 tasks per node * 12 cpus per task = 96 cores
# Then something like 600 tasks are run... however 5 tasks is the minimum
#
# If we want to start 5 tasks * 12 cpus per task = 60 cores
# I'm going to reduce this to 4 cpus per task so we can get 
# 5*4 = 20 cores and run on half of the single node gwastro has
#
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --ntasks-per-node=5
#SBATCH --mem-per-cpu=4G
#
#Don't move any environment variables into slurm-land
#SBATCH --export=NONE
#
#Kill us after a time
#SBATCH --time=10:00:00


# SLURM says we get them for free... with a job ID, which is better actually
# stdout and stderr logging (Slurm docs says we'll get these files anyway)
# SBATCH --output=slurm_12.out
# SBATCH --error=slurm_12.err

#   Run with CWD set to our rundir. I think this is breaking bind
#   SBATCH --chdir=/gscratch/gwastro/mtauraso/global_fit_LDC2a_dev

# Should be container/runscript location
apptainer="/mmfs1/home/mtauraso/ldasoft/ldasoft_apptainer.sif"
globalfit="global_fit"

DATA_DIR="/gscratch/gwastro/mtauraso"

# In-container locations since we're passing them as arguments
data="/data/LDC2_sangria_training_v2.h5"
vgb="/src/ldasoft-install/data/ldc_sangria_vgb_list.dat"

# path/to/search_sources.dat (contains starting point for MBH sampler)
mbh="/src/ldasoft-install/data/"

fmin=0.0003
samples=128
samples_max=128

#Tobs=3932160
#padding=16
#outdir=/shared/ldc/sangria/prod/sangria_training_01mo

Tobs=7864320
padding=32
#ucb=/benchmarks/ldc/sangria/sangria_training_01mo/gb_catalog.cache
outdir="$DATA_DIR/global_fit_LDC2a_dev"

#Tobs=15728640
#padding=64
#ucb=/benchmarks/ldc/sangria/sangria_training_03mo/gb_catalog.cache
#outdir=/shared/ldc/sangria/prod/sangria_training_06mo

#Tobs=31457280
#padding=128
#ucb=/benchmarks/ldc/sangria/sangria_training_06mo/gb_catalog.cache
#outdir=/benchmarks/ldc/sangria/sangria_training_12mo

Tstart=0
sources=40

#Set up whatever package we need to run with
#module load gsl-2.7.1-gcc-9.4.0-ylwugg3
#module load hdf5-1.12.2-gcc-9.4.0-f5p3zoy

global_fit_cmd="global_fit \
--h5-data ${data} \
--sangria \
--fmin ${fmin} \
--chains 24 \
--start-time ${Tstart} \
--duration ${Tobs} \
--samples ${samples} \
--padding ${padding} \
--sources ${sources} \
--rundir ${outdir} \
--mbh-search-path ${mbh} \
--known-sources ${vgb} \
"
#--catalog ${ucb} \


#APPTAINER STUFF

OMPI_VERSION="4.1.3"
UCX_VERSION="1.12.1"
CUDA_VERSION="11.6.2"
module load ompi/$OMPI_VERSION
module load ucx/$UCX_VERSION
module load cuda/$CUDA_VERSION
MPI_DIR="/sw/ompi/$OMPI_VERSION"
UCX_DIR="/sw/ucx/$UCX_VERSION"
CUDA_DIR="/sw/cuda/$CUDA_VERSION"


export APPTAINERENV_PREPEND_PATH="$MPI_DIR/bin:$UCX_DIR/bin:$CUDA_DIR/bin"
# This next line is some hackpants mchacker stuff, but essentially it is the same
# as the prepend path statment above, but for LD_LIBRARY_PATH, and allows
# executables starting on the container to find MPI shared libraries we are binding from
# the outside system.
#
# See here for more discussion and where the workaround was found:
# https://github.com/apptainer/singularity/issues/5781
#
# IMHO there ought to be APPTAINERENV_PREPEND_LD_LIBRARY_PATH, but if you look at
# implementation of APPTAINERENV_PREPEND_PATH at time of writing the bypass on
# apptainer's end is specific only to PATH, and is somewhat awkward to generalize.
export APPTAINERENV_LD_LIBRARY_PATH="$MPI_DIR/lib:$UCX_DIR/lib:$CUDA_DIR/lib64:\$LD_LIBRARY_PATH"



#global_fit_cmd="env"

#cmd="apptainer run --bind $MPI_DIR:$MPI_DIR:ro,$UCX_DIR:$UCX_DIR:ro,$CUDA_DIR:$CUDA_DIR:ro,$DATA_DIR:/data:rw ${apptainer} ${global_fit_cmd}"

cmd="apptainer run --bind $MPI_DIR:$MPI_DIR:ro,$UCX_DIR:$UCX_DIR:ro,$CUDA_DIR:$CUDA_DIR:ro,$DATA_DIR:/data:rw ${apptainer} mpirun -np $SLURM_NTASKS ${global_fit_cmd}"

echo $cmd
# Tyson's config had this at 18 w/ 12 cores available to each process (cores x 1.5). 
# Linearly we should have 6 for 4 cores
export OMP_NUM_THREADS=6
#export UCX_MEMTYPE_CACHE=n
export OMPI_MCA_mpi_cuda_support=0

$cmd

#mpirun -np $SLURM_NTASKS -x UCX_POSIX_USE_PROC_LINK=n -x UCX_TLS=tcp,self $cmd


#(time mpirun -np $SLURM_NTASKS $cmd) 1> run.out 2>&1
