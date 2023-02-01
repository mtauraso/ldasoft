#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail
#set -o xtrace

# This is intended to be a sample run script to be run on a slurm node
# It launches the apptainer with the necessary load/unload and environment setup
# Currently it should only work on Hyak/similar slurm systems
#
# TODO: put in some bypasses so that this will "just work" on an 
# appropriately configured ubuntu/OSX setup to enable development


# TODO: make this flexible on version number in an appropriate way. 
# Should accept versions which will work, reject ones that are known bad
# given the container build.
# Right now we just take whatever the slurm env gives us and run with it
module load ompi

# Second dirname call is to knock off /bin in mpirun's path
MPI_DIR=$(dirname $(dirname $(which mpirun)))


# XCXC: This maybe shouldn't be a bind at all. This might be something we pawn off on the user
#       To make sure the directories on the host they want the containerized LDASOFT to interact
#       with are actually availble.... but the best way to do this is probably the way apptainer already does it
#       host dirs are just available if there aren't conflicts in dir space.
#
# TODO: be more user friendly about this. Currently specialized to me and gwastro
# Should support passing in a parameter/env var override
DATA_DIR="/gscratch/gwastro/mtauraso"


# IS this necessary? Intent was to help the container know it was
# loaded in bind mode, but there may be a better way.
export APPTAINERENV_MPI_DIR="$MPI_DIR"

export APPTAINERENV_PREPEND_PATH="$MPI_DIR/bin"
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
export APPTAINERENV_LD_LIBRARY_PATH="$MPI_DIR/lib:\$LD_LIBRARY_PATH"

apptainer run --bind "$MPI_DIR:$MPI_DIR:ro,$DATA_DIR:/data:rw" ldasoft_apptainer.sif $@


