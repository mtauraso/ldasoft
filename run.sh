#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail
set -o xtrace

# This is intended to be a sample run script to be run on a slurm node
# It launches the apptainer with the necessary load/unload and environment setup
# Currently it should only work on Hyak
# TODO: put in some bypasses so that this will "just work" on an 
# appropriately configured ubuntu/OSX setup to enable development


# TODO: make this flexible on version number in an appropriate way. 
# Should accept versions which will work, reject ones that are known bad
# given the container build.
module load ompi/4.1.3
export MPI_DIR="/sw/ompi/4.1.3/"

# TODO: be more user friendly about this. Currently specialized to me and gwastro
# Should support passing in a parameter/env var override
DATA_DIR="/gscratch/gwastro/mtauraso"

env | grep MPI_DIR

apptainer shell --env "MPI_DIR=$MPI_DIR" --bind "$MPI_DIR:$MPI_DIR:ro,$DATA_DIR:/data:rw" ldasoft_apptainer.sif 


