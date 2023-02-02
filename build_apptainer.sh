#!/bin/bash
set -o errexit
set -o pipefail
set -o nounset

# Debug, switch between the oneshot and non oneshot versions
ONESHOT=0

#TODO take the build prefix from command line like install.sh does

# TODO, if you don't get a build directory from the command line `mktemp` one, 
# output it, and delete on success pulling the sif file back to wherever we are
BUILD_PREFIX=/tmp
BUILD_DIR=$BUILD_PREFIX/ldasoft_apptainer

# Note, this script only works when 
# invoked from the source directory 
# TODO fix so it uses its own directory rather than pwd
export SOURCE_DIR=$(pwd)
mkdir -p $BUILD_DIR

if [ -f "$SOURCE_DIR/ldasoft_apptainer.sif" ]; then
	rm $SOURCE_DIR/ldasoft_apptainer.sif
fi

if [[ $ONESHOT -eq 1 ]]; then
	echo "Doing a single stage build and binding UCX/OMPI from Hyak. If you are not on hyak this will fail"
	# Find and bind the directory we need for openmpi/ucx
	# Pin versions for now since so much is hardcoded in apptainer-oneshot.dev
	module load ucx/1.12.1
	module load ompi

	# Second dirname call is to knock off /bin in mpirun's path
	MPI_DIR=$(dirname $(dirname $(which mpirun)))
	UCX_DIR="/sw/ucx/1.12.1"
	apptainer build --bind "$MPI_DIR:$MPI_DIR:ro,$UCX_DIR:$UCX_DIR:ro"  $BUILD_DIR/ldasoft_apptainer.sif $SOURCE_DIR/apptainer-oneshot.def
else
	echo "Doing a two stage build. This should work everywhere"
	# If we have a two stage build, we don't care about binding libraries until runtime
	# Play through/carry on.
	apptainer build $BUILD_DIR/ldasoft_apptainer.sif $SOURCE_DIR/apptainer.def
fi


# For now copy back to source dir on success.
mv $BUILD_DIR/ldasoft_apptainer.sif $SOURCE_DIR

# Ideally on success we rm the build dir?
rm -r $BUILD_DIR
