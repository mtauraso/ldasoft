#!/bin/bash
set -o errexit
set -o pipefail
set -o nounset


#TODO  xcxc take the build prefix from command line like install.sh does

# TODO, if you don't get a build directory from the command line `mktemp` one, 
# output it, and delete on success pulling the sif file back to wherever we are
BUILD_PREFIX=/tmp

BUILD_DIR=$BUILD_PREFIX/ldasoft_apptainer

# Note, this script only works when 
# invoked from the source directory 
# TODO fix so it uses its own directory rather than pwd
export SOURCE_DIR=$(pwd)

mkdir -p $BUILD_DIR
# TODO find and bind the directory we need for openmpi?


if [ -f "$SOURCE_DIR/ldasoft_apptainer.sif" ]; then
	rm $SOURCE_DIR/ldasoft_apptainer.sif
fi

apptainer build $BUILD_DIR/ldasoft_apptainer.sif $SOURCE_DIR/apptainer.def

# For now copy back to source dir on success.
mv $BUILD_DIR/ldasoft_apptainer.sif $SOURCE_DIR


# Ideally on success we rm the build dir?
rm -r $BUILD_DIR
