#!/bin/bash
# create a script: setup_env.sh
#source /storage/apps/opt/spack/share/spack/setup-env.sh
spack load fftw

# Ensure library path is set
FFTW_ROOT=$(spack location -i fftw)
export LD_LIBRARY_PATH=$FFTW_ROOT/lib:$LD_LIBRARY_PATH

echo "Environment ready"
