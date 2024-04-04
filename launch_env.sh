#!/usr/bin/bash

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export MAPBOX_TOKEN=pk.eyJ1IjoiaXJhbmtvbyIsImEiOiJjbHVpa3Z4M3YwMm5rMmtsYzdxOGpid2F2In0.oX7M65zBk4WsTB2iksnxEg

if [ -z "$AGNOS_VERSION" ]; then
  export AGNOS_VERSION="9.7"
fi

export STAGING_ROOT="/data/safe_staging"
