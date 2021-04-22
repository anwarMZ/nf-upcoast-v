#!/usr/bin/env bash

set -euo pipefail

SCRIPTPATH="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
BASEDIR=$(dirname $SCRIPTPATH)

singularity pull -name ${BASEDIR}/nf-upcoastv-env.sif shub://anwarMZ/nf-upcoast-v
#singularity build ${BASEDIR}/nf-upcoastv-env.sif ${BASEDIR}/environments/Singularity

echo "Built container ${BASEDIR}/nf-upcoastv-env.sif"
