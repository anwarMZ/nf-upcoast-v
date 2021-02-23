#!/bin/bash
set -eo pipefail
#export PATH=/opt/conda/bin:$PATH
export PATH=/Users/au572806/opt/anaconda3/bin:$PATH

mkdir conda_cache_dir

# first NF run will create the conda env in the cache dir
echo run pipeline with conda --cache to create cache.. >> artifacts/test_artifact.log
export REPO=$PWD
echo REPO=$REPO >> artifacts/test_artifact.log
cd ..
echo PWD=$PWD >> $REPO/artifacts/test_artifact.log
NXF_VER=20.03.0-edge nextflow run $REPO \
       -profile conda \
       --cache $REPO/conda_cache_dir \
       --directory $REPO/.github/data/fastqs/ \
       --illumina \
       --prefix test
cp .nextflow.log $REPO/artifacts/cache_creation.conda.profile.nextflow.log

cat .nextflow.log | grep 'Conda create complete env=/Users/au572806/GitHub/nf-upcoast-v/environments/environment.yml path=/Users/au572806/GitHub/nf-upcoast-v/nf-upcoast-v/conda_cache_dir/nf-upcoastv-env' \
    && echo "Conda env created in cache dir" >> $REPO/artifacts/test_artifact.log \
	|| bash -c "echo test failed\: Conda environment not created as expected >> $REPO/artifacts/test_artifact.log && exit 1"

rm -rf results && rm -rf work && rm -rf .nextflow*
# second NF run will use the conda env created in the previous run
echo re-run pipeline with conda --cache.. >> $REPO/artifacts/test_artifact.log
NXF_VER=20.03.0-edge nextflow run $REPO \
       -profile conda \
       --cache $REPO/conda_cache_dir \
       --directory $REPO/.github/data/fastqs/ \
       --illumina \
       --prefix test
cp .nextflow.log $REPO/artifacts/cache_use.conda.profile.nextflow.log

cat .nextflow.log | grep 'Conda found local env for environment=/Users/au572806/GitHub/nf-upcoast-v/environments/environment.yml; path=/Users/au572806/GitHub/nf-upcoast-v/nf-upcoast-v/conda_cache_dir/nf-upcoastv-env' \
    && echo "Conda env found in cache dir" >> $REPO/artifacts/test_artifact.log \
	|| bash -c "echo test failed\: Conda environment not not found in cache dir >> $REPO/artifacts/test_artifact.log && exit 1"

# clean-up for following tests
rm -rf results && rm -rf work && rm -rf .nextflow*
