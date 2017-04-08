## source this file to use the environment installed at FNAL
export CONDA_DIR=/cvmfs/des.opensciencegrid.org/fnal/anaconda2/
export PATH=$CONDA_DIR/bin:$PATH
source activate lsstEnv
source /data/des40.a/data/marcelle/lsstEnv/bin/eups-setups.sh
setup lsst_sims
setup sims_maf
