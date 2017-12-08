#!/usr/bin/env bash

export INST_DIR=/global/common/cori/contrib/lsst/lsstDM
source $INST_DIR/setupStack-imsim-py2.sh
setup lsst_apps

export PYTHONPATH=~kadrlica/software/ugali/dev:$PYTHONPATH
export PYTHONPATH=~kadrlica/software/dsphsim/dev:$PYTHONPATH
