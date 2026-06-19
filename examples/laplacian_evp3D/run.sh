#!/bin/bash
export OPENBLAS_NUM_THREADS=1
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/local/Ver0.3/lib
export AFEPACK_PATH=$HOME/local/Ver0.3/include/AFEPack
export AFEPACK_TEMPLATE_PATH=$AFEPACK_PATH/template/triangle:$AFEPACK_PATH/template/tetrahedron

 ./main $1
