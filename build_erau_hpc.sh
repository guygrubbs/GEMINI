#!/bin/sh
#
# This is for Matt's machine, so he can use the older libraries he's been working with for some time.

(
module load blas/gcc/64/3.7.0
module load lapack/gcc/64/3.7.0
module load openmpi/gcc/64/1.10.3

cd objects

cmake -DMUMPS_ROOT=~/zettergmdata/lib/MUMPS_4.10.0 \
      -DSCALAPACK_ROOT=~/zettergmdata/lib/scalapack-2.0.2 \
      ..
      
cmake --build .
)
