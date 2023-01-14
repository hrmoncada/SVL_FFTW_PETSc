#!/bin/sh
clear all

## petsc 3.5.4
#export PETSC_DIR=~/Desktop/PETSC/petsc-3.5.4
#export PETSC_ARCH=linux-gnu-complex #PetscScalar is Complex

# petsc 3.7.6
export PETSC_DIR=~/Desktop/PETSC/petsc-3.7.6
export PETSC_ARCH=linux-gnu-complex #PetscScalar is Complex

make SVL_PETSC
#mpirun -np 4 out_parallel
