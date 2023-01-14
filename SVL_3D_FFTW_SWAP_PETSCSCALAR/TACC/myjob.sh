#!/bin/bash

#SBATCH -J myjob                  # Job name
#SBATCH -o myjob_out_%j           # Name of stdout output file
#SBATCH -e myjob_error_%j         # Name of stderr error file
#SBATCH -p skx-normal # skx-dev   # Queue name (SKX=Skylake)
#SBATCH -N 1                      # Total # of nodes (now required)
#SBATCH -n 1                      # Total # of mpi tasks
#SBATCH -t 00:120:00              # Total run time requested <hh:mm:ss>
#SBATCH --mail-user=hrmoncadalopez@miners.utep.edu
#SBATCH --mail-type=all           # Send email at begin and end of job
#SBATCH -A TG-DMR160140           # Allocation name (req'd if more than 1)
# Other commands must follow all #SBATCH directives ...

# Load Modules: Commmented because modules are already been preloaded
#module reset
#module load petsc/3.7-complex
#module load fftw3/3.3.6
#module load matlab/2017a 
#module load python/2.7.13
module list

# Path and date 
pwd
date

# Launch MPI application...
# export PETSC_DIR=/home1/apps/intel17/impi17_0/petsc/3.7/knightslanding-complex/lib
# export PETSC_ARCH=knightslanding-complex
make SVL_PETSC

# Use ibrun instead of mpirun or mpiexec
#ibrun ./OUTPUT_PETSC -ksp_type gmres  -ksp_converged_reason > log.txt
#ibrun ./OUTPUT_PETSC -ksp_type bicg   -ksp_converged_reason > log.txt
#ibrun ./OUTPUT_PETSC -ksp_type minres -ksp_converged_reason > log.txt
#ibrun ./OUTPUT_PETSC -ksp_type cgs -ksp_converged_reason > log.txt
#ibrun ./EXECUTABLE_OUTPUT_PETSC -ksp_type lsqr -ksp_converged_reason > log.txt
#ibrun ./OUTPUT_PETSC -ksp_type cg -pc_type jacobi -ksp_converged_reason > log.txt

ibrun ./OUTPUT_PETSC -ksp_type cg -ksp_converged_reason > log.txt
