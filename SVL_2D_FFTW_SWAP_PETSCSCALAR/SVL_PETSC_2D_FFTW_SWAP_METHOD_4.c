# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <complex.h>
# include <petsc.h>
# include "SVL_PETSC_Header.h"

# undef __FUNCT__
# define __FUNCT__ "SVL_2D_FFTW_SWAP_METHOD_4"
PetscErrorCode SVL_2D_FFTW_SWAP_METHOD_4 (PetscInt Nx, PetscInt Ny, PetscReal* U, PetscScalar* AC) {

  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nSTEP 2: SET FFTW & SWAP \n");CHKERRQ(ierr);

/*************************************************************/
/*              METHOD 4 : FORWARD FFTW & SWAP               */
/*             1. 2D FFT on the 2D array                     */
/*************************************************************/
  ierr = SVL_2D_FFTW(Nx, Ny, U, AC);CHKERRQ(ierr);
  ierr = SAVE_1D_To_2D_ARRAY_COMPLEX("OUTPUT_FFTW_REAL", "OUTPUT_FFTW_IMAG", Nx, Ny, AC);CHKERRQ(ierr);

/*************************************************************/
/*                      BACKWARD FFTW                        */
/*             FAST FOURIER TRANSFORM ON THE WEST            */
/*************************************************************/
   PetscScalar  Inv_U[Nx * Ny];
   ierr = SVL_2D_IFFTW(Nx, Ny, AC, Inv_U);CHKERRQ(ierr);
   ierr = SAVE_1D_To_2D_ARRAY_COMPLEX("OUTPUT_INV_FFTW_REAL", "OUTPUT_INV_FFTW_IMAG" ,Nx, Ny, Inv_U);CHKERRQ(ierr);

/*************************************************************/
/*      SWAP QUADRANTS (1 <--> 3, 2 <--> 4) DIAGONALLY       */
/*************************************************************/
   ierr = SVL_2D_SWAP_QUADRANTS(Nx, Ny, AC);CHKERRQ(ierr);
   ierr = SAVE_1D_To_2D_ARRAY_COMPLEX("OUTPUT_SWAP_REAL", "OUTPUT_SWAP_IMAG" ,Nx, Ny, AC);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}/* END FUNCTION SVL SWAP QUADRANTS */
