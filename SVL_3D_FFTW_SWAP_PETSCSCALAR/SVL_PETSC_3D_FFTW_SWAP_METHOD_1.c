/***************************************************/
/*          FAST FOURIER TRANSFORM ON THE WEST     */
/*                       METHOD 1                  */
/***************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <complex.h>
# include <petsc.h>
# include "SVL_PETSC_3D_HEADER.h"

# undef __FUNCT__
# define __FUNCT__ "SVL_3D_FFTW_SWAP_METHOD_1"
PetscErrorCode SVL_3D_FFTW_SWAP_METHOD_1(PetscInt Nx, PetscInt Ny, PetscInt Nz, PetscScalar* U, PetscScalar* AC) {

  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nSTEP 2 : METHOD 1, SET 3D FFTW, SWAP AND IFFTW \n");CHKERRQ(ierr);

/*************************************************************/
/*                       FORWARD 3D FFTW                     */
/*************************************************************/
   ierr = SVL_3D_FFTW(Nx, Ny, Nz, U, AC);CHKERRQ(ierr);
   ierr = SAVE_1D_TO_3D_ARRAY_COMPLEX("OUTPUT_FFTW_REAL","OUTPUT_FFTW_IMAG", Nx, Ny, Nz, AC);CHKERRQ(ierr);

/*************************************************************/
/*                      BACKWARD 3D FFTW                     */
/*************************************************************/
   PetscScalar  Inv_U[Nx * Ny * Nz]; // Inv_U = INV_FFTW(U)
   ierr = SVL_3D_IFFTW(Nx, Ny, Nz, AC, Inv_U);CHKERRQ(ierr);
   ierr = SAVE_1D_TO_3D_ARRAY_COMPLEX("OUTPUT_INV_FFTW_REAL","OUTPUT_INV_FFTW_IMAG", Nx, Ny, Nz, Inv_U);CHKERRQ(ierr);

/*************************************************************/
/*                           SWAP                            */
/*************************************************************/
/*           TOP              3D           BOTTOM            */
/*         1    2                          5    6            */ 
/*         4    3                          8    7            */
/*      SWAP QUADRANTS (1 <--> 7, 2 <--> 8) DIAGONALLY       */
/*      SWAP QUADRANTS (3 <--> 5, 4 <--> 6) DIAGONALLY       */
/*                   --------------------                    */
/*                            2D                             */
/*                          1    2                           */ 
/*                          4    3                           */
/*      SWAP QUADRANTS (1 <--> 3, 2 <--> 4) DIAGONALLY       */
/*************************************************************/
   ierr = SVL_3D_SWAP_QUADRANTS(Nx, Ny, Nz, AC);CHKERRQ(ierr);
   ierr = SAVE_1D_TO_3D_ARRAY_COMPLEX("OUTPUT_SWAP_REAL","OUTPUT_SWAP_IMAG", Nx, Ny, Nz, AC);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}/* END FUNCTION SVL SWAP QUADRANTS */
