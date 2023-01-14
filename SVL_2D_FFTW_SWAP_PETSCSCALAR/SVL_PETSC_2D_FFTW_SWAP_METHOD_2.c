/***************************************************/
/*              TRANSPOSE XY TO YX                 */
/***************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <complex.h>
# include <petsc.h>
# include "SVL_PETSC_Header.h"

# undef __FUNCT__
# define __FUNCT__ "SVL_2D_FFTW_SWAP_METHOD_2"
PetscErrorCode SVL_2D_FFTW_SWAP_METHOD_2(PetscInt Nx, PetscInt Ny, PetscReal* U, PetscScalar* AC) {

  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nSTEP 3: SET FFTW & SWAP \n");CHKERRQ(ierr);

/*************************************************************/
/*              METHOD 2 : FORWARD FFTW & SWAP               */
/*             1. Transpose 2D Matrix                        */ 
/*             2. 1D FFT on each Horizontal Row              */
/*             3. Transpose 2D Matrix                        */   
/*             4. 1D FFT on each Horizontal Row              */
/*************************************************************/
  PetscReal   U_array[Ny];
  PetscScalar AC_array[Ny];
  PetscInt    i, j;

/*************************************************************/
/*                   TRANSPOSE XY TO YX                      */
/*************************************************************/
  ierr = SVL_2D_TRANSPOSE_REAL(Nx, Ny, U);CHKERRQ(ierr); 

/*************************************************************/
/*   Nx x 1D Horizontal rows arrangment - FFTW along x-axis  */
/*************************************************************/
  for (i = 0; i < Nx; i++) { // Loop through the rows.
// Build a vector of size Ny  
      for (j = 0; j < Ny; j++) { // Loop through the columns. 
           U_array[j] = U[i * Ny + j] ;   // Horizontal rows arrangment
      }           
// Get FFTW for the vector         
      ierr = SVL_1D_FFTW(Nx, Ny, U_array, AC_array);CHKERRQ(ierr);  
       
// Store the results in a new vector array.
      for (j = 0; j < Ny; j++) { // Loop through the rows.
         AC[i * Ny  + j] = AC_array[j];
      }
  }

/*************************************************************/
/*                   TRANSPOSE YX TO XY                      */
/*************************************************************/
  ierr = SVL_2D_TRANSPOSE_COMPLEX(Nx, Ny, AC);CHKERRQ(ierr); 

/*************************************************************/
/*   Nx x 1D Horizontal rows arrangment - FFTW along y-axis */
/*************************************************************/
  for (i = 0; i < Nx; i++) { // Loop through the rows.
// Build a vector of size Ny  
      for (j = 0; j < Ny; j++) { // Loop through the columns.  
           U_array[j] = AC[i * Ny + j] ;   // Horizontal rows arrangment
      }           
// Get FFTW for the vector         
      ierr = SVL_1D_FFTW(Nx, Ny, U_array, AC_array);CHKERRQ(ierr);  
       
// Store the results in a new vector array.
      for (j = 0; j < Ny; j++) { // Loop through the rows.
         AC[i * Ny  + j] = AC_array[j];
      }
  }

/*************************************************************/
/*             For a total of (Nx x Ny) x 1D FFTW            */
/*************************************************************/
   ierr = SAVE_1D_To_2D_ARRAY_COMPLEX("OUTPUT_FFTW_REAL", "OUTPUT_FFTW_IMAG", Nx, Ny, AC);CHKERRQ(ierr);

/*************************************************************/
/*                      BACKWARD FFTW                        */
/*             FAST FOURIER TRANSFORM ON THE WEST            */
/*************************************************************/
   PetscScalar Inv_U[Nx * Ny];
   ierr = SVL_2D_IFFTW(Nx, Ny, AC, Inv_U);CHKERRQ(ierr);
   ierr = SAVE_1D_To_2D_ARRAY_COMPLEX("OUTPUT_INV_FFTW_REAL", "OUTPUT_INV_FFTW_IMAG" ,Nx, Ny, Inv_U);CHKERRQ(ierr);

/*************************************************************/
/*                           SWAP                            */
/*************************************************************/
/*                            2D                             */
/*                          1    2                           */ 
/*                          4    3                           */
/*      SWAP QUADRANTS (1 <--> 3, 2 <--> 4) DIAGONALLY       */
/*                          3    4                           */ 
/*                          2    1                           */
/*************************************************************/
   ierr = SVL_2D_SWAP_QUADRANTS(Nx, Ny, AC);CHKERRQ(ierr);
   ierr = SAVE_1D_To_2D_ARRAY_COMPLEX("OUTPUT_SWAP_REAL","OUTPUT_SWAP_IMAG", Nx, Ny, AC);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}/* END FUNCTION SVL SWAP QUADRANTS */
