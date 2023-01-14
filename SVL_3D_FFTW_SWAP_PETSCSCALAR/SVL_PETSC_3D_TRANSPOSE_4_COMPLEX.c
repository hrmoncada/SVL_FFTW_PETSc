/***************************************************/
/*           TRANSPOSE 2 FROM XZY TO XYZ           */
/***************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <complex.h>
# include <petsc.h>
# include "SVL_PETSC_3D_HEADER.h"

# undef __FUNCT__
# define __FUNCT__ "SVL_3D_TRANSPOSE_4_COMPLEX"
PetscErrorCode SVL_3D_TRANSPOSE_4_COMPLEX(PetscInt Nx, PetscInt Ny, PetscInt Nz, PetscScalar* U) {
  PetscInt i, j, k;  
  PetscScalar T[Nx * Ny * Nz];// Transpose matrix
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n	STEP 3.1 : TRANSPOSE 4\n");CHKERRQ(ierr);

  for (i = 0; i < Nx; i++) { // Loop through the height.
      for (j = 0; j < Ny; j++) { // Loop through the rows.   
          for (k = 0; k < Nz; k++) { // Loop through the columns.
              T[i + j * Nx * Nz + k * Nx] = U[i * Ny * Nz + j * Nz + k]; //also work
              //T[i + j * Nx + k * Nx * Ny] = U[i * Ny * Nz + j * Nz + k]; 
          }
      }
  }

  for (i = 0; i < Nx; i++) {
     for (j = 0; j < Ny; j++) {
         for (k = 0; k < Nz; k++) {
             U[i * Ny * Nz + j * Nz + k] = T[i * Ny * Nz + j * Nz + k];
         }
     }	
  }


  PetscFunctionReturn(0);
}/* END FUNCTION SVL SWAP QUADRANTS */


