/*************************************************************/
/*                BUILD 3D UNIT CELL ARRAY                   */
/*                       i -> Height                         */
/*                       j -> Columns                        */
/*                       k -> Rows                           */
/*************************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <petsc.h>
# include "SVL_PETSC_3D_HEADER.h"

# undef __FUNCT__
# define __FUNCT__ "SVL_3D_UNIT_CELL"
PetscErrorCode SVL_3D_UNIT_CELL(PetscInt Nx, PetscInt Ny, PetscInt Nz, PetscScalar* U, PetscReal a, PetscReal dx, PetscReal dy, PetscReal dz) {
  PetscErrorCode     ierr;
  PetscInt           i, j, k;
  PetscReal          pi = 4.0 * atan(1.0);
  PetscReal          K = 2*pi/a;  /* a = shrink : make smaller in size or amount; to cause to shrink or contract; reduce. */    
  PetscReal          array_X[Nx][Ny][Nz], array_Y[Nx][Ny][Nz], array_Z[Nx][Ny][Nz], EU[Nx][Ny][Nz];

  PetscFunctionBeginUser;

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\nSTEP 1: BUILD UNI CELL-DEVICE\n");CHKERRQ(ierr);

// Initialize the 3D Unit Cell Array
  for (i = 0; i < Nx; ++i) {
      for (j = 0; j < Ny; ++j) { 
          for (k = 0; k < Nz; ++k) {  
	    array_X[i][j][k] = (j - Ny/2) * dx;//(i * Ny * Nz +  j * Nz + k);
              array_Y[i][j][k] = (k - Nz/2) * dy;//(i * Ny * Nz +  j * Nz + k);
              array_Z[i][j][k] = (i - Nx/2) * dz;//(i * Ny * Nz +  j * Nz + k);*/
              U[i * Ny * Nz +  j * Nz + k] = 0.0 + 0.0 * I;  
          }   
      }
  }

/* Save Zero Unit Cell Array */
  //ierr = SAVE_1D_TO_3D_ARRAY_REAL("OUTPUT_UNIT_CELL_ZERO", Nx, Ny, Nz, real(U));CHKERRQ(ierr);
   ierr = SAVE_1D_TO_3D_ARRAY_COMPLEX("OUTPUT_UNIT_CELL_ZERO_REAL","OUTPUT_UNIT_CELL_ZERO_IMAG", Nx, Ny, Nz, U);CHKERRQ(ierr);
  for (i = 0; i < Nx; ++i) {
      for (j = 0; j < Ny; ++j) { 
          for (k = 0; k < Nz; ++k) {  
	    EU[i][j][k] = EU[i][j][k] + cos(K * array_X[i][j][k]);
	    EU[i][j][k] = EU[i][j][k] + cos(K * array_Y[i][j][k]);
	    EU[i][j][k] = EU[i][j][k] + cos(K * array_Z[i][j][k]);
          }   
      }
  }

  for (i = 0; i < Nx; ++i) {
      for (j = 0; j < Ny; ++j) { 
          for (k = 0; k < Nz; ++k) { 
              if (EU[i][j][k] > 0.8) { 
  	        EU[i][j][k] = 1;
              }else{
                  EU[i][j][k] = 0;
              }
          }   
      }
  }

/* Fill the Unit Cell array OUTPUT */
  for (i = 0; i < Nx; ++i) {
      for (j = 0; j < Ny; ++j) { 
          for (k = 0; k < Nz; ++k) {  
	    U[i * Ny * Nz +  j * Nz + k] = EU[i][j][k] * (1.0 + 0.0 * I);;
          }   
      }
  }

/* Save Unit Cell Array */
   //ierr = SAVE_1D_TO_3D_ARRAY_REAL("OUTPUT_UNIT_CELL_REAL", Nx, Ny, Nz, U);CHKERRQ(ierr);
   ierr = SAVE_1D_TO_3D_ARRAY_COMPLEX("OUTPUT_UNIT_CELL_REAL","OUTPUT_UNIT_CELL_IMAG", Nx, Ny, Nz, U);CHKERRQ(ierr);

   PetscFunctionReturn(0);  

} /* END FUNCTION SVL_3D_UNIT_CELL  */
