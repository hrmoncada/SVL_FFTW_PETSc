/*************************************************************/
/*                  ADD TRIANGLE                             */
/* INSERTING AN EQUILATERAL TRIANGLE INTO A SQUARE UNIT CELL */
/*************************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <petsc.h>
# include "SVL_PETSC_Header.h"

# undef __FUNCT__
# define __FUNCT__ "SVL_2D_TRIANGLE"
PetscErrorCode SVL_2D_TRIANGLE(PetscInt Nx, PetscInt Ny, PetscReal* U) {
   PetscInt nx1, nx2, ny1, ny2;
   PetscInt nx, ny, nxa, nxb;
   PetscInt i, j;
   PetscReal f;
   PetscErrorCode ierr;

   PetscFunctionBeginUser;

   ierr = PetscPrintf(PETSC_COMM_WORLD, "\n	STEP 1.2: BUILD TRIANGLE-DEVICE\n");CHKERRQ(ierr);

  for (i = 0; i < Nx; ++i) {
      for (j = 0; j < Ny; ++j) { 
          U[i * Ny +  j] = 0.0;  // Initialize the unit cell 
      }
  }

/* Compute the position indices */
   nx1 = round(0.1*Nx);
   nx2 = round(0.9*Nx);
   ny1 = round(0.1*Ny);
   ny2 = round(0.9*Ny);

  for (ny = ny1; ny < ny2; ++ny) {
      f  =  (ny - ny1)/(double)(ny2 - ny1); // f must be declare as double
      nx  = round(f*(nx2 - nx1 + 1)); // we can not round
      nxa = round((Nx - nx)/2 + 1); // me
      //nxa = round((Nx - nx)/2); // raymond
      nxb = nxa + nx - 1; 
      for (nx = nxa-1; nx <= nxb; ++nx) {
	 U[nx * Ny + ny] = 1;
      }	
  }

  ierr = SAVE_1D_To_2D_ARRAY_REAL("OUTPUT_UNIT_CELL_REAL",Nx, Ny, U);CHKERRQ(ierr);
   
  PetscFunctionReturn(0);  
} /* END FUNCTION SVL_TRIANGLE  */
