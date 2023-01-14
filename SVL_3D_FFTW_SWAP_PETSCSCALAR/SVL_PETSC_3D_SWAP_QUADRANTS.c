/***************************************************/
/* SWAP QUADRANTS (1 <--> 3, 2 <--> 4) DIAGONALLY  */
/***************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <complex.h>
# include <petsc.h>
# include "SVL_PETSC_3D_HEADER.h"

# undef __FUNCT__
# define __FUNCT__ "SVL_3D_SWAP_QUADRANTS"
PetscErrorCode SVL_3D_SWAP_QUADRANTS(PetscInt Nx, PetscInt Ny, PetscInt Nz, PetscScalar* U) {
  PetscInt i, j, k;
  PetscScalar tmp17, tmp28, tmp35, tmp46;
  PetscInt N2x = Nx/2;
  PetscInt N2y = Ny/2;
  PetscInt N2z = Nz/2;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n   STEP 2.3: SET SWAP\n");CHKERRQ(ierr);

  PetscScalar A_re[Nx * Ny * Nz];

  for (i = 0; i < Nx; i++) {
      for (j = 0; j < Ny; j++) {
          for (k = 0; k < Nz; ++k) {  
              A_re[i * Ny * Nz +  j * Nz + k] =  U[i * Ny * Nz +  j * Nz + k]; 
          }
       }
   }

  for(i = 0; i < N2x; ++i) {
     for(j = 0; j < N2y; ++j) {
          for (k = 0; k < N2z; ++k) {  
//17
	       tmp17 = A_re[i*Ny*Nz +  j*Nz + k];
	       A_re[i*Ny*Nz +  j*Nz + k] = A_re[(i + N2x)*Ny*Nz + (j + N2y)*Nz + (k + N2z)] ;
	       A_re[(i + N2x)*Ny*Nz + (j + N2y)*Nz + (k + N2z)] = tmp17;
//28
	       tmp28 = A_re[i*Ny*Nz + j*Ny + (k + N2z)];
	       A_re[i*Ny*Nz+ j*Ny + (k + N2z)] = A_re[(i + N2x)*Ny*Nz + (j + N2y)*Nz + k];
	       A_re[(i + N2x)*Ny*Nz + (j + N2y)*Nz + k] = tmp28;
//35
	       tmp35 = A_re[i*Ny*Nz +  (j + N2y)*Nz + k];
	       A_re[i*Ny*Nz +  (j + N2y)*Nz + k] = A_re[(i + N2x)*Ny*Nz + (2*j+1)*N2y + k] ;
	       A_re[(i + N2x)*Ny*Nz + (2*j+1)*N2y + k] = tmp35;
//46
	       tmp46 = A_re[i*Ny*Nz +  (2*(N2y+j)+1)*N2z + k];
	       A_re[i*Ny*Nz +  (2*(N2y+j)+1)*N2z + k] = A_re[(i + N2x)*Ny*Nz + j*Ny + k];
	       A_re[(i + N2x)*Ny*Nz + j*Ny + k] = tmp46;
           }
      }
  }

  for (i = 0; i < Nx; i++) {
       for (j = 0; j < Ny; j++) {
           for (k = 0; k < Nz; ++k) { 
               U[i*Ny*Nz +  j*Nz + k] = A_re[i*Ny*Nz +  j*Nz + k] ;
           }
       }	
  }

  PetscFunctionReturn(0);
}/* END FUNCTION SVL_SWAPQUADRANTS */
