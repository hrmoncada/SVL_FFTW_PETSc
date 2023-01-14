/***************************************************/
/*          SVL_ORIENTATION_VECTOR                     */
/***************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <petsc.h>
# include "SVL_PETSC_3D_HEADER.h"

# undef __FUNCT__
# define __FUNCT__ "SVL_3D_RHS"
PetscErrorCode SVL_3D_RHS(PetscInt New_Nx, PetscInt New_Ny, PetscInt New_Nz, PetscReal* Kx, PetscReal* Ky, PetscReal* Kz, PetscReal* f) {  
  PetscInt i, j, k;
  PetscInt M = New_Nx * New_Ny * New_Nz;

  PetscFunctionBeginUser;

  for (i = 0; i < New_Nx; i++) {
      for (j = 0; j < New_Ny; j++) {
          for (k = 0; k < New_Nz; k++) { 
              f[i * New_Nx * New_Ny + j * New_Nz + k]       = Kx[i * New_Nx * New_Ny + j * New_Nz + k];
	    f[i * New_Nx * New_Ny + j * New_Nz + k + M]   = Ky[i * New_Nx * New_Ny + j * New_Nz + k];  
              f[i * New_Nx * New_Ny + j * New_Nz + k + 2*M] = Kz[i * New_Nx * New_Ny + j * New_Nz + k];   
           }
       }
   }

   PetscFunctionReturn(0);  
}/*   END GRADING VECTOR     */

