/***************************************************/
/*          SVL_ORIENTATION_VECTOR                     */
/***************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <petsc.h>
# include "SVL_PETSC_3D_HEADER.h"

# undef __FUNCT__
# define __FUNCT__ "SVL_3D_ORIENTATION_VECTOR"
PetscErrorCode SVL_3D_ORIENTATION_VECTOR(PetscInt New_Nx, PetscInt New_Ny, PetscInt New_Nz, PetscReal* Kx, PetscReal* Ky, PetscReal* Kz, PetscReal KX, PetscReal KY, PetscReal KZ) {  
  PetscInt           i, j, k;
  //PetscErrorCode     ierr;

  PetscFunctionBeginUser;
/* Writting Grading vectors in Column array fashion, KX(:) and KY(:) */

  for (i = 0; i < New_Nx; i++) {
      for (j = 0; j < New_Ny; j++) {
          for (k = 0; k < New_Nz; k++) { 
                Kx[i * New_Nx * New_Ny + j * New_Nz + k] = KX;//KX[nk];
	      Ky[i * New_Nx * New_Ny + j * New_Nz + k] = KY;//KY[nk];  
                Kz[i * New_Nx * New_Ny + j * New_Nz + k] = KZ;//KZ[nk];   
           }
       }
   }

   /*ierr = SAVE_1D_TO_3D_ARRAY_REAL("OUTPUT_Kx", New_Nx, New_Ny, New_Nz, Kx);CHKERRQ(ierr);
   ierr = SAVE_1D_TO_3D_ARRAY_REAL("OUTPUT_Ky", New_Nx, New_Ny, New_Nz, Ky);CHKERRQ(ierr);
   ierr = SAVE_1D_TO_3D_ARRAY_REAL("OUTPUT_Kz", New_Nx, New_Ny, New_Nz, Kz);CHKERRQ(ierr);
*/
   PetscFunctionReturn(0);  
}/*   END ORIENTATION VECTOR     */

