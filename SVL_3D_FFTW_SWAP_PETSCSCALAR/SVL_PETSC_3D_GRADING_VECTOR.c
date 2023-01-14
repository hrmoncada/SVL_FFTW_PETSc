/***************************************************/
/*          SVL_GRADING_VECTOR                     */
/***************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <petsc.h>
# include "SVL_PETSC_3D_HEADER.h"

# undef __FUNCT__
# define __FUNCT__ "SVL_3D_GRADING_VECTOR"
PetscErrorCode SVL_3D_GRADING_VECTOR(PetscInt NM, PetscInt NN, PetscInt NP, PetscReal  Sx, PetscReal Sy, PetscReal Sz, PetscReal* KX, PetscReal* KY, PetscReal* KZ){
  PetscInt i, j, k;  
  PetscReal  KX_3D_plane[NM * NN * NP], KY_3D_plane[NM * NN * NP], KZ_3D_plane[NM * NN * NP];
  PetscReal pi = 4.0 * atan(1.0);   
  PetscErrorCode ierr;

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n   STEP 4: GRADING VECTOR (KX, KY, KZ) - WAVE VECTOR IN CARTESIAN COORDINATES\n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n   	  Writing KX(:), KY(:) and KZ(:) in a column array fashion\n");CHKERRQ(ierr);
  
  PetscFunctionBeginUser;

/* Writting Grading vectors in Column array fashion, KX(:) and KY(:) */

/* 3D planes of grading vector */
   for (i = 0; i < NM; i++) {
       for (j = 0; j < NN; j++) {
           for (k = 0; k < NP; k++) {
                  KX_3D_plane[i * NN * NP +  j * NP + k]  =  (i - NM/2); 
	        KY_3D_plane[i * NN * NP +  j * NP + k]  =  (j - NN/2);
	        KZ_3D_plane[i * NN * NP +  j * NP + k]  =  (k - NP/2);
            }             
        } 
    }

   ierr = SAVE_1D_TO_3D_ARRAY_REAL("OUTPUT_KX_3D_plane", NM, NN, NP, KX_3D_plane);CHKERRQ(ierr);
   ierr = SAVE_1D_TO_3D_ARRAY_REAL("OUTPUT_KY_3D_plane", NM, NN, NP, KY_3D_plane);CHKERRQ(ierr);
   ierr = SAVE_1D_TO_3D_ARRAY_REAL("OUTPUT_KZ_3D_plane", NM, NN, NP, KZ_3D_plane);CHKERRQ(ierr);

/* 3D grading vector */
   for (i = 0; i < NM; i++) {
       for (j = 0; j < NN; j++) {
           for (k = 0; k < NP; k++) {
                  KX[i * NN * NP +  j * NP + k]  = (2*pi/Sx) * (i - NM/2); 
	        KY[i * NN * NP +  j * NP + k]  = (2*pi/Sy) * (j - NN/2);
	        KZ[i * NN * NP +  j * NP + k]  = (2*pi/Sz) * (k - NP/2);
            }             
        } 
    }

   ierr = SAVE_1D_TO_3D_ARRAY_REAL("OUTPUT_KX", NM, NN, NP, KX);CHKERRQ(ierr);
   ierr = SAVE_1D_TO_3D_ARRAY_REAL("OUTPUT_KY", NM, NN, NP, KY);CHKERRQ(ierr);
   ierr = SAVE_1D_TO_3D_ARRAY_REAL("OUTPUT_KZ", NM, NN, NP, KZ);CHKERRQ(ierr);

     PetscFunctionReturn(0);  
}/*   END GRADING VECTOR     */

