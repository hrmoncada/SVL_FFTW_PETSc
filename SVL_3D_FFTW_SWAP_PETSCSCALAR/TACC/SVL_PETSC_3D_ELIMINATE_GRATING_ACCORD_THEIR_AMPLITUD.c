/**********************************************************************************************************************************************************************/
/*                                                     ELIMINATION GRATING ACCORDING TO THEIR AMPLITUD                                                                */
/* AMNP array store the spatial harmonics and KAMNP is mask array of 0's and 1's, the 1's at KAMNP represent the harmonics values that pass the threshold.            */  
/* The threshold is used to find the harmonics that are below and above the threshold. This condition is used to build the mask KAMNP array.                          */   
/* On the KAMP array, the 1's represent Max[AMNP_{ijk}] spatial harmonics above the threshold on the AMNP array, and the 0's represent the  values below the threshold*/
/**********************************************************************************************************************************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <complex.h>
# include <petsc.h>
# include "SVL_PETSC_3D_HEADER.h"

# undef __FUNCT__
# define __FUNCT__ "SVL_3D_ELIMINATE_GRATING_ACCORD_THEIR_AMPLITUD"
PetscErrorCode SVL_3D_ELIMINATE_GRATING_ACCORD_THEIR_AMPLITUD(PetscInt NM, PetscInt NN, PetscInt NP, PetscScalar* AMNP, PetscReal* KX, PetscReal* KY, PetscReal* KZ, PetscReal* KAMNP) {
  PetscInt   i, j, k;
  PetscInt   count = 0;
  PetscReal  Threshold; // Spatial Harmonic Thershold 
  PetscErrorCode ierr;
  PetscFunctionBeginUser;

/* Find the threshold : Highest harmonics value */
   ierr = PetscPrintf(PETSC_COMM_WORLD,"\nIMPROVEMENT 1 : ELIMINATION GRATING ACCORDING TO THEIR AMPLITUD\n");CHKERRQ(ierr);

/* Initialize the mask KAMNP with 0's */ 
   for (i = 0; i < NM; i++) {
       for (j = 0; j < NN; j++) {
           for (k = 0; k < NP; k++) {
               KAMNP[i * NN * NP + j * NP + k] = 0.0;
           }
       }
   }

/* Find the Threshold */
  Threshold = cabs(AMNP[0]);   // Initialize the Thershold 

  for (i = 0; i < NM; i++) {
       for (j = 0; j < NN; j++) {
           for (k = 0; k < NP; k++) {
                if (cabs(AMNP[i * NN * NP + j * NP + k]) > Threshold) {
                   Threshold  = cabs(AMNP[i * NN * NP + j * NP + k]);  // Pick the Max[AMNP_{ijk}]
                }
           }
       }
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n        GET THRESHOLD AMPLITUD = %7.5f\n",Threshold);CHKERRQ(ierr);

/* Build the KAMNP mask array */
   count = 0;
   for (i = 0; i < NM; i++) {
       for (j = 0; j < NN; j++) {
           for (k = 0; k < NP; k++) {
               if (cabs(AMNP[i * NN * NP + j * NP + k]) > (0.02 * Threshold)){
                   KAMNP[i * NN * NP + j * NP + k] = 1.0;     // Set the value to 1's for those than pass the threshold 
                   count++;                                    // count how many are not eliminated 
               }
           }   
      }	
   }

   ierr = PetscPrintf(PETSC_COMM_WORLD,"\n        NUMBER OF HARMONICS THAN PASS THE THRESHOLD = %d\n", count);CHKERRQ(ierr);
   ierr = SAVE_1D_TO_3D_ARRAY_REAL("OUTPUT_ELIMINATION_GRATING_AMPLITUD", NM, NN, NP, KAMNP);CHKERRQ(ierr);

   PetscFunctionReturn(0);
} //END FUNCTION ELIMINATION_GRATING_ACCORD_AMPLITUD
