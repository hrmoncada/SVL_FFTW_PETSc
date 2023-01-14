/*************************************************************/
/*                   IMPLEMENT_IMPROVEMENTS                  */
/*************************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <complex.h>
# include <petsc.h>
# include "SVL_PETSC_3D_HEADER.h"

# undef __FUNCT__
# define __FUNCT__ "SVL_3D_IMPROVEMENTS"
PetscErrorCode SVL_3D_IMPROVEMENTS(PetscInt NM, PetscInt NN, PetscInt NP, PetscInt* NK, PetscScalar* AMNP, PetscReal* KX, PetscReal* KY, PetscReal* KZ, PetscReal* KAMNP, PetscReal* KPQ, PetscReal* KC) {

  PetscInt  i, j, k;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

/* *NK parameters is passing by reference*/
/* Array AMNP  TRUCATE THE SPATIAL HARMONICS  */
/* Array KAMNP IMPROVEMENT 1: ELIMINATION GRATING ACCORDING TO THEIR AMPLITUD*/
/* Array KPQ   IMPROVEMENT 2: IDENTIFIED COLLINEAR PLANAR GRATING */
/* Array KC    IMPLEMENT IMPROVEMENTS 1 & 2 */
/* This file must also contain the IMPROMENT 3 which is not include on the file yet*/
/* Array K???  IMPROVEMENT 3: PHASE FUNCTION DIFFER BY ONLY A CONSTANT ANGLE ALPHA_i >>> STILL NO IMPLEMENTED */
/* STILL Need to implement the third improvement to be combine with 1st and 2nd  */

/* 1: IMPROMENT 1 (ELIMINATIONING GRATINGS ACCORDING TO THEIR AMPLITUD - KAMNP)*/
/* 2: IMPROMENT 2 (IDENTIFIED COLLINEAR PLANAR GRATING - KPQ)*/
/* 3: IMPROMENTS 1 AND 2  (KAMNP * KPQ) */
/* 4: NO IMPROMENTS */

   PetscInt improvement = 1;

   switch (improvement) {
	case 1:
	     ierr = PetscPrintf(PETSC_COMM_WORLD,"\n   STEP 5.4: IMPLEMENT IMPROVEMENT 1\n");CHKERRQ(ierr);
	     PetscInt  count_KAMNP = 0;
               for (i = 0; i < NM; i++) {
                   for (j = 0; j < NN; j++) {
                       for (k = 0; k < NP; k++) {
                           if (KAMNP[i * NN * NP +  j * NP + k] != 0.0) {
                              KX[count_KAMNP] = KAMNP[i * NN * NP +  j * NP + k] * KX[i * NN * NP +  j * NP + k];
                              KY[count_KAMNP] = KAMNP[i * NN * NP +  j * NP + k] * KY[i * NN * NP +  j * NP + k];
                              KZ[count_KAMNP] = KAMNP[i * NN * NP +  j * NP + k] * KZ[i * NN * NP +  j * NP + k];
                            AMNP[count_KAMNP] = KAMNP[i * NN * NP +  j * NP + k] * AMNP[i * NN * NP +  j * NP + k] ; // NEW AMNP ARRAY OF SIZE NK << NM*NN*NP
//ierr = PetscPrintf(PETSC_COMM_WORLD,"\ncount = %d, (PX, PY, PZ) = (%d,  %d, %d), (KX, KY, KZ) = (%f,  %f, %f), AMNP[%d] = %f\n", count_KAMNP, i, j, k, KX[count_KAMNP], KY[count_KAMNP], KZ[count_KAMNP], count_KAMNP, AMNP[count_KAMNP]);CHKERRQ(ierr);
                            count_KAMNP++;
                          }
                       }   
                    }
                }
                ierr = PetscPrintf(PETSC_COMM_WORLD,"\n        Number of iteration = %d\n", count_KAMNP);CHKERRQ(ierr);

/* NK get update, *NK parameters is passing by reference, and NK << NM*NN : NK count the number of Spatial Harmonics than pass the threshold and do not fall below that threshold */
               *NK =  count_KAMNP;
	     break;

	case 2:
	     ierr = PetscPrintf(PETSC_COMM_WORLD,"\n   STEP 5.4: IMPLEMENT IMPROVEMENT 2\n");CHKERRQ(ierr);
               PetscInt  count_KPQ = 0;
               for (i = 0; i < NM; i++) {
                   for (j = 0; j < NN; j++) {
                       for (k = 0; k < NP; k++) {
                           if (KPQ[i * NN * NP +  j * NP + k] != 0.0) {
                              KX[count_KPQ]  = KX[i * NN * NP +  j * NP + k] * KPQ[i * NN * NP +  j * NP + k];
                              KY[count_KPQ]  = KY[i * NN * NP +  j * NP + k] * KPQ[i * NN * NP +  j * NP + k];
                              KZ[count_KPQ]  = KZ[i * NN * NP +  j * NP + k] * KPQ[i * NN * NP +  j * NP + k];
                              AMNP[count_KPQ]  = AMNP[i * NN * NP +  j * NP + k] * KPQ[i * NN * NP +  j * NP + k]; // NEW AMNP ARRAY OF SIZE NK << NM*NN*NP
                        //ierr = PetscPrintf(PETSC_COMM_WORLD,"\ncount = %d, (PX, PY, PZ) = (%d,  %d, %d), (KX, KY, KZ) = (%f,  %f, %f)\n", count_KPQ, i, j, k, KX[count_KPQ], KY[count_KPQ], KZ[count_KPQ]);CHKERRQ(ierr);
                             count_KPQ++;
                          }
                      }   
                   }
                }
                ierr = PetscPrintf(PETSC_COMM_WORLD,"\n        Number of iteration = %d\n", count_KPQ);CHKERRQ(ierr);

 /* NK get update, *NK parameters is passing by reference, and  NK <<< NM*NN : NK count the number of grating vector after improvement 2 */
	     *NK =  count_KPQ;
               break;

	case 3:
	     ierr = PetscPrintf(PETSC_COMM_WORLD,"\n   STEP 5.4: IMPLEMENT IMPROVEMENTS 1 AND 2\n");CHKERRQ(ierr);
	     PetscInt  count_KC = 0;
	     for (i = 0; i < NM; i++) {
	         for (j = 0; j < NN; j++) {
                        for (k = 0; k < NP; k++) {
	                 KC[i * NN * NP +  j * NP + k]  = KAMNP[i * NN * NP +  j * NP + k] * KPQ[i * NN * NP +  j * NP + k]; /* COMBINE IMPROMENTS 1 and 2 */
                           if (KC[i * NN * NP +  j * NP + k] == 1) {
                              KX[count_KC]  = KX[i * NN * NP +  j * NP + k] * KC[i * NN * NP +  j * NP + k];
                              KY[count_KC]  = KY[i * NN * NP +  j * NP + k] * KC[i * NN * NP +  j * NP + k];
                              KZ[count_KC]  = KZ[i * NN * NP +  j * NP + k] * KC[i * NN * NP +  j * NP + k];
                              AMNP[count_KC]  = AMNP[i * NN * NP +  j * NP + k] * KC[i * NN * NP +  j * NP + k]; // NEW AMNP ARRAY OF SIZE NK << NM*NN*NP
                              //ierr = PetscPrintf(PETSC_COMM_WORLD,"\ncount = %d, (PX, PY, PZ) = (%d,  %d, %d), (KX, KY, KZ) = (%f,  %f, %f)\n", count_KC, i, j, k, KX[count_KC], KY[count_KC], KZ[count_KC]);CHKERRQ(ierr);
                             count_KC++;
                          }
                       }
                    }
                }
                ierr = PetscPrintf(PETSC_COMM_WORLD,"\n        Number of iteration = %d\n", count_KC);CHKERRQ(ierr);

/* NK get update, *NK parameters is passing by reference, and  NK < NM*NN : NK count the number grating vector remain after combine improments 1 & 2 */
	     *NK =  count_KC; 
               break;

	default:
	     ierr = PetscPrintf(PETSC_COMM_WORLD,"\n   STEP 5.4: NO IMPROVEMENTS ARE IMPLEMENTED \n");CHKERRQ(ierr);

   } // end switch

   PetscFunctionReturn(0);
} //END FUNCTION IMPROVEMENTS
