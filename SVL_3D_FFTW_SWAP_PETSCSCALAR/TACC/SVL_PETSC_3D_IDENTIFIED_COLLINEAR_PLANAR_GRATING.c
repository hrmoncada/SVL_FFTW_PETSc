/*************************************************************/
/*           IDENTIFIED COLLINEAR PLANAR GRATING             */
/*************************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <complex.h>
# include <petsc.h>
# include "SVL_PETSC_3D_HEADER.h"

# undef __FUNCT__
# define __FUNCT__ "SVL_3D_IDENTIFIED_COLLINEAR_PLANAR_GRATING"
PetscErrorCode SVL_3D_IDENTIFIED_COLLINEAR_PLANAR_GRATING(PetscInt NM, PetscInt NN, PetscInt NP, PetscReal* KX, PetscReal* KY, PetscReal* KZ, PetscReal* KPQ){  
  PetscInt  i, j, k, n ;
  PetscInt  count, cc;
  PetscReal angle, numerator, denominator;
  PetscReal pi = 4.0 * atan(1.0);
  PetscErrorCode ierr;
  PetscFunctionBeginUser;

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n   STEP 5.2: IDENTIFIED COLLINEAR PLANAR GRATING\n");CHKERRQ(ierr);

/* Initilize the mask KPQ with 1's */ 
   for (i = 0; i < NM; i++) {
       for (j = 0; j < NN; j++) {
           for (k = 0; k < NP; k++) {
                KPQ[i * NN * NP +  j * NP + k] = 1.0;
           }
       }
   }

/*A and B are array vectors, and |A| and |B| are the norm

  1. A.B = |A|*|B|*cos(angle) ==> angle  =  ArcCos (A.B)/(|A|*|B|)

  2. Colinear vector (parallel), angle  =  0 or 180 

  3. We save one colinear vector (1's) and eliminate the rest (0's).

  4. The mask array KPQ has a size NMxNNxNP with 1's and 0's 

Since the value od the spatial harmonic is important to reconstruct the unit cell, 
Save the colinear values that has a correspondent higest spatial harmonic value on the AMNP array.
The highest harmionis are on the center,  use a spiral loop which start the loop at the center of
the tructated array AMNP.
 */
    count = 0;
    for (n = 0; n < NM * NN * NP; n++) {
        cc = 0;
        for (i = 0; i < NM; i++) {
            for (j = 0; j < NN; j++) {
                for (k = 0; k < NP; k++) {
                    if ((i * NN * NP +  j * NP + k) > k) {
    	                numerator = KX[n] * KX[i * NN * NP +  j * NP + k] + KY[n] * KY[i * NN * NP +  j * NP + k] 
                                  + KZ[n] * KZ[i * NN * NP +  j * NP + k];
                        denominator = sqrt(KX[n] * KX[n] + KY[n] * KY[n] + KZ[n] * KZ[n]) *
                                      sqrt(KX[i * NN * NP +  j * NP + k] * KX[i * NN * NP +  j * NP + k] + 
                                           KY[i * NN * NP +  j * NP + k] * KY[i * NN * NP +  j * NP + k] +
                                           KZ[i * NN * NP +  j * NP + k] * KZ[i * NN * NP +  j * NP + k]); 
                        angle  = acos(numerator/denominator);
                        if (angle == 0.0 || angle == pi) {   // Eliminate everyone that is parallel (angle 0 or 180)                           
                           KPQ[i * NN * NP +  j * NP + k] = 0.0;
                           /*ierr = PetscPrintf(PETSC_COMM_WORLD,"\n        cc = %d, n = %d count = %d, angle = %f (PX, PY, PZ) = (%d,  %d, %d), KPQ = %f\n",  cc, n, count, angle, i, j, k, KPQ[i * NN * NP +  j * NP + k]);CHKERRQ(ierr);*/
	                   cc++;  /* Count how many are collinear to a particular grating vector */
                           count++; /* Count how many are parallel */
                        }
                   }                
                }
            }  
        }
    }

   ierr = PetscPrintf(PETSC_COMM_WORLD,"\n        KPQ count = %d\n", count);CHKERRQ(ierr);

   /*mid = NN*(NN/2)+1;
   for (i = mid; i < NM*NN*NP; i++) {
       KPQ[i] = 0.0;
   }*/

   PetscFunctionReturn(0);  

} //END FUNCTION IDENTIFIED COLLINEAR PLANAR GRATING
