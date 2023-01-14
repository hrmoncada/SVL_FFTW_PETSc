/*****************************************************************************************************************/
/*                                            IMPLEMENT_IMPROVEMENTS                                             */
/*****************************************************************************************************************/
/* This file must also contain the IMPROMENT 3 which is not include on the file yet                              */
/* Array AMNP[NK]  TRUCATE THE SPATIAL HARMONICS                                                                 */
/* Array KAMNP[NK] IMPROVEMENT 1: ELIMINATION GRATING ACCORDING TO THEIR AMPLITUD                                */
/* Array KPQ[NK]   IMPROVEMENT 2: IDENTIFIED COLLINEAR PLANAR GRATING                                            */
/* Array KC[NK]    IMPROVEMENT 3: IMPLEMENT IMPROVEMENTS 1 + 2 - Mask collinear planar gratings KC = KAMNP * KPQ */
/* Array KPQA[NK]  IMPROVEMENT 4:(No implemented yet) PHASE FUNCTION DIFFER BY ONLY A CONSTANT ANGLE ALPHA_i     */
/* KX[NK], KY[NK], KZ[NK]   Store Grading vector values                                                          */
/* PX[NK], PY[NK], PZ[NK]   Store X-Position, Y-Position, Z-Position Highest Spatial Harmonics                   */
/*****************************************************************************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <complex.h>
# include <petsc.h>
# include "SVL_PETSC_3D_HEADER.h"

# undef __FUNCT__
# define __FUNCT__ "SVL_3D_IMPLEMENT_IMPROVEMENTS"
PetscErrorCode SVL_3D_IMPLEMENT_IMPROVEMENTS(PetscInt NM, PetscInt NN, PetscInt NP, PetscInt* NK, PetscReal* KX, PetscReal* KY, PetscReal* KZ, PetscScalar* AMNP) {
  PetscInt      i, j, k;
  PetscInt      count_KAMNP = 0;                   // Counter IMPROVEMENT 1
  PetscReal     KAMNP[NM * NN * NP];               // ARRAY IMPROVEMENT 1: Mask - Highest FFTW Amplitudes
// PetscReal    KPQ[NM * NN * NP];                 // ARRAY IMPROVEMENT 2: Mask - Collinear Planar Gratings 
// PetscReal    KC[NM * NN * NP];                  // ARRAY IMPROVEMENT 3: Mask - IMPLEMENT IMPROVEMENTS 1 + 2  (KC = KAMNP * KPQ)
// PetscInt     PX[NM * NN * NP], PY[NM * NN * NP], PZ[NM * NN * NP];  // Store x-Position, y-Position, z-Position Highest Spatial Harmonics 
  PetscErrorCode ierr;
 
  PetscFunctionBeginUser;  

/* (*NK) parameters is passing by reference. NK << NM*NN*NP, On each improvement NK get Update to a smaller value.
   &Threshold and &NK parameters are passing by reference, mean is passed as pointers (*NK).*/

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nSTEP 5 : IMPLEMENT IMPROVEMENTS\n");CHKERRQ(ierr);
/**********************************************************************************************/
/* IMPROVEMENT 1: ELIMINATION GRATING ACCORDING TO THEIR AMPLITUD         NC = AMNP*KAMNP     */
/* IMPROVEMENT 2: IDENTIFIED COLLINEAR PLANAR GRATING                     NC = AMNP*KPQ       */
/* IMPROVEMENT 3: IMPLEMENT IMPROVEMENTS 1 + 2                            NC = AMNP*KAMNP*KPQ */
/* IMPROVEMENT 4: PHASE FUNCTION DIFFER BY ONLY A CONSTANT ANGLE ALPHA_i (No implemented yet) */
/**********************************************************************************************/
   int operator = 1;
   switch(operator) {
       case 1:
           ierr = SVL_3D_ELIMINATE_GRATING_ACCORD_THEIR_AMPLITUD(NM, NN, NP, AMNP, KX, KY, KZ, KAMNP);CHKERRQ(ierr); 

	   for (i = 0; i < NM; i++) {
	       for (j = 0; j < NN; j++) {
		 for (k = 0; k < NP; k++) {
		     if (KAMNP[i * NN * NP + j * NP + k] != 0.0) {
		         KX[count_KAMNP] = KAMNP[i * NN * NP + j * NP + k] * KX[i * NN * NP + j * NP + k];
		         KY[count_KAMNP] = KAMNP[i * NN * NP + j * NP + k] * KY[i * NN * NP + j * NP + k];
		         KZ[count_KAMNP] = KAMNP[i * NN * NP + j * NP + k] * KZ[i * NN * NP + j * NP + k];
		       AMNP[count_KAMNP] = KAMNP[i * NN * NP + j * NP + k] * AMNP[i * NN * NP + j * NP + k] ; // NEW AMNP ARRAY OF SIZE NK << NM*NN*NP
	               //ierr = PetscPrintf(PETSC_COMM_WORLD,"\ncount = %d, (PX, PY, PZ) = (%d,  %d, %d), (KX, KY, KZ) = (%f,  %f, %f), AMNP[%d] = %f\n", count_KAMNP, i, j, k, KX[count_KAMNP], KY[count_KAMNP], KZ[count_KAMNP], count_KAMNP, AMNP[count_KAMNP]);CHKERRQ(ierr);
		         count_KAMNP++;
		     }
		  }   
	         }
	     }
/* count_KAMNP count, and NK is the total number of Spatial Harmonics than pass the threshold and don't fall below that threshold.
   The Update NK will fix AMNP new size, since the array got smaller NK <<< NM * NN * NP
   In short, We eliminate all the gratings (Kx, KY, KZ) from the expansion with FFT amplitudes that fall below the threshold.*/
             *NK =  count_KAMNP; 
	     ierr = PetscPrintf(PETSC_COMM_WORLD,"\n(Number of iteration, array size) =( %d, %d) \n", count_KAMNP, *NK);CHKERRQ(ierr);      
             break;
       case 2:
      /*    ierr = PetscPrintf(PETSC_COMM_WORLD,"\nIMPROVEMENT 2 : IDENTIFIED COLLINEAR PLANAR GRATINGS \n");CHKERRQ(ierr);
	    ierr = SVL_3D_IDENTIFIED_COLLINEAR_PLANAR_GRATING(NM, NN, NP, KX, KY, KZ, KPQ); CHKERRQ(ierr);
	    ierr = SAVE_1D_TO_3D_ARRAY_REAL("OUTPUT_COLLINEAR_PLANAR_GRATING", NM, NN, NP, KPQ);CHKERRQ(ierr);
            *NK  =  count_KPQ; // This will fix AMNP new size, since the array got smaller NK <<< NM * NN * NP
	    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n(Number of iteration, array size) =( %d, %d) \n", count_KPQ, *NK);CHKERRQ(ierr);   
	*/
           break;
       case 3:
        /*  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nIMPROVEMENT 3: IMPLEMENT IMPROVEMENTS 1 + 2\n");CHKERRQ(ierr);
	    ierr = SVL_3D_IMPROVEMENTS_12(NM, NN, NP, NK, AMNP, KX, KY, KZ, KAMNP, KPQ, KC);CHKERRQ(ierr); 
	    ierr = SAVE_1D_TO_3D_ARRAY_REAL("OUTPUT_IMPLEMENT_IMPROVEMENTS_1_AND_2", NM, NN, NP, KC);CHKERRQ(ierr);
            *NK  =  count_KPQ; // This will fix AMNP new size, since the array got smaller NK <<< NM * NN * NP
	    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n(Number of iteration, array size) =( %d, %d) \n", count_KC, *NK);CHKERRQ(ierr); 
	*/
           break;
       case 4:
        /*  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nMPROVEMENT 4: PHASE FUNCTION DIFFER BY ONLY A CONSTANT ANGLE ALPHA_i\n");CHKERRQ(ierr);
	    ierr = SVL_3D_PHASE_FUNCTION_DIFFER_CONSTANT_ANGLE(NM, NN, NP, NK, AMNP, KX, KY, KZ, KAMNP, KPQ, KC);CHKERRQ(ierr); 
	    ierr = SAVE_1D_TO_3D_ARRAY_REAL("OUTPUT_IMPLEMENT_IMPROVEMENTS_1_AND_2", NM, NN, NP, KC);CHKERRQ(ierr);
            *NK  =  count_KPQ; // This will fix AMNP new size, since the array got smaller NK <<< NM * NN * NP
	    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n(Number of iteration, array size) =( %d, %d) \n", count_KC, *NK);CHKERRQ(ierr); 
	*/
           break;
       default:
           ierr = PetscPrintf(PETSC_COMM_WORLD,"\nTRY AGAIN!!!, YOU DID NOT IMPLEMENTED ANY IMPROVEMENT \n");CHKERRQ(ierr);
           break;
   }

  PetscFunctionReturn(0);
} //END FUNCTION ELIMINATION_GRATING_ACCORD_AMPLITUD
