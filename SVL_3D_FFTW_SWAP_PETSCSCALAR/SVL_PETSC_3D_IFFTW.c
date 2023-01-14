/***************************************************/
/*                  FFTW                           */
/*     FAST FOURIER TRANSFORM ON THE WEST          */
/***************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <fftw3.h>
# include <complex.h>
# include <petsc.h>
# include "SVL_PETSC_3D_HEADER.h"

# undef __FUNCT__
# define __FUNCT__ "SVL_3D_IFFTW"
PetscErrorCode SVL_3D_IFFTW(PetscInt Nx, PetscInt Ny, PetscInt Nz, PetscScalar* U, PetscScalar* Inv_U) {
/* Input and Output*/ 
    PetscInt i, j, k; 
    fftw_complex    *U_in2, *U_out2; 
    PetscInt         N = Nx * Ny * Nz;
    Vec              y2, z; 
    PetscScalar      a;
    PetscBool        view = PETSC_FALSE;
    PetscErrorCode   ierr;
  
    PetscFunctionBeginUser;

    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n   STEP 2.2: BACKWARD FFTW\n");CHKERRQ(ierr);

/* FFTW backward plan */ 
    fftw_plan  b_plan;

/* Set Input and Output */ 
    U_in2  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
    U_out2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);

    ierr = VecCreateSeqWithArray(MPI_COMM_SELF, 2,(PetscInt) N,(const PetscScalar*) U_in2, &y2);CHKERRQ(ierr); 
    ierr = PetscObjectSetName((PetscObject) y2, "Frequency space vector");CHKERRQ(ierr);

    ierr = VecCreateSeqWithArray(MPI_COMM_SELF, 2,(PetscInt) N,(const PetscScalar*) U_out2, &z);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) z, "Reconstructed vector");CHKERRQ(ierr);

/*Set up plan  */
    //b_plan = fftw_plan_dft_3d(Nx, Ny, Nz, U_in2, U_out2, FFTW_BACKWARD, FFTW_ESTIMATE);
    PetscInt         dim[3];
    dim[0] = Nx; dim[1] = Ny; dim[2] = Nx;
    b_plan = fftw_plan_dft(3, dim, U_in2, U_out2, FFTW_BACKWARD, FFTW_ESTIMATE);

/* Initialize Real Space Vector (Input data): x_array -> x */
   PetscScalar  *y2_array;
   ierr = VecGetArray(y2, &y2_array);CHKERRQ(ierr);

/* x_array -> U */
   for (i = 0; i < Nx; ++i) {
       for (j = 0; j < Ny; ++j) {
           for (k = 0; k < Nz; ++k) {  
	     y2_array[i * Ny * Nz +  j * Nz + k] = U[i * Ny * Nz +  j * Nz + k]; 
           } 
       } 
   }

/* x -> x_array -> U */
   ierr = VecRestoreArray(y2, &y2_array);CHKERRQ(ierr);

   if (view) {
     ierr = VecView(y2, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
   }

/* Execute backward fftw plan, output file U_out2 */
    if (view){ierr = VecView(y2,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);}
    fftw_execute(b_plan);

/* FFTW computes an unnormalized DFT, thus z = N*x */
    a = 1.0/(PetscReal) N;
    ierr = VecScale(z,a);CHKERRQ(ierr);

/* Free spaces */
    fftw_destroy_plan(b_plan);
    fftw_free(U_in2); 
    fftw_free(U_out2);

/*************************************************************/
/*                END FUNCTION BACKWARD FFTW                 */
/*************************************************************/
/*************************************************************/
/*                           VIEW                            */
/*************************************************************/
/* view complex output*/
    PetscScalar  *ya2, *za; 
    ierr = VecGetArray(y2,&ya2);CHKERRQ(ierr);
    ierr = VecGetArray(z,&za);CHKERRQ(ierr);

// output data forward fftw  
   for (i = 0; i < Nx; i++) {
       for (j = 0; j < Ny; j++) {
           for (k = 0; k < Nz; ++k) {  
	     Inv_U[i * Ny * Nz +  j * Nz + k] = za[i * Ny * Nz +  j * Nz + k]; 
           } 
       }
   }

    ierr = VecRestoreArray(y2,&ya2);CHKERRQ(ierr);
    ierr = VecRestoreArray(z,&za);CHKERRQ(ierr);

    ierr = VecDestroy(&y2);CHKERRQ(ierr);
    ierr = VecDestroy(&z);CHKERRQ(ierr);

  PetscFunctionReturn(0); 
} /* END FUNCTION SL_IFFTW */

