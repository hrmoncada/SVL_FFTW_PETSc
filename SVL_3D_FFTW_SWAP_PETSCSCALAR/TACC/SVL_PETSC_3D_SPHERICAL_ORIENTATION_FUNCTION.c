/***************************************************/
/*          SVL ORIENTATION FUNCTION               */
/***************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <petsc.h>
 #include "petsctime.h"
# include "SVL_PETSC_3D_HEADER.h"

# undef __FUNCT__
# define __FUNCT__ "SVL_3D_SPHERICAL_ORIENTATION_FUNCTION"
PetscErrorCode SVL_3D_SPHERICAL_ORIENTATION_FUNCTION(PetscInt New_Nx, PetscInt New_Ny, PetscInt New_Nz, PetscReal New_dx, PetscReal New_dy, PetscReal New_dz, PetscReal* X, PetscReal* Y, PetscReal* Z, PetscReal* RSQ, PetscReal* THETA, PetscReal* VARPHI) {
  PetscErrorCode  ierr;
  PetscReal       xa, ya, za, mean_xa, mean_ya, mean_za, SHIFT_X, SHIFT_Y, SHIFT_Z;
  //PetscReal       X[New_Nx * New_Ny * New_Nz], Y[New_Nx * New_Ny * New_Nz], Z[New_Nx * New_Ny * New_Nz];
  PetscInt        i, j, k;
  PetscReal       pi = 4.0 * atan(1.0);

  PetscPrintf(PETSC_COMM_WORLD,"\n   STEP 6 : RSQ, THETA and VARPHI\n");

  PetscFunctionBeginUser;

  for (i = 0; i < New_Nx; i++) xa = xa + i * New_dy; 
  for (i = 0; i < New_Ny; i++) ya = ya + i * New_dy; 
  for (i = 0; i < New_Ny; i++) za = za + i * New_dz; 

// Average along axis
  mean_xa =  (double) xa/New_Nx;
  mean_ya =  (double) ya/New_Ny;
  mean_za =  (double) za/New_Nz;
 
// meshgrid, X, Y and Z are equal array. The arrays (X,Y,Z) are the components of the 3D cartesian coordenates plane (x,y,z)
  for (i = 0; i < New_Nx; ++i){  
      for (j = 0; j < New_Ny; ++j){
          for (k = 0; k < New_Nz; ++k){   
	       X[i * New_Ny * New_Nz + j * New_Nz + k] = i * (New_dx); // Meshgrid X, stripes along x-axis
               Y[i * New_Ny * New_Nz + j * New_Nz + k] = j * (New_dy); // Meshgrid Y, stripes along y-axis
               Z[i * New_Ny * New_Nz + j * New_Nz + k] = k * (New_dz); // Meshgrid Z, stripes along z-axis
               THETA[i * New_Ny * New_Nz + j * New_Nz + k]  = 0.0; // Initialize Rotation array
               VARPHI[i * New_Ny * New_Nz + j * New_Nz + k] = 0.0; // Initialize Rotation array
               RSQ[i * New_Ny * New_Nz + j * New_Nz + k] = 0.0;
           }
      }
  }

  ierr = SAVE_1D_TO_3D_ARRAY_REAL("OUTPUT_X", New_Nx, New_Ny, New_Nz, X);CHKERRQ(ierr);
  ierr = SAVE_1D_TO_3D_ARRAY_REAL("OUTPUT_Y", New_Nx, New_Ny, New_Nz, Y);CHKERRQ(ierr);
  ierr = SAVE_1D_TO_3D_ARRAY_REAL("OUTPUT_Z", New_Nx, New_Ny, New_Nz, Z);CHKERRQ(ierr);

/* RSQ is part of 3D Gaussian function to smooth or blur 3D Array */ 
   for (i = 0; i < New_Nx; i++){ 
       for (j = 0; j < New_Ny; j++){
          for (k = 0; k < New_Nz; ++k){   
              SHIFT_X = X[i * New_Ny * New_Nz + j * New_Nz + k] - mean_xa; // Shifted along x-axis
              SHIFT_Y = Y[i * New_Ny * New_Nz + j * New_Nz + k] - mean_ya; // Shifted along y-axis
              SHIFT_Z = Z[i * New_Ny * New_Nz + j * New_Nz + k] - mean_za; // Shifted along z-axis
              RSQ[i * New_Ny * New_Nz + j * New_Nz + k] = SHIFT_X * SHIFT_X + SHIFT_Y * SHIFT_Y + SHIFT_Z * SHIFT_Z; /* RSQ = (X-mean(X))^2 + (Y-mean(Y))^2 + (Z-mean(Z))^2 */ 
              if (RSQ[i * New_Ny * New_Nz + j * New_Nz + k] < 1) {
                  THETA[i * New_Ny * New_Nz + j * New_Nz + k]  = pi/4.0;  // Just a small part of the array is fill with differents values than zero
                  VARPHI[i * New_Ny * New_Nz + j * New_Nz + k] = pi/4.0; 
              }
          }
      }
   }

   ierr = SAVE_1D_TO_3D_ARRAY_REAL("OUTPUT_THETA", New_Nx, New_Ny, New_Nz, THETA);CHKERRQ(ierr);
   ierr = SAVE_1D_TO_3D_ARRAY_REAL("OUTPUT_VARPHI", New_Nx, New_Ny, New_Nz, VARPHI);CHKERRQ(ierr);
   ierr = SAVE_1D_TO_3D_ARRAY_REAL("OUTPUT_RSQ", New_Nx, New_Ny, New_Nz, RSQ);CHKERRQ(ierr);
   PetscFunctionReturn(0); 
} /*    END SVL ORIENTATION FUNCTION    */
