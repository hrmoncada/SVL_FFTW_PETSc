/***************************************************/
/*                    FFTW-SWAP                    */
/***************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <complex.h>
# include <petsc.h>
# include "SVL_PETSC_3D_HEADER.h"

# undef __FUNCT__
# define __FUNCT__ "SVL_3D_FFTW_SWAP_METHOD_2"
PetscErrorCode SVL_3D_FFTW_SWAP_METHOD_2(PetscInt Nx, PetscInt Ny, PetscInt Nz, PetscScalar* U, PetscScalar* AC) {
  PetscInt i, j, k;  
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nSTEP 3: SET FFTW & SWAP - METHOD 1\n");CHKERRQ(ierr);

/******************************************************************************/
/*                         METHOD 2 : FORWARD FFTW & SWAP                     */
/* for (i = 0; i < Nx; i++) { // Loop through the height (X AXIS)             */
/*   for (j = 0; j < Ny; j++) { // Loop through the rows (Y AXIS).            */ 
/*     for (k = 0; k < Nz; k++) { // Loop through the columns (Z AXIS)        */
/*        U[i * Ny * Nz + j * Nz + k] = {{{0 1 2},{3 4 5},{6 7 8}},           */  
/*                                       {{9 10 11},{12 13 14;},{15 16 17}},  */
/*                                       {{18 19 20},{21 22 23},{24 25 26}}}; */
/*                                                                              */
/* 1. Transpose 1: ALONG Y                                                    */
/*        T[i * Ny * Nz + j + k * Ny] = U[i * Ny * Nz + j * Nz + k];          */ 
/*        T[i * Ny * Nz + j + k * Ny] = {{{0 3 6},{1 4 7},{2 5 8}},           */
/*                                       {{9 12 15},{10 13 16},{11 14 17}},   */
/*                                       {{18 21 24},{19 22 25},{20 23 26}}}; */
/*                                                                            */ 
/*    1D FFT ALONG Y                                                          */
/*                                                                            */
/* 2. Transpose 2: ALONG Z                                                    */
/*        T[i * Ny * Nz + j + k * Ny] = U[i * Ny * Nz + j * Nz + k];          */ 
/*        T[i * Ny * Nz + j + k * Ny] = {{{0 1 2},{3 4 5},{6 7 8}},           */  
/*                                       {{9 10 11},{12 13 14;},{15 16 17}},  */
/*                                       {{18 19 20},{21 22 23},{24 25 26}}}; */
/*    1D FFT ALONG Z                                                          */
/*                                                                            */ 
/* 3. Transpose 3: ALONG X                                                    */
/*       T[i + j * Nx * Nz + k * Nx] = U[i * Ny * Nz + j * Nz + k];           */
/*       T[i + j * Nx * Nz + k * Nx] = {{{0 9 18},{1 10 19},{2 11 20}},       */
/*                                      {{3 12 21},{4 13 22},{5 14 23}},      */
/*                                      {{6 15 24},{7 16 25},{8 17 26}}};     */
/*                                                                            */
/* OR   T[i + j * Nx + k * Nx * Ny] = U[i * Ny * Nz + j * Nz + k];            */ 
/*      T[i + j * Nx + k * Nx * Ny] = {{{0 9 18},{3 12 21},{6 15 24}},        */
/*                                     {{1 10 19},{4 13 22},{7 16 25}},       */
/*                                     {{2 11 20},{5 14 23},{8 17 26}}};      */
/*                                                                            */
/*    1D FFT ALONG X                                                          */
/*    Repeat Transpose 3: ALONG Z                                             */ 
/*      T[i + j * Nx * Nz + k * Nx] = U[i * Ny * Nz + j * Nz + k];            */
/*      T[i + j * Nx * Nz + k * Nx] = {{{0 1 2},{3 4 5},{6 7 8}},             */  
/*                                     {{9 10 11},{12 13 14;},{15 16 17}},    */
/*                                     {{18 19 20},{21 22 23},{24 25 26}}};   */
/* OR   T[i + j * Nx + k * Nx * Ny] = U[i * Ny * Nz + j * Nz + k];            */
/*      T[i + j * Nx + k * Nx * Ny] = {{{0 1 2},{3 4 5},{6 7 8}},             */  
/*                                     {{9 10 11},{12 13 14;},{15 16 17}},    */
/*                                     {{18 19 20},{21 22 23},{24 25 26}}};   */  
/******************************************************************************/
  PetscScalar U_array[Nz];
  PetscScalar AC_array[Nz];
/*************************************************************/
/*             FFTW 1: along Y axis - Columns arrangment     */
/*************************************************************/
/*                     TRANSPOSE 1: Z TO Y                   */
/*************************************************************/
  //ierr = SVL_3D_TRANSPOSE_1_REAL(Nx, Ny, Nz, U);CHKERRQ(ierr);
  ierr = SVL_3D_TRANSPOSE_1_COMPLEX(Nx, Ny, Nz, U);CHKERRQ(ierr);
/* Input columns */
  for (i = 0; i < Nx; i++) { // Loop through the height.
      for (j = 0; j < Ny; j++) { // Loop through the rows. 
// Build a vector of size Nz    
          for (k = 0; k < Nz; k++) { // Loop through the columns.
              U_array[k] = U[i * Ny * Nz + j * Nz + k];   
          }           
// Get 1D FFTW of an array of size Nz         
          ierr = SVL_1D_FFTW(Nx, Ny, Nz, U_array, AC_array);CHKERRQ(ierr);         
// Store the results in a new vector array.
          for (k = 0; k < Nz; k++) { // Loop through the rows.
              AC[i * Ny * Nz + j * Nz + k] = AC_array[k];
          }
      }
  }

  ierr = SAVE_1D_TO_3D_ARRAY_COMPLEX("OUTPUT_FFTW_REAL_1", "OUTPUT_FFTW_IMAG_1", Nx, Ny, Nz, AC);CHKERRQ(ierr);

/*************************************************************/
/*              FFTW 2: along Z axis  - Rows arrangment      */
/*************************************************************/
/*                     TRANSPOSE 2 : Y to Z                  */
/*************************************************************/  
 ierr = SVL_3D_TRANSPOSE_2_COMPLEX(Nx, Ny, Nz, AC);CHKERRQ(ierr);

/* Input rows */
  for (i = 0; i < Nx; i++) { // Loop through the height.
      for (j = 0; j < Ny; j++) { // Loop through the rows. 
// Build a vector of size Nz    
          for (k = 0; k < Nz; k++) { // Loop through the columns.
              U_array[k] = AC[i * Ny * Nz + j * Nz + k ]; 
          }  
// Get FFTW for the vector         
           ierr = SVL_1D_FFTW(Nx, Ny, Nz, U_array, AC_array);CHKERRQ(ierr);      
// Store the results in a new vector array.
          for (k = 0; k < Nz; k++) { // Loop through the rows.
              AC[i * Ny * Nz + j * Nz + k] = AC_array[k];
          }
      }
  }

  ierr = SAVE_1D_TO_3D_ARRAY_COMPLEX("OUTPUT_FFTW_REAL_2", "OUTPUT_FFTW_IMAG_2", Nx, Ny, Nz, AC);CHKERRQ(ierr);
/*************************************************************/
/*          FFTW 3: along X axis - Height arrangment         */
/*************************************************************/
/*                    TRANSPOSE 3 : Z to X                   */
/*************************************************************/
 ierr = SVL_3D_TRANSPOSE_3_COMPLEX(Nx, Ny, Nz, AC);CHKERRQ(ierr);

/* Input heights */
  for (i = 0; i < Nx; i++) { // Loop through the height.
      for (j = 0; j < Ny; j++) { // Loop through the rows. 
// Build a vector of size Nz    
          for (k = 0; k < Nz; k++) { // Loop through the columns.
              U_array[k] = AC[i * Ny * Nz + j * Nz + k];  
          }        
// Get FFTW for the vector       
          ierr = SVL_1D_FFTW(Nx, Ny, Nz, U_array, AC_array);CHKERRQ(ierr);           
// Store the results in a new vector array.
          for (k = 0; k < Nz; k++) { // Loop through the rows.
                AC[i * Ny * Nz + j * Nz + k]= AC_array[k];
         }
      }
  }

/*************************************************************/
/*                 REPEAD TRANSPOSE 3 : X TO Z               */
/*************************************************************/
  ierr = SVL_3D_TRANSPOSE_3_COMPLEX(Nx, Ny, Nz, AC);CHKERRQ(ierr);

/*************************************************************/
/*                        FFTW OUTPUT                        */
/*************************************************************/
  ierr = SAVE_1D_TO_3D_ARRAY_COMPLEX("OUTPUT_FFTW_REAL", "OUTPUT_FFTW_IMAG", Nx, Ny, Nz, AC);CHKERRQ(ierr);

/*************************************************************/
/*                      BACKWARD FFTW                        */
/*             FAST FOURIER TRANSFORM ON THE WEST            */
/*************************************************************/
   PetscScalar  Inv_U[Nx * Ny * Nz]; // Inv_U = INV_FFTW(U)
   ierr = SVL_3D_IFFTW(Nx, Ny, Nz, AC, Inv_U);CHKERRQ(ierr);
   ierr = SAVE_1D_TO_3D_ARRAY_COMPLEX("OUTPUT_INV_FFTW_REAL","OUTPUT_INV_FFTW_IMAG", Nx, Ny, Nz, Inv_U);CHKERRQ(ierr);

/*************************************************************/
/*                           SWAP                            */
/*************************************************************/
/*           TOP              3D           BOTTOM            */
/*         1    2                          5    6            */ 
/*         4    3                          8    7            */
/*      SWAP QUADRANTS (1 <--> 7, 2 <--> 8) DIAGONALLY       */
/*      SWAP QUADRANTS (3 <--> 5, 4 <--> 6) DIAGONALLY       */
/*                   --------------------                    */
/*                            2D                             */
/*                          1    2                           */ 
/*                          4    3                           */
/*      SWAP QUADRANTS (1 <--> 3, 2 <--> 4) DIAGONALLY       */
/*************************************************************/
   ierr = SVL_3D_SWAP_QUADRANTS(Nx, Ny, Nz, AC);CHKERRQ(ierr);
   ierr = SAVE_1D_TO_3D_ARRAY_COMPLEX("OUTPUT_SWAP_REAL","OUTPUT_SWAP_IMAG", Nx, Ny, Nz, AC);CHKERRQ(ierr);
 
  PetscFunctionReturn(0);
}/* END FUNCTION SVL SWAP QUADRANTS */
