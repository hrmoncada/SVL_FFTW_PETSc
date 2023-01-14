/*************************************************************/
/*                     START MAIN LOOP                       */
/*************************************************************/

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <complex.h>
# include <petsc.h>
# include <petscksp.h>
# include <petscmat.h>
# include <petscvec.h>
# include "SVL_PETSC_3D_HEADER.h"

# undef __FUNCT__
# define __FUNCT__ "SVL_3D_LOOP"
PetscErrorCode SVL_3D_LOOP(char Switch_Cood_Sys, PetscInt M, PetscInt New_Nx, PetscInt New_Ny, PetscInt New_Nz, PetscInt NK, PetscScalar* AMNP, PetscReal* KX, PetscReal* KY, PetscReal* KZ, PetscReal* THETA, PetscReal* VARPHI, PetscReal* PER, Mat D) {
  PetscErrorCode     ierr;
  PetscInt           i, nk;
  PetscReal          Kx[New_Nx * New_Ny * New_Nz], Ky[New_Nx * New_Ny * New_Nz], Kz[New_Nx * New_Ny * New_Nz];
  PetscReal          f[3 * New_Nx * New_Ny * New_Nz]; // Ax = f, f is a column vector, dimension number: 2 for 2D and 3 for 3D 
  Vec                f_vec;
  Vec                B;
  Mat                C;
  KSP                ksp;

  PetscFunctionBeginUser;

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\nSTEP 13 :SPATIAL VARIANCE LOOP\n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n          TOTAL NUMBER OF SPATIAL HARMONICS ON THE 3D TRUNCATION ARRAY = %4d\n",NK);CHKERRQ(ierr);

/* OUTPUT */
  Vec PHI;
  ierr = VecCreate(PETSC_COMM_WORLD,&PHI);CHKERRQ(ierr);
  ierr = VecSetSizes(PHI,PETSC_DECIDE, M);CHKERRQ(ierr);
  ierr = VecSetFromOptions(PHI);CHKERRQ(ierr);

  Vec UC;
  ierr = VecCreate(PETSC_COMM_WORLD,&UC);CHKERRQ(ierr);
  ierr = VecSetSizes(UC,PETSC_DECIDE, M);CHKERRQ(ierr);
  ierr = VecSetFromOptions(UC);CHKERRQ(ierr);

/************************************************************************/
/*                           View  Matrix                               */
/************************************************************************/
   /*PetscInt          m, n;
   ierr = PetscPrintf(PETSC_COMM_WORLD,"\nView Matrix D on LOOP\n----------------\n");
   ierr = MatView(D,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
   ierr = MatGetSize(D,&m,&n);CHKERRQ(ierr);
   ierr = PetscPrintf(PETSC_COMM_WORLD,"\nD[cols, rows] = [%d %d]\n", m, n);
   ierr = PetscPrintf(PETSC_COMM_WORLD,"\nM = %d\n", M);
   ierr = PetscPrintf(PETSC_COMM_WORLD,"\nNK = %d\n", NK);*/

/*************************************************************/
/*                     START MAIN LOOP                       */
/*************************************************************/
   for (nk = 0; nk < NK; nk++) { 
   ierr = PetscPrintf(PETSC_COMM_WORLD,"\nnk = %d\n", nk);
/****************************************************************/
/*  STEP 1:     BUILDING THE GRADING VECTOR FOR EACH nk         */
/*  K(r) = 2*pi/Delta(r)[a_x*cos(THETA(r)) + a_y*sin(THETA(r))] */
/****************************************************************/  
     ierr = SVL_3D_ORIENTATION_VECTOR(New_Nx, New_Ny, New_Nz, Kx, Ky, Kz, KX[nk], KY[nk], KZ[nk]); CHKERRQ(ierr);

/*************************************************************/
/* C = CYLINDRICAL (RHO, THETA. ZH)                          */
/* S = SPHERICAL   (RHO, THETA, VARPHI)                      */
/*************************************************************/
/*             CARTESIAN TO SPHERICAL/CYLINDRIC              */
/*                            AND                            */
/*             SPHERICAL/CYLINDRIC  TO CARTESIAN             */
/*************************************************************/ 
/* Switch Coordinate System CYLINDRIC (C) & SPHERICAL(S) */
    switch(Switch_Cood_Sys) {
         case 'C':
            ierr = SVL_3D_CYLINDRICAL_TRANSLATION(New_Nx, New_Ny, New_Nz, Kx, Ky, Kz, THETA, PER); CHKERRQ(ierr);
           break;
         case 'S':
            ierr = SVL_3D_SPHERICAL_TRANSLATION(New_Nx, New_Ny, New_Nz, Kx, Ky, Kz, THETA, VARPHI, PER); CHKERRQ(ierr);
           break;
        default:
           ierr = PetscPrintf(PETSC_COMM_WORLD,"\nTRY AGAIN!!!, YOU DID NOT IMPLEMENTED ANY ATTRIBUTES \n");CHKERRQ(ierr);
           break;
   }

/*************************************************************/
/* Attribute 3 -  Fill Fraction Function                     */ 
/* WRITING f AS COLUMN ARRAY, f = [Kx(:); Ky(:); Kz(:)]      */
/* each Kx column is order one under the other,              */
/* next Ky and Kz columns are order on the same way          */   
/* RSH : Define Rigth Hand Side Linear System    (A x = f)   */   
/*************************************************************/
      ierr = SVL_3D_RHS(New_Nx, New_Ny, New_Nz, Kx, Ky, Kz, f);CHKERRQ(ierr);

/************************************************************************/
/*                           Creat  Vector  f_vec                       */
/************************************************************************/
 /*  Create vector f_vect using f[]*/
     ierr = VecCreate(PETSC_COMM_WORLD,&f_vec);CHKERRQ(ierr);
     ierr = VecSetSizes(f_vec,PETSC_DECIDE, 3*M);CHKERRQ(ierr); /* change the number: 2 for 2D and 3 for 3D */
     ierr = VecSetFromOptions(f_vec);CHKERRQ(ierr);

/* Join the elements of vector array "f_vec" into a vector "vec"  */
     for (i = 0; i < 3*M; i++) {/* Careful here set the appropied dimensions: "2" for 2D and "3" for 3D */
         PetscScalar v = f[i];
         ierr = VecSetValues(f_vec,1,&i,&v,INSERT_VALUES);CHKERRQ(ierr); 
     }

  /* Assemble vector */
      ierr = VecAssemblyBegin(f_vec);CHKERRQ(ierr);
      ierr = VecAssemblyEnd(f_vec);CHKERRQ(ierr);

/************************************************************************/
/*                           View  Vector  f_vec                       */
/************************************************************************/
     /*PetscInt vsize;
     ierr = PetscPrintf(PETSC_COMM_WORLD,"\nView vector f_vec\n----------------\n");
     ierr = VecView(f_vec,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
     ierr = VecGetSize(f_vec,&vsize);CHKERRQ(ierr);
     ierr = PetscPrintf(PETSC_COMM_WORLD, "\nf_vec size = %d \n", vsize);CHKERRQ(ierr);*/

/************************************************************************/
/*                           Creat  Vector  B                           */
/************************************************************************/
     ierr = VecCreate(PETSC_COMM_WORLD,&B);CHKERRQ(ierr);
     ierr = VecSetSizes(B,PETSC_DECIDE,M);CHKERRQ(ierr);
     ierr = VecSetFromOptions(B);CHKERRQ(ierr);

/************************************************************************/
/*                           View  Vector                               */
/************************************************************************/
     /*PetscInt Bsize;
     ierr = PetscPrintf(PETSC_COMM_WORLD,"\nView vector B\n----------------\n");
     ierr = VecView(B,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
     ierr = VecGetSize(B,&Bsize);CHKERRQ(ierr);
     ierr = PetscPrintf(PETSC_COMM_WORLD, "\nB size = %d \n\n", Bsize);CHKERRQ(ierr);
     ierr = PetscPrintf(PETSC_COMM_WORLD, "\nNK 4 = %d \n", nk);CHKERRQ(ierr);*/

/************************************************************************/
/*          Matrix-vector multiplication B = D^T*f                      */
/************************************************************************/     
     ierr =  MatMultTranspose(D,f_vec,B);CHKERRQ(ierr);

/* Extended processs: 1. Find transpose matrix D^T, 2. Next multiple D^T*f */ 
     //Mat Dtrans;
     //ierr = MatTranspose(D,MAT_INITIAL_MATRIX, &Dtrans);CHKERRQ(ierr);
     //ierr = MatMult(Dtrans,f_vec,B);CHKERRQ(ierr);

/************************************************************************/
/*                           View  Matrix                               */
/************************************************************************/
     //ierr = PetscPrintf(PETSC_COMM_WORLD,"\nView Matrix D^T\n----------------\n");
     //ierr = MatView(Dtrans,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
     //ierr = MatGetSize(Dtrans,&m,&n);CHKERRQ(ierr);
     //ierr = PetscPrintf(PETSC_COMM_WORLD,"\nD^T[cols, rows] = [%d %d]\n", m, n);
     /*  PetscInt a = 1; 
   if( a == 1) {  
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Exiting the program....\n");CHKERRQ(ierr);
      exit(1);
   } */
/************************************************************************/
/*          Matrix-Matrix multiplication C = D^T*D                      */
/************************************************************************/
     MatTransposeMatMult(D, D, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &C);

     //ierr = PetscPrintf(PETSC_COMM_WORLD, "\nNK 5 = %d \n", nk);CHKERRQ(ierr);
/************************************************************************/
/*                         View  Matrix  C = D^T*D                      */
/************************************************************************/
    /* ierr = PetscPrintf(PETSC_COMM_WORLD,"\nView Matrix C = D^T*D\n----------------\n");
     ierr = MatView(C,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
     ierr = MatGetSize(C,&m,&n);CHKERRQ(ierr);
     ierr = PetscPrintf(PETSC_COMM_WORLD,"\nC[cols, rows] = [%d %d]\n", m, n);
*/
/************************************************************************/  
/* SOLVE LINEAR EQUATION   C X = B  , D  is mxn over determined matrix  */
/*                                   (more rows than columns )          */ 
/* Least square:  C X = B => (D^T D) X = (D^T b)                        */ 
/*                                   X = (D^T D)^(-1) (D^T b)           */
/************************************************************************/    

/* Create linear solver context */
    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);

/* Set operators. Here the matrix that defines the linear system  also serves as the preconditioning matrix. ksp = C */
    ierr = KSPSetOperators(ksp,C,C);CHKERRQ(ierr);

/* Set runtime options, e.g., -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
   These options will override those specified above as long as KSPSetFromOptions() is called _after_ any other customization routines.*/
    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

/* Solve linear system.  Here we explicitly call KSPSetUp() for more detailed performance monitoring of certain preconditioners, such
   as ICC and ILU.  This call is optional, as KSPSetUp() will  automatically be called within KSPSolve() if it hasn't been called already. */
    ierr = KSPSetUp(ksp);CHKERRQ(ierr);
    
/************************************************************************/ 
/* Least square:  C X = B => (D^T D) X = (D^T b) ,  ksp = C             */ 
/*                                                                      */
/* Solve : ksp * PHI = B   ==>  PHI = ksp\B                             */
/************************************************************************/ 
    ierr = KSPSolve(ksp,B,PHI);CHKERRQ(ierr); 

/* View solver info; we could instead use the option -ksp_view to print this info to the screen at the conclusion of KSPSolve(). */
  //ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  
/* OUTPUT PHI */
    /*if (nk == 0) { 
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\nOUTPUT_PHI_0\n");
        PetscViewer viewer;
        ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "OUTPUT_PHI_0.mat", &viewer);CHKERRQ(ierr);
       // ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_DENSE);CHKERRQ(ierr);
        ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_DENSE);CHKERRQ(ierr);
        ierr = VecView(PHI,viewer);CHKERRQ(ierr);
        ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&viewer);
     }*/
    if (nk == NK-1) { 
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\nOUTPUT_PHI\n");
        PetscViewer viewer;
        ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "OUTPUT_PHI.mat", &viewer);CHKERRQ(ierr);
       // ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_DENSE);CHKERRQ(ierr);
        ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_DENSE);CHKERRQ(ierr);
        ierr = VecView(PHI,viewer);CHKERRQ(ierr);
        ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&viewer);
     }

/*************************************************************/    
/*           SOLVE the exponential S = AMN[nk]*exp(i*PHI)    */
/*************************************************************/

/* PETSC_i : the imaginary number i */
     PetscScalar imag_num = 1.0*I;

/* VecScale : Scales a vector  x[i] = alpha * x[i] => PHI_i = PETS_i * PHI_i = I*PHI_i, for i=1,...,n.. */
     VecScale(PHI,imag_num);

/* VecExp: Replaces each component of a vector by  e^x_i => PHI_i =  e^(I*PHI_i) */
     VecExp(PHI);

/* VecScale : Scales a vector  x[i] = alpha * x[i]  => S = PHI = AMN[nk] * e^(I*PHI),    for i=1,...,n.. */
     VecScale(PHI, AMNP[nk]); //AMN[nk] is a complex number; 
 
/* OUTPUT S */
     /*if (nk == 0) {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"OUTPUT_S_0\n");
        PetscViewer viewer;
        ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "OUTPUT_S_0.mat", &viewer);CHKERRQ(ierr);
       // ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_DENSE);CHKERRQ(ierr);
        ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_DENSE);CHKERRQ(ierr);
        ierr = VecView(PHI,viewer);CHKERRQ(ierr);
        ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&viewer);
     }*/
     if (nk == NK-1) {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"OUTPUT_S\n");
        PetscViewer viewer;
        ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "OUTPUT_S.mat", &viewer);CHKERRQ(ierr);
       // ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_DENSE);CHKERRQ(ierr);
        ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_DENSE);CHKERRQ(ierr);
        ierr = VecView(PHI,viewer);CHKERRQ(ierr);
        ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&viewer);
     }

/*************************************************************/    
/*           SOLVE the exponential UC = UC + S               */
/*************************************************************/
     PetscScalar alpha = 1.0;

/* VecAXPY: Computes  y = alpha x + y => UC =  alpha * (AMN[nk]*PHI) + UC = alpha * S + UC */
     VecAXPY(UC, alpha , PHI);

/* VecAYPX: Computes  y = x + alpha y => UC = (AMN[nk]*PHI) +  alpha * UC = S + alpha* UC */
     //VecAYPX(UC, alpha , PHI);

     ierr = VecAssemblyBegin(UC);CHKERRQ(ierr);
     ierr = VecAssemblyEnd(UC);CHKERRQ(ierr);
   
/* OUTPUT UC */
     /*if (nk == 0) { 
        ierr = PetscPrintf(PETSC_COMM_WORLD,"OUTPUT_UC_0\n");
        PetscViewer viewer;
        ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "OUTPUT_UC_0.mat", &viewer);CHKERRQ(ierr);
       // ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_DENSE);CHKERRQ(ierr);
        ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_DENSE);CHKERRQ(ierr);
        ierr = VecView(UC,viewer);CHKERRQ(ierr);
        ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&viewer);
     }*/
     if (nk == NK-1) { 
        ierr = PetscPrintf(PETSC_COMM_WORLD,"OUTPUT_UC\n");
        PetscViewer viewer;
        ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "OUTPUT_UC.mat", &viewer);CHKERRQ(ierr);
       // ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_DENSE);CHKERRQ(ierr);
        ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_DENSE);CHKERRQ(ierr);
        ierr = VecView(UC,viewer);CHKERRQ(ierr);
        ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&viewer);
     }

   } /* End nk loop array */

    ierr = PetscPrintf(PETSC_COMM_WORLD, "\nnk = %d \n\n", nk);CHKERRQ(ierr);

/*************************************************************/
/*                      END MAIN LOOP                       */
/*************************************************************/

/* Free spaces : Clean up memory usage.*/ 
   //VecDestroy(f_vec);
   //VecDestroy(B);
    ierr = VecDestroy(&PHI);CHKERRQ(ierr);
    //ierr = VecDestroy(&UC);CHKERRQ(ierr);
    ierr = MatDestroy(&D);CHKERRQ(ierr);
    ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
   //MatDestroy(C);


   PetscFunctionReturn(0);  
}
