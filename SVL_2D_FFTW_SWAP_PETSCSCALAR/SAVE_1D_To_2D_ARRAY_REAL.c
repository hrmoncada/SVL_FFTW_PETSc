/*************************************************************/
/*      INPUT 2D ARRAY (==> PRINT ==>) OUTPUT 2D ARRAY       */
/*************************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <petsc.h>
# include "SVL_PETSC_Header.h"

#undef __FUNCT__
#define __FUNCT__ "SAVE_1D_To_2D_ARRAY_REAL"
PetscErrorCode SAVE_1D_To_2D_ARRAY_REAL (const char* desc, PetscInt m, PetscInt n, PetscReal* a) {
  PetscInt       i, j;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

 /* open file and write header */
  FILE *fp; 
  fp = fopen(desc, "w"); // save mat into a file as Square Unit Cell of zeros

  //ierr = PetscPrintf(PETSC_COMM_WORLD,"%s\n",desc);CHKERRQ(ierr);

  for (i = 0; i < m; ++i) {
       for (j = 0; j < n; ++j) { 
            //ierr = PetscPrintf(PETSC_COMM_WORLD,"%1.0f   ", a[i*n + j]);CHKERRQ(ierr);
	    ierr = PetscFPrintf(PETSC_COMM_WORLD,fp,"%f   ", a[i*n + j]);CHKERRQ(ierr);
	
       }
       //ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);	
       ierr = PetscFPrintf(PETSC_COMM_WORLD,fp,"\n");CHKERRQ(ierr);
  }
  fclose(fp);

  PetscFunctionReturn(0);
}
