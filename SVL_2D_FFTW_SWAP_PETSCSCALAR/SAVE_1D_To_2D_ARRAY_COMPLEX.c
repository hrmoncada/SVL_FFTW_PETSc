/*************************************************************/
/*      INPUT 2D ARRAY (==> PRINT ==>) OUTPUT 2D ARRAY       */
/*************************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <petsc.h>
# include "SVL_PETSC_Header.h"

#undef __FUNCT__
#define __FUNCT__ "SAVE_1D_To_2D_ARRAY_COMPLEX"
PetscErrorCode SAVE_1D_To_2D_ARRAY_COMPLEX(const char* desc_real, const char* desc_imag, PetscInt m, PetscInt n, PetscScalar* a) {
  PetscInt       i, j;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

 /* open file and write header */
  FILE *fp1, *fp2; 
  fp1 = fopen(desc_real, "w"); // save mat into a file as Square Unit Cell of zeros
  fp2 = fopen(desc_imag, "w"); // save mat into a file as Square Unit Cell of zeros

  //ierr = PetscPrintf(PETSC_COMM_WORLD,"%s, and %s\n", desc_real,desc_imag);CHKERRQ(ierr);

  for (i = 0; i < m; ++i) {
       for (j = 0; j < n; ++j) { 
            ierr = PetscPrintf(PETSC_COMM_WORLD,"%f + i %f  ", creal(a[i * n + j]), cimag(a[i * n + j]));CHKERRQ(ierr);
	    ierr = PetscFPrintf(PETSC_COMM_WORLD,fp1,"%f   ", creal(a[i * n + j]));CHKERRQ(ierr);
	    ierr = PetscFPrintf(PETSC_COMM_WORLD,fp2,"%f   ", cimag(a[i * n + j]));CHKERRQ(ierr);
	
       }
       ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);	
       ierr = PetscFPrintf(PETSC_COMM_WORLD,fp1,"\n");CHKERRQ(ierr);
       ierr = PetscFPrintf(PETSC_COMM_WORLD,fp2,"\n");CHKERRQ(ierr);
  }
  fclose(fp1);
  fclose(fp2);

  PetscFunctionReturn(0);
}
