#include <petsc/private/matimpl.h>        /*I "petscmat.h" I*/

#undef __FUNCT__
#define __FUNCT__ "MatGetRowSumAbs"
/*@
   MatGetRowSumAbs - Gets the sum of the absolute values of each row of the matrix

   Logically Collective on Mat and Vec

   Input Parameters:
.  mat - the matrix

   Output Parameter:
.  v - the vector for storing the sum of rows

   Level: intermediate

   Notes: This code is slow since it is not currently specialized for different formats

   Concepts: matrices^getting row sums

.seealso: MatGetDiagonal(), MatGetSubMatrices(), MatGetSubmatrix(), MatGetRowMax(), MatGetRowMin()
@*/
PETSC_EXTERN PetscErrorCode  MatGetRowSumAbs(Mat mat, Vec v)
{
  PetscInt       start = 0, end = 0, row;
  PetscScalar    *array;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(mat,MAT_CLASSID,1);
  PetscValidType(mat,1);
  PetscValidHeaderSpecific(v,VEC_CLASSID,2);
  if (!mat->assembled) SETERRQ(PetscObjectComm((PetscObject)mat),PETSC_ERR_ARG_WRONGSTATE,"Not for unassembled matrix");
  MatCheckPreallocated(mat,1);
  ierr = MatGetOwnershipRange(mat, &start, &end);CHKERRQ(ierr);
  ierr = VecGetArray(v, &array);CHKERRQ(ierr);
  for (row = start; row < end; ++row) {
    PetscInt          ncols, col;
    const PetscInt    *cols;
    const PetscScalar *vals;

    array[row - start] = 0.0;

    ierr = MatGetRow(mat, row, &ncols, &cols, &vals);CHKERRQ(ierr);
    for (col = 0; col < ncols; col++) {
      array[row - start] += abs(vals[col]);
    }
    ierr = MatRestoreRow(mat, row, &ncols, &cols, &vals);CHKERRQ(ierr);
  }
  ierr = VecRestoreArray(v, &array);CHKERRQ(ierr);
  ierr = PetscObjectStateIncrease((PetscObject) v);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

