
static char help[] = "help !";

#include <slepceps.h>
#include <petscblaslapack.h>
#include <assert.h>
#include <sys/time.h>

/*
   User-defined routines
*/

static int
countbits(long v)
{
    int c;
    for(c = 0; v; v >>= 1)
    {
      c += v & 1;
    }
    return c;
}


//
// From:
// https://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetNaive
//

static const unsigned char BitsSetTable256[256] = 
{
#   define B2(n) n,     n+1,     n+1,     n+2
#   define B4(n) B2(n), B2(n+1), B2(n+1), B2(n+2)
#   define B6(n) B4(n), B4(n+1), B4(n+1), B4(n+2)
    B6(0), B6(1), B6(1), B6(2)
};

static int
countbits_fast(long v0)
{
    uint32_t v;
    assert(0xffffffffL&v0 == v0);
    v = (uint32_t)v0;

    uint32_t c; // c is the total bits set in v
    c = BitsSetTable256[v & 0xff] + 
        BitsSetTable256[(v >> 8) & 0xff] + 
        BitsSetTable256[(v >> 16) & 0xff] + 
        BitsSetTable256[v >> 24]; 
    return c;
}

#include "body.h"

#undef __FUNCT__
#define __FUNCT__ "MatMult_Hamiltonian"
/*
    Computes y <-- A*x, where A is the block tridiagonal matrix
 */
PetscErrorCode MatMult_Hamiltonian(Mat A,Vec x,Vec y)
{
  void              *ctx;
  long               nx,lo,i,j;
  const PetscScalar *px;
  PetscScalar       *py;
  PetscErrorCode    ierr;

  PetscFunctionBeginUser;
  ierr = MatShellGetContext(A,&ctx);CHKERRQ(ierr);
  nx = *(long*)ctx;
  ierr = VecGetArrayRead(x,&px);CHKERRQ(ierr);
  ierr = VecGetArray(y,&py);CHKERRQ(ierr);

    printf("/"); fflush(stdout);
    matmult(py, px, nx);
    printf("\\"); fflush(stdout);

  ierr = VecRestoreArrayRead(x,&px);CHKERRQ(ierr);
  ierr = VecRestoreArray(y,&py);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__

void
dump(nx)
{
    PetscScalar *px, *py;
    long i, j;
    px = (PetscScalar *)alloca(sizeof(PetscScalar)*nx);
    py = (PetscScalar *)alloca(sizeof(PetscScalar)*nx);
    for(i=0; i<nx; i++)
    {
        memset(px, 0, sizeof(PetscScalar)*nx);
        memset(py, 0, sizeof(PetscScalar)*nx);
        px[i] = 1.;

        matmult(py, px, nx);
        for(j=0; j<nx; j++)
            printf("%.0f ", py[j]);
        printf("\n");
    }
}


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  Mat            A;               /* operator matrix */
  EPS            eps;             /* eigenproblem solver context */
  EPSType        type;
  PetscMPIInt    size;
  PetscInt       n,nev;
  PetscBool      terse;
  PetscErrorCode ierr;

    n = DIMS;

  SlepcInitialize(&argc,&argv,(char*)0,help);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  if (size != 1) SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor example only");

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nHamiltonian Eigenproblem (matrix-free version), n=%D\n\n",n);CHKERRQ(ierr);

//    dump(DIMS);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the operator matrix that defines the eigensystem, Ax=kx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = MatCreateShell(PETSC_COMM_WORLD,n,n,n,n,&n,&A);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatShellSetOperation(A,MATOP_MULT,(void(*)())MatMult_Hamiltonian);CHKERRQ(ierr);
//  ierr = MatShellSetOperation(A,MATOP_MULT_TRANSPOSE,(void(*)())MatMult_Hamiltonian);CHKERRQ(ierr);
//  ierr = MatShellSetOperation(A,MATOP_GET_DIAGONAL,(void(*)())MatGetDiagonal_Hamiltonian);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create eigensolver context
  */
  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);

  /*
     Set operators. In this case, it is a standard eigenvalue problem
  */
  ierr = EPSSetOperators(eps,A,NULL);CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRQ(ierr);

  /*
     Set solver parameters at runtime
  */
  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);

//    printf("EPSSetInitialSpace\n");
//    Vec iv;
//    VecCreateSeq(PETSC_COMM_SELF, n, &iv);
//    ierr = EPSSetInitialSpace(eps, 1, &iv);CHKERRQ(ierr);
//    VecDestroy(&iv);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    printf("EPSSolve\n");
  ierr = EPSSolve(eps);CHKERRQ(ierr);

  /*
     Optional: Get some information from the solver and display it
  */
  ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
  ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* show detailed info unless -terse option is given by user */
  ierr = PetscOptionsHasName(NULL,NULL,"-terse",&terse);CHKERRQ(ierr);
  if (terse) {
    ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL);CHKERRQ(ierr);
  } else {
    ierr = PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);CHKERRQ(ierr);
    ierr = EPSReasonView(eps,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }
  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

