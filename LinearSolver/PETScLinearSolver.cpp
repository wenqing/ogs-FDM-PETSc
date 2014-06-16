/*!
   \brief Definition of member functions of class PETScLinearSolver

   10~11.2011. WW

*/

#include<iostream>

#include "PETScLinearSolver.h"

#include <vector>

PETScLinearSolver:: ~PETScLinearSolver()
{
   VecDestroy(&b);
   VecDestroy(&x);
   MatDestroy(&A);
   if(lsolver) KSPDestroy(&lsolver);
   // if(prec) PCDestroy(&prec);

   PetscPrintf(PETSC_COMM_WORLD,"\n>>Number of Unknows: %d", m_size);
   PetscPrintf(PETSC_COMM_WORLD,"\n>>Elapsed time in linear solver: %fs", time_elapsed);
}

/*!
  \brief KSP and PC type

 KSPRICHARDSON "richardson"
 KSPCHEBYCHEV  "chebychev"
 KSPCG         "cg"
 KSPCGNE       "cgne"
 KSPNASH       "nash"
 KSPSTCG       "stcg"
 KSPGLTR       "gltr"
 KSPGMRES      "gmres"
 KSPFGMRES     "fgmres"
 KSPLGMRES     "lgmres"
 KSPDGMRES     "dgmres"
 KSPTCQMR      "tcqmr"
 KSPBCGS       "bcgs"
 KSPIBCGS        "ibcgs"
 KSPBCGSL        "bcgsl"
 KSPCGS        "cgs"
 KSPTFQMR      "tfqmr"
 KSPCR         "cr"
 KSPLSQR       "lsqr"
 KSPPREONLY    "preonly"
 KSPQCG        "qcg"
 KSPBICG       "bicg"
 KSPMINRES     "minres"
 KSPSYMMLQ     "symmlq"
 KSPLCD        "lcd"
 KSPPYTHON     "python"
 KSPBROYDEN    "broyden"
 KSPGCR        "gcr"
 KSPNGMRES     "ngmres"
 KSPSPECEST    "specest"

 PCNONE            "none"
 PCJACOBI          "jacobi"
 PCSOR             "sor"
 PCLU              "lu"
 PCSHELL           "shell"
 PCBJACOBI         "bjacobi"
 PCMG              "mg"
 PCEISENSTAT       "eisenstat"
 PCILU             "ilu"
 PCICC             "icc"
 PCASM             "asm"
 PCGASM            "gasm"
 PCKSP             "ksp"
 PCCOMPOSITE       "composite"
 PCREDUNDANT       "redundant"
 PCSPAI            "spai"
 PCNN              "nn"
 PCCHOLESKY        "cholesky"
 PCPBJACOBI        "pbjacobi"
 PCMAT             "mat"
 PCHYPRE           "hypre"
 PCPARMS           "parms"
 PCFIELDSPLIT      "fieldsplit"
 PCTFS             "tfs"
 PCML              "ml"
 PCPROMETHEUS      "prometheus"
 PCGALERKIN        "galerkin"
 PCEXOTIC          "exotic"
 PCHMPI            "hmpi"
 PCSUPPORTGRAPH    "supportgraph"
 PCASA             "asa"
 PCCP              "cp"
 PCBFBT            "bfbt"
 PCLSC             "lsc"
 PCPYTHON          "python"
 PCPFMG            "pfmg"
 PCSYSPFMG         "syspfmg"
 PCREDISTRIBUTE    "redistribute"
 PCSACUSP          "sacusp"
 PCSACUSPPOLY      "sacusppoly"
 PCBICGSTABCUSP    "bicgstabcusp"
 PCSVD             "svd"
 PCAINVCUSP        "ainvcusp"
 PCGAMG            "gamg"

*/
void PETScLinearSolver::Config(const float tol, const KSPType lsol, const PCType prec_type)
{
   ltolerance = tol;

   KSPCreate(PETSC_COMM_WORLD,&lsolver);
   KSPSetOperators(lsolver, A, A,DIFFERENT_NONZERO_PATTERN);
   KSPSetType(lsolver,lsol);

   KSPGetPC(lsolver, &prec);
   PCSetType(prec, prec_type); //  PCJACOBI); //PCNONE);
   KSPSetTolerances(lsolver,ltolerance,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
   KSPSetFromOptions(lsolver);

}
//-----------------------------------------------------------------
void PETScLinearSolver::VectorCreate(PetscInt m)
{
   //PetscErrorCode ierr;  // returned value from PETSc functions
   VecCreate(PETSC_COMM_WORLD, &b);
   VecSetSizes(b, PETSC_DECIDE, m);
   VecSetFromOptions(b);
   VecDuplicate(b, &x);

   VecGetLocalSize(x, &m_size_loc);
}


void PETScLinearSolver::MatrixCreate( PetscInt m, PetscInt n)
{
   // Number of nonzeros per row in DIAGONAL portion of
   // local submatrix (same value is used for all local rows)
   int d_nz;
   // Number of nonzeros per row in the OFF-DIAGONAL portion of
   // local submatrix (same value is used for all local rows).
   int o_nz;
   // Number of nonzeros per row (same for all rows)
   int nz;

   d_nz = 5;
   o_nz = 5;
   nz = 5;

   MatCreate(PETSC_COMM_WORLD, &A);
   MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m,n);
   MatSetFromOptions(A);
   MatSetType(A,MATMPIAIJ);
   MatMPIAIJSetPreallocation(A,d_nz,PETSC_NULL, o_nz,PETSC_NULL);
   MatSeqAIJSetPreallocation(A, nz ,PETSC_NULL);
   MatGetOwnershipRange(A,&i_start,&i_end);

   //TEST
   // std::cout<<"sub_a  "<<i_start<<";   sub_d "<<i_end<<std::endl;

}

void PETScLinearSolver::Solver()
{

   //TEST
   //PetscViewer viewer;
   //PetscViewerASCIIOpen(PETSC_COMM_WORLD, "x.txt", &viewer);
   //PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
   //PetscObjectSetName((PetscObject)x,"Solution");
   //VecView(x, viewer);



   int its;
   PetscLogDouble v1,v2;
   KSPConvergedReason reason;

   PetscGetTime(&v1);

   KSPSetOperators(lsolver, A, A,DIFFERENT_NONZERO_PATTERN);

   KSPSolve(lsolver, b, x);


   KSPGetConvergedReason(lsolver,&reason); //CHKERRQ(ierr);
   if (reason==KSP_DIVERGED_INDEFINITE_PC)
   {
      PetscPrintf(PETSC_COMM_WORLD,"\nDivergence because of indefinite preconditioner;\n");
      PetscPrintf(PETSC_COMM_WORLD,"Run the executable again but with -pc_factor_shift_positive_definite option.\n");
   }
   else if (reason<0)
   {
      PetscPrintf(PETSC_COMM_WORLD,"\nOther kind of divergence: this should not happen.\n");
   }
   else
   {
      const char *slv_type;
      const char *prc_type;
      KSPGetType(lsolver, &slv_type);
      PCGetType(prec, &prc_type);

      PetscPrintf(PETSC_COMM_WORLD,"\n================================================");
      PetscPrintf(PETSC_COMM_WORLD, "\nLinear solver %s with %s preconditioner",
                  slv_type, prc_type);
      KSPGetIterationNumber(lsolver,&its); //CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_WORLD,"\nConvergence in %d iterations.\n",(int)its);
      PetscPrintf(PETSC_COMM_WORLD,"\n================================================");
   }
   PetscPrintf(PETSC_COMM_WORLD,"\n");

   VecAssemblyBegin(x);
   VecAssemblyEnd(x);

   PetscGetTime(&v2);
   time_elapsed += v2-v1;

//#define  EXIT_TEST
#ifdef EXIT_TEST
   //TEST
   PetscViewer viewer;
   EQSV_Viewer("petsc", viewer);
//    PetscFinalize();
//    exit(0);
#endif
}

void PETScLinearSolver::AssembleRHS_PETSc()
{
   VecAssemblyBegin(b);
   VecAssemblyEnd(b);
}
void PETScLinearSolver::AssembleUnkowns_PETSc()
{
   VecAssemblyBegin(x);
   VecAssemblyEnd(x);
}
void PETScLinearSolver::AssembleMatrixPETSc()
{
   MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
}

void PETScLinearSolver::UpdateSolutions(PetscScalar *u)
{

#ifdef TEST_MEM_PETSC
   //TEST
   PetscLogDouble mem1, mem2;
   PetscMemoryGetCurrentUsage(&mem1);
#endif

   PetscScalar *xp = NULL; //nullptr;
   VecGetArray(x, &xp);

   gatherLocalVectors(xp, u);

   //MPI_Barrier(PETSC_COMM_WORLD);

   VecRestoreArray(x, &xp);

   //TEST
#ifdef TEST_MEM_PETSC
   PetscMemoryGetCurrentUsage(&mem2);
   PetscPrintf(PETSC_COMM_WORLD, "### Memory usage by Updating. Before :%f After:%f Increase:%d\n", mem1, mem2, (int)(mem2 - mem1));
#endif

}

void PETScLinearSolver::gatherLocalVectors( PetscScalar local_array[],
      PetscScalar global_array[])
{
   // Collect vectors from processors.
   int size_rank;
   MPI_Comm_size(PETSC_COMM_WORLD, &size_rank);

   // number of elements to be sent for each rank
   std::vector<PetscInt>  i_cnt(size_rank);
   // offset in the receive vector of the data from each rank
   std::vector<PetscInt>  i_disp(size_rank);

   MPI_Allgather(&m_size_loc, 1, MPI_INT, &i_cnt[0], 1, MPI_INT, PETSC_COMM_WORLD);

   // colloect local array
   PetscInt offset = 0;
   for(PetscInt i=0; i<size_rank; i++)
   {
      i_disp[i] = offset;
      offset += i_cnt[i];
   }

   MPI_Allgatherv(local_array, m_size_loc, MPI_DOUBLE,
                  global_array, &i_cnt[0], &i_disp[0], MPI_DOUBLE, PETSC_COMM_WORLD);

}

void PETScLinearSolver::set_bVectorEntry(const int i, const double value )
{

   VecSetValues(b,1,&i,&value,INSERT_VALUES);
}
void PETScLinearSolver::set_xVectorEntry(const int i, const double value)
{

   VecSetValues(x,1,&i,&value,INSERT_VALUES);
}


void PETScLinearSolver::add_bVectorEntry(const int i, const double value,InsertMode mode )
{

   VecSetValue(b, i, value, mode);
}
void PETScLinearSolver::add_xVectorEntry(const int i, const double value, InsertMode mode)
{

   VecSetValue(x, i, value,mode);
}

void PETScLinearSolver::Initialize( )
{

   VecSet(b, 0.0);
   VecSet(x, 0.0);
   MatZeroEntries(A);
}

void PETScLinearSolver::addMatrixEntry(const int i, const int j, const double value)
{

   MatSetValue(A, i, j, value, ADD_VALUES);
}

void PETScLinearSolver::addMatrixEntries(const int m,const int idxm[], const int n,
      const int idxn[],const PetscScalar v[])
{

   MatSetValues(A, m, idxm, n, idxn, v, ADD_VALUES);
}

void PETScLinearSolver::zeroRows_in_Matrix(const int nrows, const  PetscInt *rows)
{
   PetscScalar one = 1.0;
   // Each process indicates only rows it owns that are to be zeroed
   // MatSetOption(A, MAT_NO_OFF_PROC_ZERO_ROWS,PETSC_TRUE);

//#define TEST_SYM_MAT
#ifdef TEST_SYM_MAT
   if(nrows>0)
      MatZeroRowsColumns(A, nrows, rows, one, PETSC_NULL, PETSC_NULL);
   else
      MatZeroRowsColumns(A, 0, PETSC_NULL, one, PETSC_NULL, PETSC_NULL);
#else
   if(nrows>0)
      MatZeroRows (A, nrows, rows, one, PETSC_NULL, PETSC_NULL);
   else
      MatZeroRows(A, 0, PETSC_NULL, one, PETSC_NULL, PETSC_NULL);
#endif
}

void PETScLinearSolver::EQSV_Viewer(std::string file_name, PetscViewer viewer)
{
   std::string fname = file_name + "_eqs_dump.txt";
   PetscViewerASCIIOpen(PETSC_COMM_WORLD, fname.c_str(), &viewer);
   PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
   PetscObjectSetName((PetscObject)A,"Stiffness_matrix");
   PetscObjectSetName((PetscObject)b,"RHS");
   PetscObjectSetName((PetscObject)x,"Solution");
   MatView(A,viewer);
   VecView(b, viewer);
   VecView(x, viewer);
}
