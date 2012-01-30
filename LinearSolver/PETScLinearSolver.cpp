/*!
   \brief Definition of member functions of class PETScLinearSolver

   10~11.2011. WW

*/

#include<iostream>

#include "PETScLinearSolver.h"

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
   sol_type = lsol;
   pc_type = prec_type; 

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
  //MatSeqAIJSetPreallocation(A, nz ,PETSC_NULL);
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

      PetscPrintf(PETSC_COMM_WORLD,"\n================================================");         
      PetscPrintf(PETSC_COMM_WORLD, "\nLinear solver %s with %s preconditioner",
                                    sol_type.c_str(), pc_type.c_str() );         
      KSPGetIterationNumber(lsolver,&its); //CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_WORLD,"\nConvergence in %d iterations.\n",(int)its);
      PetscPrintf(PETSC_COMM_WORLD,"\n================================================");           }
   PetscPrintf(PETSC_COMM_WORLD,"\n");

   VecAssemblyBegin(x);
   VecAssemblyEnd(x);

   PetscGetTime(&v2);
   time_elapsed += v2-v1;

  /*
  //TEST
   PetscViewer viewer;
   PetscViewerASCIIOpen(PETSC_COMM_WORLD, "x2.txt", &viewer);
   PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
   PetscObjectSetName((PetscObject)x,"Solution");
   VecView(x, viewer);   
  */


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


void PETScLinearSolver::UpdateSolutions(PetscScalar *u0, PetscScalar *u1)
{
  int i, j;
  PetscScalar *xp;
 
  int receivecount;
  PetscInt low,high,otherlow;
  MPI_Status status; 
  PetscInt count;
  int tag = 9999;
  VecGetOwnershipRange(x, &low, &high);
  VecGetLocalSize(x, &count);


  VecGetArray(x, &xp);
  for(i=0; i<count; i++)
    u1[i] = xp[i];


  double *u_temp = new double[m_size];


  // Collect solution from processes.
  for(j=0; j<count; j++)
    u_temp[low+j] = u1[j];
  for(i=0;i<mpi_size;i++) 
  {
     if(i != rank)
     {

       MPI_Sendrecv( &count, 1, MPI_INT, i,tag, 
                     &receivecount,1,MPI_INT,i,tag, PETSC_COMM_WORLD ,&status);
       MPI_Sendrecv( &low, 1, MPI_INT, i,tag,
                 &otherlow,1,MPI_INT,i,tag,PETSC_COMM_WORLD,&status );
       MPI_Sendrecv( u1, count, MPI_DOUBLE, i,tag,
                     u0,receivecount,MPI_DOUBLE,i,tag, PETSC_COMM_WORLD,&status  );
       for(j=0;j<receivecount;j++)
         u_temp[otherlow+j] = u0[j];
     }
  }


  //MPI_Barrier(PETSC_COMM_WORLD);
  // Copy the collected solution to the array for the previous solution
  for(i=0;i<m_size;i++)
  {
    u1[i] = u_temp[i];
    u0[i] = u_temp[i];
  }
 

  delete [] u_temp;

  VecRestoreArray(x, &xp);

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
  MatZeroRows (A, nrows, rows, one, NULL, NULL);

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
