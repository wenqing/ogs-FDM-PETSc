/*!
   \brief Declaration of class PETScLinearSolver

   11.2011. WW

*/
#ifndef PETSC_LSOLVER_INC
#define PETSC_LSOLVER_INC

#include <string>

#include "petscmat.h"
#include "petscksp.h"

typedef Mat PETSc_Mat;
typedef Vec PETSc_Vec;

class PETScLinearSolver
{
   public:
      PETScLinearSolver (const int size)
         :lsolver(NULL),prec(NULL)
      {
         ltolerance = 1.e-10;
         m_size = size;
         m_size_loc = PETSC_DECIDE;
         time_elapsed = 0.0;
      }
      ~PETScLinearSolver();

      void Config(const float tol, const KSPType lsol, const PCType prec_type);

      void Init()
      {

         VectorCreate(m_size);
         MatrixCreate(m_size, m_size);
      }

      void Solver();
      void AssembleRHS_PETSc();
      void AssembleUnkowns_PETSc( );
      void AssembleMatrixPETSc();

      void UpdateSolutions(PetscScalar *u);


      int Size() const
      {
         return m_size;
      }

      void set_xVectorEntry(const int i, const double value);
      void set_bVectorEntry(const int i, const double value);

      void add_xVectorEntry(const int i, const double value, InsertMode mode);
      void add_bVectorEntry(const int i, const double value, InsertMode mode);
      void addMatrixEntry(const int i, const int j, const double value);
      void addMatrixEntries(const int m,const int idxm[], const int n,
                            const int idxn[],const PetscScalar v[]);

      void Initialize();

      void zeroRows_in_Matrix(const int nrow, const  PetscInt *rows);
      void zeroMatrix()
      {
         MatZeroEntries(A);
      }


      void set_rank_size(const int m_rank, const int size)
      {
         mpi_size = size;
         rank = m_rank;
      }


      PetscInt getStartRow() const
      {
         return i_start;
      }
      PetscInt getEndRow() const
      {
         return i_end;
      }

      void EQSV_Viewer(std::string file_name, PetscViewer viewer);

   private:
      PETSc_Mat  A;
      PETSc_Vec b;
      PETSc_Vec x;
      KSP lsolver;
      PC prec;
      PetscInt i_start;
      PetscInt i_end;

      PetscLogDouble time_elapsed;

      PetscInt m_size;
      PetscInt m_size_loc;

      float ltolerance;

      int mpi_size;
      int rank;

      void VectorCreate(PetscInt m);
      void MatrixCreate(PetscInt m, PetscInt n);

      /*!
         \brief  collect local vectors
         \param  local_array  local array
         \param  global_array global array
      */
      void gatherLocalVectors(PetscScalar local_array[],
                              PetscScalar global_array[]);


};


#endif
