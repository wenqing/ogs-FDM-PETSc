/*========================================================================
 GeoSys - class Matrix, Sparse matrix (Declaration)   
 Task:       Matrix object for full matrices.  If the size of matrix is
             small, this class can do efficient matrix operation. 
 Function:   See the declaration below
 Design and programm WW
==========================================================================*/
#ifndef matrix_class_INC

#define matrix_class_INC

#include<iostream>
#include<fstream>

class FiniteDifference;
namespace Math_Group{
  using namespace std;
//

class Matrix
{
   public:
     Matrix(const int rows, const int cols=1);
     Matrix();
     explicit Matrix(const Matrix& m);
     //     
     void resize(const int rows, const int cols=1);
     //
     virtual ~Matrix();
     void ReleaseMemory(); //06.2010. WW
//----------------------------------------------
#ifdef OverLoadNEW_DELETE
     // Allocate memory
     void* operator new(size_t sz);
     // Release memory
     void operator delete(void* m);
#endif     
//----------------------------------------------
     // Operators
     virtual void operator = (const double a);
     virtual void operator *= (const double a);
     virtual void operator += (const double a);
     void operator = (const Matrix& m);
     void operator += (const Matrix& m);
     void operator -= (const Matrix& m);

     void GetTranspose(Matrix& m);

	 // vec_result = This*vec. vec_result must be initialized
     void multi(const double *vec, double *vec_result, const double fac=1.0);
     // m_result = this*m. m_result must be initialized
	 void multi(const Matrix& m, Matrix& m_result, const double fac=1.0);
	 // m_result = this*m1*m2. m_result must be initialized
     void multi(const Matrix& m1, const Matrix& m2, Matrix& m_result);

     // Access to members     
     virtual double& operator() (const int i, const int j=0) const;  
     void LimitSize(const int nRows, const int nCols=1);  

     int Rows() const {return nrows;}
     int Cols() const {return ncols;}
     int Size() const {return size;}

     // Print
     void Write(ostream& os=cout);
     void Write_BIN(fstream& os);
     void Read_BIN(fstream& is);
   protected:
     double *data;     
     int nrows, nrows0;
     int ncols, ncols0;  
     int size; 
     bool Sym;
};

// Symmetrical matrix. 12-01-2005. WW
class SymMatrix:public Matrix
{
   public:
     SymMatrix(const int dim);
     SymMatrix();
     explicit SymMatrix(const SymMatrix& m);
     
     void resize(const int dim);

     ~SymMatrix() {}

//----------------------------------------------
#ifdef OverLoadNEW_DELETE
     // Allocate memory
     void* operator new(size_t sz);
#endif     
//----------------------------------------------

     // Operators
     void operator = (const double a);
     void operator *= (const double a);
     void operator += (const double a);
     void operator = (const SymMatrix& m);
     void operator += (const SymMatrix& m);
     void operator -= (const SymMatrix& m);
     void LimitSize(const int dim);  

     // Access to members     
     double& operator() (const int i, const int j) const;  
};
typedef Matrix Vec;


/// Sparse matrix storage type
enum StorageType { CRS, JDS};
class SparseTable
{
    public:
      SparseTable(FiniteDifference *fdm, bool quadratic, bool symm=false, StorageType stype = JDS);
      ~SparseTable();   
      void Write(ostream &os=cout);    
    private:
      bool symmetry;
      // Topology mapping from data array to matrix
      long *entry_column;
      long *num_column_entries;     // number of entries of each columns in sparse table
      long *row_index_mapping_n2o;  // Row index of sparse table to row index of matrix 
      long *row_index_mapping_o2n;  // Inverse of last 
      long *diag_entry;             // Global index to the index of  entry_column
      long size_entry_column;
      long max_columns;  
      long rows;

      StorageType storage_type;  
      // Quick vector sort
      void Quicksort(long top, long bottom);
      long PartitionArray(long top, long bottom);
      // 
      friend class CSparseMatrix;
};

//08.2007 WW
//
class CSparseMatrix
{
   public:
     CSparseMatrix(const SparseTable &sparse_table, const int dof);
     CSparseMatrix(const CSparseMatrix& m); //17.03.2008
     ~CSparseMatrix();
     // Preconditioner
     void Precond_Jacobi(double *vec_s, double *vec_r);

     // Operator
     void operator = (const double a);
     void operator *= (const double a);
     void operator += (const double a);
     void operator = (const CSparseMatrix& m);
     void operator += (const CSparseMatrix& m);
     void operator -= (const CSparseMatrix& m);
     // Vector pass through augment and bring results back.
     void multiVec(double *vec_s, double *vec_r);
     void Trans_MultiVec(double *vec_s, double *vec_r);
     void Diagonize(const long idiag, const double b_given, double *b); 
     //
     // Access to members     
     double& operator() (const long i, const long j=0) const;  
     //
     long Dim() const {return DOF*rows;}
     int Dof() const {return DOF;}
     void SetDOF(const int dof_n) {DOF = dof_n;} //_new. 10/2008. WW
     long Size() const {return rows;}
     StorageType GetStorageType() const {return storage_type;}

     // Print
     void Write(ostream &os=cout);   
     void Write_BIN(ostream &os);   

   private:
     // Data
     double *entry;  
     mutable double zero_e; 

     /// 0. 03.2011. WW
     StorageType storage_type; 
     // 
     bool symmetry;
     // Topology mapping from data array to matrix. All are only pointers to the 
     // correpinding members in SparseTable, and no memory are allocated for them 
     long *entry_column;
     long *num_column_entries;  // number of entries of each columns in sparse table
     long *row_index_mapping_n2o;  // Row index of sparse table to row index of matrix 
     long *row_index_mapping_o2n;  // Inverse of last 
     long *diag_entry;
     long size_entry_column;
     long max_columns;  
     long rows;
     //
     int DOF;
     //
};
}
//==========================================================================

#endif
