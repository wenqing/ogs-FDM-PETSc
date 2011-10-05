/*========================================================================
 GeoSys - class Matrix, Sparse matrix (Declaration)   
 Task:       Matrix object for full matrices.  If the size of matrix is
             small, this class can do efficient matrix operation. 
 Function:   See the declaration below
 Design and programm WW
==========================================================================*/
#ifndef SparseMatrix_INC

#define SparseMatrix_INC

#include<iostream>
#include<fstream>



namespace Math_Group{

  using namespace std;
//

//08.2007 WW
//
class SparseMatrix
{
   public:
     SparseMatrix(const SparseTable &sparse_table, const int dof);
     SparseMatrix(const SparseMatrix& m); //17.03.2008
     ~SparseMatrix();
     // Preconditioner
     void Precond_Jacobi(double *vec_s, double *vec_r);

     // Operator
     void operator = (const double a);
     void operator *= (const double a);
     void operator += (const double a);
     void operator = (const SparseMatrix& m);
     void operator += (const SparseMatrix& m);
     void operator -= (const SparseMatrix& m);
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
