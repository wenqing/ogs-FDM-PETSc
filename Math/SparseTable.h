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


namespace Math_Group{

  using _FDM::FiniteDifference;
  using namespace std;
//
  class SparseMatrix;

/// Sparse matrix storage type
enum StorageType { CRS, JDS};
class SparseTable
{
    public:
      SparseTable(FiniteDifference *fdm);
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
      // 
      friend class SparseMatrix;
};

};
}
//==========================================================================

#endif
