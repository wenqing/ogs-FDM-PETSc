/*========================================================================
 GeoSys - class Matrix (Definition)
          class vec  
 Task:       Matrix object for full matrices.  If the size of matrix is
             small, this class can do efficient matrix operation. 
 Function:   See the definition below
 programming:
  22/08/2004  WW  
==========================================================================*/

/// Matrix
#include <iomanip>
#include <float.h>
#include <cmath>
//
#include "matrix_class.h"

///PDE related classes:
#include "fdm.h"
#include "geo.h"


namespace Math_Group{

////PDE related class:
//using _FDM::Point;

// Constructors
Matrix::Matrix(const int rows, const int cols)
{    
  if(rows*cols>0)
  {
    Sym = false;
    nrows = rows;
    ncols = cols;
    nrows0 = rows;
    ncols0 = ncols;
    size = nrows*ncols;
    data = new double[size];
    for(int i=0; i<size; i++) data[i] = 0.0;
  }
}
Matrix::Matrix()
{
   Sym = false;
   nrows = 0;
   ncols = 0;
   nrows0 = 0;
   ncols0 = 0;
   size = 0;
   data = 0;
}
Matrix::Matrix(const Matrix& m)
{
   Sym = m.Sym;
   nrows = m.nrows;
   ncols = m.ncols;
   nrows0 = m.nrows0;
   ncols0 = m.ncols0;
   size = m.size;
   data = new double[size];
   for(int i=0; i<size; i++) data[i] = 0.0;
}

void Matrix::resize(const int rows, const int cols)
{
 
   if(size>0) { 
      delete [] data;
      data = NULL;
   }
     
   if(rows*cols>0)
   {
      Sym = false;
      nrows = rows;
      ncols = cols;
      nrows0 = rows;
      ncols0 = ncols;
      size = nrows*ncols;
      data = new double[size];
      for(int i=0; i<size; i++) data[i] = 0.0;
   }
}

Matrix::~Matrix()
{
    delete [] data;
    data = NULL;
}
// 06.2010. WW
void Matrix::ReleaseMemory()
{
    delete [] data;
    data = NULL;
}

//----------------------------------------------
#ifdef OverLoadNEW_DELETE
void* Matrix::operator new(size_t sz) {
  //printf("operator new: %d Bytes\n", sz);
  void* m = malloc(sz);
  if(!m) puts("out of memory");
  return m;
}

void Matrix::operator delete(void* m) {

  Matrix* mm = static_cast<Matrix*>(m); 
  free(mm); 
}
#endif
//----------------------------------------------
//
void Matrix::operator = (const double a)
{
    for(int i=0; i<size; i++) data[i] = a;
}
void Matrix::operator *= (const double a)
{
    for(int i=0; i<size; i++) data[i] *= a;
}
void Matrix::operator += (const double a)
{
    for(int i=0; i<size; i++) data[i] += a;
}
//
void Matrix::operator = (const Matrix& m)
{
#ifdef gDEBUG    
    if(nrows!=m.Rows()||ncols!=m.Cols())
    {
        cout<<"\n The sizes of the two matrices are not matched"<<endl;
        abort();
    }
#endif
    for(int i=0; i<nrows; i++) 
      for(int j=0; j<ncols; j++) 
         data[i*ncols+j] = m(i,j);
}

//
void Matrix::operator += (const Matrix& m)
{
#ifdef gDEBUG    
   if(nrows!=m.Rows()||ncols!=m.Cols())
    {
        cout<<"\n The sizes of the two matrices are not matched"<<endl;
        abort();
    }
#endif
    for(int i=0; i<nrows; i++) 
      for(int j=0; j<ncols; j++) 
         data[i*ncols+j] += m(i,j);
}

//
void Matrix::operator -= (const Matrix& m)
{
#ifdef gDEBUG    
    if(nrows!=m.Rows()||ncols!=m.Cols()) //Assertion, will be removed
    {
        cout<<"\n The sizes of the two matrices are not matched"<<endl;
        abort();
    }
#endif
    for(int i=0; i<nrows; i++) 
      for(int j=0; j<ncols; j++) 
         data[i*ncols+j] -= m(i,j);
}
//
void Matrix::GetTranspose(Matrix& m)
{
 #ifdef gDEBUG    
    if(ncols!=m.Rows()&&nrows!=m.Cols())
    {
        cout<<"\n The sizes of the two matrices are not matched"<<endl;
        abort();
    }
 #endif

   for(int i=0; i<m.Rows(); i++)
       for(int j=0; j<m.Cols(); j++)
       {
//          m(i,j) = data[j*ncols+i];
           m(i,j) = (*this)(j,i);
       }

}
//
// m_results = this*m. m_results must be initialized
void Matrix::multi(const Matrix& m, Matrix& m_result, const double fac)
{
 #ifdef gDEBUG    
    if(ncols!=m.Rows()&&nrows!=m_result.Rows()&&m.Cols()!=m_result.Cols())
    {
        cout<<"\n The sizes of the two matrices are not matched"<<endl;
        abort();
    }
 #endif
    for(int i=0; i<m_result.Rows(); i++) 
    {
       for(int j=0; j<m_result.Cols(); j++) 
       { 
           if(Sym&&(j>i)) continue;
           // m_result(i,j) = 0.0;
           for(int k=0; k<ncols; k++)
//            m_result(i,j) += fac*data[i*ncols+k]*m(k,j);
              m_result(i,j) += fac*(*this)(i,k)*m(k,j);
           
       }
    }
}

//
// m_results = this*m1*m2. m_results must be  initialized
void Matrix::multi(const Matrix& m1, const Matrix& m2, Matrix& m_result)
{
  #ifdef gDEBUG    
    if(ncols!=m1.Rows()&&m1.Cols()!=m2.Rows()
       &&m2.Cols()!=m_result.Cols()&&nrows!=m_result.Rows())
    {
        cout<<"\n The sizes of the two matrices are not matched"<<endl;
        abort();
    }
 #endif
    for(int i=0; i<m_result.Rows(); i++) 
    {
       for(int j=0; j<m_result.Cols(); j++) 
       { 
           if(Sym&&(j>i)) continue;
           //m_result(i,j) = 0.0;
           for(int k=0; k<ncols; k++)
           {
              for(int l=0; l<m2.Rows(); l++) 
//                m_result(i,j) += data[i*ncols+k]*m1(k,l)*m2(l,j);
                m_result(i,j) += (*this)(i,k)*m1(k,l)*m2(l,j);
           }
       }
    }
}
// vec_result = This*vec. vec_result must be  initialized
void Matrix::multi(const double *vec, double *vec_result, const double fac)
{
    for(int i=0; i<nrows; i++)
	{
       for(int j=0; j<ncols; j++)
//         vec_result[i] += fac*data[i*ncols+j]*vec[j];
         vec_result[i] += fac*(*this)(i,j)*vec[j];
    }
}

double& Matrix::operator() (const int i, const int j) const
{
 #ifdef gDEBUG    
    if(i>=nrows||j>=ncols)
    {
        cout<<"\n Index exceeds the size of the matrix"<<endl;
        abort();
    }
 #endif
    return data[i*ncols+j]; 
} 
void  Matrix::LimitSize(const int nRows, const int nCols)
{
 #ifdef gDEBUG    
    if(nRows>nrows0||nCols>ncols0)
    {
        cout<<"\n Given size exceeds the original size of the matrix"<<endl;
        abort();
    }
 #endif
    nrows = nRows;
	ncols = nCols;
    size = nrows*ncols;
}


/**************************************************************************
MathLib-Method: 
Task: 
Programing:
08/2004 WW Implementation
02/2005 WW Change name
**************************************************************************/
void Matrix::Write(ostream& os)
{
    //os<<"============================================="<<endl;
    //os<<"Rows: "<<Rows()<<"  Columns: "<<Cols()<<endl;
 
    os.setf(ios::scientific,ios::floatfield);
    os.precision(12);

    for(int i=0; i<nrows; i++)
    {
       os<< "| ";
       for(int j=0; j<ncols; j++)
         os<<(*this)(i,j)<<" ";
       os<< "| "<<endl;
    }
    os<<endl;
    //os<<"============================================="<<endl;
    //os<<endl;     
}

/**************************************************************************
MathLib-Method: 
Task: 
Programing:
01/2006 WW Implementation
**************************************************************************/
void Matrix::Write_BIN(fstream& os)
{
    for(int i=0; i<size; i++) 
      os.write((char*)(&data[i]), sizeof(data[i]));
}
/**************************************************************************
MathLib-Method: 
Task: 
Programing:
01/2006 WW Implementation
**************************************************************************/
void Matrix::Read_BIN(fstream& is)
{
    for(int i=0; i<size; i++) 
      is.read((char*)(&data[i]), sizeof(data[i]));
}


/*\!
********************************************************************
   Create sparse matrix table
   01/2006 WW
   08/2007 WW
   10/2007 WW
   03/2010 WW: CRS storage
********************************************************************
*/
SparseTable::SparseTable(FiniteDifference *fdm)
{
  
   long i=0, j=0;
   long **larraybuffer;
   larraybuffer = NULL;
   storage_type = CRS;
   symmetry = false;  
 
   //
   rows = (long)fdm->grid_point_in_use.size();

   size_entry_column = 0;
   diag_entry = new long[rows]; 

   row_index_mapping_n2o = NULL; 
   row_index_mapping_o2n = NULL;
 
   
   /// CRS storage
   /// num_column_entries saves vector ptr of CRS 
   num_column_entries = new long[rows+1];

   vector<long> A_index;
   long col_index;

   for(i=0; i<rows; i++)
   {
      num_column_entries[i] = (long)A_index.size();
 
      for(j=0; j<fdm->grid_point_in_use[i]->GetNumNeighborPoints(); j++)
      {
         col_index = fdm->grid_point_in_use[i]->GetNeighborIndex(j);
                        
         if(i == col_index)
            diag_entry[i] = (long)A_index.size();
         A_index.push_back(col_index);
      }
   }
      
   size_entry_column = (long)A_index.size(); 
   num_column_entries[rows] = size_entry_column;

   entry_column = new long[size_entry_column]; 
   for(i=0; i<size_entry_column; i++)
      entry_column[i] = A_index[i];

}
/*\!
********************************************************************
   Create sparse matrix table
   08/2007 WW
   10/2007 WW
********************************************************************/
void SparseTable::Write(ostream &os)
{
   long i, k, counter=0;

   os.width(10);
   os<<"Symmetry: "<<symmetry<<endl;
   os<<"\n*** Row index  "<<endl;
 
   if(storage_type == CRS)
   {
      os<<"\n*** Sparse entry  "<<endl;
      for (i = 0; i < rows; i++)
      {
          for (k = num_column_entries[i]; k < num_column_entries[i+1]; k++)
             os<<entry_column[k]+1<<" ";
          os<<endl; 
      }
   } 
   else if(storage_type == JDS)
   {
      for (i = 0; i < rows; i++)
        os<<row_index_mapping_n2o[i]+1<<endl;
      // 
      os<<"\n*** Sparse entry  "<<endl;
      for (k = 0; k < max_columns; k++)
      {
         os<<"--Column: "<<k+1<<endl;
         for (i = 0; i < num_column_entries[k]; i++)
         {          
            os<<entry_column[counter]+1<<endl;;
            counter++;
         } 
         os<<endl;        
      } 
   } 
}    

/*\!
********************************************************************
   Create sparse matrix table
   08/2007 WW
   10/2007 WW
********************************************************************/
SparseTable::~SparseTable()
{
  if(entry_column) delete [] entry_column;
  if(num_column_entries) delete [] num_column_entries; 
  if(row_index_mapping_n2o) delete [] row_index_mapping_n2o;    
  if(row_index_mapping_o2n) delete [] row_index_mapping_o2n;    
  if(diag_entry) delete [] diag_entry;
  entry_column = NULL;
  num_column_entries = NULL; 
  row_index_mapping_n2o = NULL;    
  row_index_mapping_o2n = NULL;    
  diag_entry = NULL;
}
/*\!
********************************************************************
   Constructor of sparse matrix
   Arguments:  
      sparse_table: Sparse graph
      dof:  Degree of freedom given by PDE        
   08/2007 WW
   10/2007 WW
********************************************************************/
CSparseMatrix::CSparseMatrix(const SparseTable &sparse_table, const int dof):DOF(dof)
{
  symmetry = sparse_table.symmetry;
  size_entry_column = sparse_table.size_entry_column;
  max_columns = sparse_table.max_columns;  
  rows = sparse_table.rows;
  storage_type = sparse_table.storage_type;
  // Topology mapping from data array to matrix
  // Only refer address
  entry_column = sparse_table.entry_column;
  num_column_entries = sparse_table.num_column_entries;  
  row_index_mapping_n2o = sparse_table.row_index_mapping_n2o;        
  row_index_mapping_o2n = sparse_table.row_index_mapping_o2n;        
  diag_entry = sparse_table.diag_entry;
  // Values of all sparse entries    
  entry = new double[dof*dof*size_entry_column+1];
  entry[dof*dof*size_entry_column] = 0.;
  zero_e = 0.;
  //
}
/*\!
********************************************************************
   Constructor of sparse matrix
   Arguments:  
      sparse_table: Sparse graph
      dof:  Degree of freedom given by PDE        
   03/2008 WW
********************************************************************/
CSparseMatrix::CSparseMatrix(const CSparseMatrix& m)
{

  symmetry = m.symmetry;
  size_entry_column = m.size_entry_column;
  max_columns = m.max_columns;  
  rows = m.rows;
  storage_type = m.storage_type;

  // Topology mapping from data array to matrix
  // Only refer address
  entry_column = m.entry_column;
  num_column_entries = m.num_column_entries;  
  row_index_mapping_n2o = m.row_index_mapping_n2o;        
  row_index_mapping_o2n = m.row_index_mapping_o2n;        
  diag_entry = m.diag_entry;
  DOF = m.DOF;
  // Values of all sparse entries    
  entry = new double[DOF*DOF*size_entry_column+1];
  for(long i=0; i<DOF*DOF*size_entry_column; i++)
    entry[i] = m.entry[i];
  //
  entry[DOF*DOF*size_entry_column] = 0.;
  zero_e = 0.;   
}

/*\!
********************************************************************
   Desstructor of sparse matrix
   08/2007 WW
   10/2007 WW
   03/2011 WW 
********************************************************************/
CSparseMatrix::~CSparseMatrix()
{
  delete [] entry;
  entry = NULL;
}
/*\!
********************************************************************
   Operator of sparse matrix
   08/2007 WW
   10/2008 WW
********************************************************************/
double& CSparseMatrix::operator() (const long i, const long j) const
{
 #ifdef gDEBUG    
  if(i>=rows*DOF||j>=rows*DOF)
  {
    cout<<"\n Index exceeds the dimension of the matrix"<<endl;
    abort();
  }
 #endif    
  long ii, jj, ir, jr, k;
  ii = i;
  jj = j;
  if(symmetry)
  {
    if(i>j)
    {
      ii = j;
      jj = i;
    }       
  }
  ir = ii%rows;
  jr = jj%rows;
  ii /= rows;
  jj /= rows;
  //
  k = -1;

  if(storage_type==JDS)
  { 
     long row_in_parse_table, counter;
     row_in_parse_table = row_index_mapping_o2n[ir];
     counter = row_in_parse_table;
     for (k = 0; k < max_columns; k++)
     {
        if(row_in_parse_table>=num_column_entries[k]) 
          return zero_e;
        if(entry_column[counter]==jr)
          break;  // Found the entry  
        counter += num_column_entries[k]; 
     }
     if(counter>=size_entry_column)
       return zero_e;
     //  Zero entry;  
     k = (ii*DOF+jj)*size_entry_column+counter;
  }
  else if(storage_type==CRS)
  {
     /// Left boundary of this row: num_column_entries[ir]
     /// Right boundary of this row: num_column_entries[ir+1]
     /// Search target is jr
     k = binarySearch(entry_column, jr, num_column_entries[ir], num_column_entries[ir+1]); 
     if(k==-1)
       return zero_e; 

     k = (ii*DOF+jj)*size_entry_column+k;
  }
  
  return entry[k]; // 
} 

/*\!
********************************************************************
   Desstructor of sparse matrix
   08/2007 WW
   10/2007 WW
********************************************************************/
void CSparseMatrix::operator = (const double a)
{
   long size = DOF*DOF*size_entry_column;
   for(long i=0; i<size; i++)
     entry[i] = a; 
}
/*\!
********************************************************************
   Desstructor of sparse matrix
   08/2007 WW
   10/2007 WW
********************************************************************/
void CSparseMatrix::operator *= (const double a)
{
   long size = DOF*DOF*size_entry_column;
   for(long i=0; i<size; i++)
     entry[i] *= a; 
}
/*\!
********************************************************************
   Desstructor of sparse matrix
   08/2007 WW
   10/2007 WW
********************************************************************/
void CSparseMatrix::operator += (const double a)
{
   long size = DOF*DOF*size_entry_column;
   for(long i=0; i<size; i++)
     entry[i] += a; 
}
/*\!
********************************************************************
   Desstructor of sparse matrix
   08/2007 WW
   10/2007 WW
********************************************************************/
void CSparseMatrix::operator = (const CSparseMatrix& m)
{
   long size = DOF*DOF*size_entry_column;
 #ifdef gDEBUG    
   if(size!=m.DOF*m.DOF*m.size_entry_column)
   {
    cout<<"\n Dimensions of two matrices do not match"<<endl;
    abort();
   }
#endif
   for(long i=0; i<size; i++)
     entry[i] = m.entry[i]; 
}
/*\!
********************************************************************
   Desstructor of sparse matrix
   08/2007 WW
   10/2007 WW
********************************************************************/
void CSparseMatrix::operator += (const CSparseMatrix& m)
{
   long size = DOF*DOF*size_entry_column;
 #ifdef gDEBUG    
   if(size!=m.DOF*m.DOF*m.size_entry_column)
   {
    cout<<"\n Dimensions of two matrices do not match"<<endl;
    abort();
   }
#endif
   for(long i=0; i<size; i++)
     entry[i] += m.entry[i]; 
}
/*\!
********************************************************************
   Desstructor of sparse matrix
   08/2007 WW
   10/2007 WW
********************************************************************/
void CSparseMatrix::operator -= (const CSparseMatrix& m)
{
   long size = DOF*DOF*size_entry_column;
 #ifdef gDEBUG    
   if(size!=m.DOF*m.DOF*m.size_entry_column)
   {
    cout<<"\n Dimensions of two matrices do not match"<<endl;
    abort();
   }
#endif
   for(long i=0; i<size; i++)
     entry[i] -= m.entry[i]; 
}
/*!
********************************************************************
   Output sparse matrix
   08/2007 WW
   10/2007 WW
   03/2011 WW  CRS
********************************************************************/
void CSparseMatrix::Write(ostream &os)
{
  //
  long i, k, ii, jj, row_in_parse_table, counter;
  os<<"*** Non-zero entries of matrix:  "<<endl;
  os.width(14);
  os.precision(8); 
  // 
  if(storage_type == CRS )
  {
     for(ii=0; ii<DOF; ii++)
     {
       for(i=0; i<rows; i++)
       {
         for(jj=0; jj<DOF; jj++)
         {
           for (k = num_column_entries[i]; k < num_column_entries[i+1]; k++)
           {
//TEST
                // if(fabs(entry[(ii*DOF+jj)*size_entry_column+counter])>DBL_MIN) //DBL_EPSILON)
               os<<setw(10)<<ii*rows+i<<" "
                 <<setw(10)<< jj*rows+entry_column[k]<<" "
                 <<setw(15)<<entry[(ii*DOF+jj)*size_entry_column+k]<<endl;  
            }
         }
       }
     }
  }
  else if(storage_type == JDS )
  {
     for(ii=0; ii<DOF; ii++)
     {
       for(i=0; i<rows; i++)
       {
         row_in_parse_table = row_index_mapping_o2n[i];
         for(jj=0; jj<DOF; jj++)
         {
           counter = row_in_parse_table;
           for (k = 0; k < max_columns; k++)
           {
             if(row_in_parse_table<num_column_entries[k])
             {
//TEST
                // if(fabs(entry[(ii*DOF+jj)*size_entry_column+counter])>DBL_MIN) //DBL_EPSILON)
               os<<setw(10)<<ii*rows+i<<" "
                 <<setw(10)<< jj*rows+entry_column[counter]<<" "
                 <<setw(15)<<entry[(ii*DOF+jj)*size_entry_column+counter]<<endl;  
               counter += num_column_entries[k];
             }
             else
               break; 
           }
         }
       }
     }
   }
}
//--------------------------------------------------------------
/*!
   \brief Write matrix to a binary file
   
   03.2011. WW
*/
void CSparseMatrix::Write_BIN(ostream &os)
{
  if(storage_type == JDS )
     return; 
  //
  if(DOF == 1)
  {
     os.write((char*) &rows, sizeof(long));
     os.write((char*) num_column_entries, (rows+1)*sizeof(long));
     os.write((char*) entry_column, (num_column_entries[rows]-1)*sizeof(long));
     os.write((char*) entry, (num_column_entries[rows]-1)*sizeof(double));
  }
  else
  {
     long i, k, ii, jj, size;
     long *ptr;
     long *A_index;
     double *A_value;
     
     ptr = new long[DOF*rows+1];
     size = DOF*DOF*num_column_entries[rows];
     A_index = new long[size];
     A_value = new double[size];

     long counter = 0;
     
     for(ii=0; ii<DOF; ii++)
     {
       for(i=0; i<rows; i++)
       {
         ptr[ii*rows+i] = counter;
         for(jj=0; jj<DOF; jj++)
         {
           for (k = num_column_entries[i]; k < num_column_entries[i+1]; k++)
           {
               A_index[counter] = jj*rows+entry_column[k];
               A_value[counter] = entry[(ii*DOF+jj)*size_entry_column+k];
               counter++; 
            }
         }
       }
     }
     ptr[DOF*rows] = counter;

     ii = DOF*rows;
     os.write((char*) &ii, sizeof(long));
     os.write((char*) ptr, (DOF*rows+1)*sizeof(long));
     os.write((char*) A_index, size*sizeof(long));
     os.write((char*) A_value, size*sizeof(double));

     delete [] ptr;
     delete [] A_index;
     delete [] A_value;

  }
}

/*!
********************************************************************
   Perform A*x
   Arguments:  
      vec_sr: M*vec_s-->vec_r
   01/2006 WW
   08/2007 WW
   10/2007 WW      //Jagged diagonal storage is not suitable for 
                     multithread simulation
   03/2011 WW      CRS storage
********************************************************************/
void CSparseMatrix::multiVec(double *vec_s, double *vec_r)
{
  long i, j, k, ii, jj, kk,ll,idof, jdof, counter;
  for(i=0; i<rows*DOF; i++)
    vec_r[i] = 0.0;
  //
  counter=0;
  if(DOF>1)
  {
    // Although this piece of code can deal with the case
    // of DOF = 1, we also prepare a special piece of code for
    // the case of DOF = 1 just for efficiency
    if(storage_type==CRS)
    {
       /// ptr is num_column_entries
       for (ii = 0; ii < rows; ii++)
       {
          for (j = num_column_entries[ii]; j < num_column_entries[ii+1]; j++)
          {          
             jj=entry_column[j];
             for(idof=0; idof<DOF; idof++)
             {
                kk = idof*rows+ii;
                for(jdof=0; jdof<DOF; jdof++)
                {
                  ll = jdof*rows+jj; 
                  k = (idof*DOF+jdof)*size_entry_column+j;
                  vec_r[kk] += entry[k]*vec_s[ll];
                  if(symmetry&(kk!=ll))
                     vec_r[ll] += entry[k]*vec_s[kk];
                }
             }
          }         
       }
        
    }
    else if(storage_type==JDS)
    { 
       for (k = 0; k < max_columns; k++)
       {
          for (i = 0; i < num_column_entries[k]; i++)
          {          
             ii = row_index_mapping_n2o[i];  
             jj=entry_column[counter];
             for(idof=0; idof<DOF; idof++)
             {
                kk = idof*rows+ii;
                for(jdof=0; jdof<DOF; jdof++)
                {
                  ll = jdof*rows+jj; 
                  j = (idof*DOF+jdof)*size_entry_column+counter;
                  vec_r[kk] += entry[j]*vec_s[ll];
                  if(symmetry&(kk!=ll))
                     vec_r[ll] += entry[j]*vec_s[kk];
                }
             }
             counter++;
          }         
       }
    }

  }
  else  // DOF = 1
  {
    if(storage_type==CRS)
    {
       /// ptr is num_column_entries
       for (ii = 0; ii < rows; ii++)
       {
          for (j = num_column_entries[ii]; j < num_column_entries[ii+1]; j++)
          {          
             jj=entry_column[j];
             vec_r[ii] += entry[j]*vec_s[jj];
             if(symmetry&(ii!=jj))
                 vec_r[jj] += entry[j]*vec_s[ii];
          }         
       }
        
    }
    else if(storage_type==JDS)
    {
       for (k = 0; k < max_columns; k++)
       {
          for (i = 0; i < num_column_entries[k]; i++)
          {          
             ii = row_index_mapping_n2o[i];  
             jj=entry_column[counter];
             vec_r[ii] += entry[counter]*vec_s[jj];
             if(symmetry&(ii!=jj))
                vec_r[jj] += entry[counter]*vec_s[ii];
             counter++;
          }         
       }
     }
  }
}

/*\!
********************************************************************
   Perform A^T*x
   Arguments:  
      vec_sr: M^T*vec_s-->vec_r
   08/2008 WW
********************************************************************/
void CSparseMatrix::Trans_MultiVec(double *vec_s, double *vec_r)
{
  long i, j, k, ii, jj, kk,ll,idof, jdof, counter;
  for(i=0; i<rows*DOF; i++)
    vec_r[i] = 0.0;
  //
  counter=0;
  if(DOF>1)
  {
    // Although this piece of code can deal with the case
    // of DOF = 1, we also prepare a special piece of code for
    // the case of DOF = 1 just for efficiency
    if(storage_type==CRS)
    {
       /// ptr is num_column_entries
       for (ii = 0; ii < rows; ii++)
       {
          for (j = num_column_entries[ii]; j < num_column_entries[ii+1]; j++)
          {          
             jj=entry_column[j];
             for(idof=0; idof<DOF; idof++)
             {
                kk = idof*rows+ii;
                for(jdof=0; jdof<DOF; jdof++)
                {
                  ll = jdof*rows+jj; 
                  k = (idof*DOF+jdof)*size_entry_column+j;
                  vec_r[ll] += entry[k]*vec_s[kk];
                  if(symmetry&(kk!=ll))
                     vec_r[kk] += entry[k]*vec_s[ll];
                }
             }
          }         
       }
        
    }
    else if(storage_type==JDS)
    {
       for (k = 0; k < max_columns; k++)
       {
          for (i = 0; i < num_column_entries[k]; i++)
          {          
             ii = row_index_mapping_n2o[i];  
             jj=entry_column[counter];
             for(idof=0; idof<DOF; idof++)
             {
                kk = idof*rows+ii;
                for(jdof=0; jdof<DOF; jdof++)
                {
                  ll = jdof*rows+jj; 
                  j = (idof*DOF+jdof)*size_entry_column+counter;
                  vec_r[ll] += entry[j]*vec_s[kk];
                  if(symmetry&(kk!=ll))
                     vec_r[kk] += entry[j]*vec_s[ll];
                }
             }
             counter++;
          }         
       }
     }
  }
  else  // DOF = 1
  {
    if(storage_type==CRS)
    {
       /// ptr is num_column_entries
       for (ii = 0; ii < rows; ii++)
       {
          for (j = num_column_entries[ii]; j < num_column_entries[ii+1]; j++)
          {          
             jj=entry_column[j];
             vec_r[jj] += entry[j]*vec_s[ii];
             if(symmetry&(ii!=jj))
                 vec_r[ii] += entry[j]*vec_s[jj];
          }         
       }
        
    }
    else if(storage_type==JDS)
    {
       for (k = 0; k < max_columns; k++)
       {
          for (i = 0; i < num_column_entries[k]; i++)
          {          
             ii = row_index_mapping_n2o[i];  
             jj=entry_column[counter];
             vec_r[jj] += entry[counter]*vec_s[ii];
             if(symmetry&(ii!=jj))
                vec_r[ii] += entry[counter]*vec_s[jj];
             counter++;
          }         
       }
    }
  }
}
//
/*\!
********************************************************************
   Set
        A(ii,ii) = x_i, 
        A(ii, j) = 0., j!=ii
        A(i, ii) = 0., i!=ii
        b_i -= A(i,k)b_k  // b_k is given
   Programm:  
   10/2007 WW
   03/2011 WW  CRS storage
********************************************************************/
void CSparseMatrix::Diagonize(const long idiag, const double b_given, double *b)
{
  //
  double vdiag = 0.;
  long i, j, k, ii, jj, j0;
  long id = idiag%rows;
  ii = idiag/rows;

  if(storage_type == CRS)
  {
     /// Diagonal entry and the row where the diagonal entry exists
     j = diag_entry[id]; 
     vdiag = entry[(ii*DOF+ii)*size_entry_column+j];
     /// Row where the diagonal entry exists 
     for(jj=0; jj<DOF; jj++)
     {
        for(k=num_column_entries[id]; k<num_column_entries[id+1]; k++)
        {
           j0=entry_column[k];
           if(id==j0&&jj==ii)  // Diagonal entry
              continue; 
           entry[(ii*DOF+jj)*size_entry_column+k] = 0.; 
        }
     }    
     
     /// Clean column id                
     for (i = 0; i < rows; i++)
     {
        j = binarySearch(entry_column, id, num_column_entries[i], num_column_entries[i+1]); 
        if(j == -1)
           continue;
        j0=entry_column[j];

        for(jj=0; jj<DOF; jj++)
        {
           if(i == j0&&ii==jj) continue; 
           k = (jj*DOF+ii)*size_entry_column+j;
           b[jj*rows+i] -= entry[k]*b_given; 
           entry[k] = 0.; 
           // Room for symmetry case
        }           
     }      
  }
  else if(storage_type == JDS)
  {
     long i0, row_in_parse_table, counter;

     // Row is zero
     row_in_parse_table = row_index_mapping_o2n[id];
     counter = row_in_parse_table;
     for (k = 0; k < max_columns; k++)
     {
        if(row_in_parse_table<num_column_entries[k])
        {
           j0=entry_column[counter];
           for(jj=0; jj<DOF; jj++)
           {
              if(id==j0&&jj==ii)   
                vdiag = entry[(ii*DOF+jj)*size_entry_column+counter];
              else  
                entry[(ii*DOF+jj)*size_entry_column+counter] = 0.;
           }
           counter += num_column_entries[k];
        }
        else
           break;
     }
     //
     counter=0;
     for (k = 0; k < max_columns; k++)
     {
       for (i = 0; i < num_column_entries[k]; i++)
       {          
         i0 = row_index_mapping_n2o[i]; 
         /*
         if(i0 == id)
         {
           counter++;
           continue; 
         }
         */
         j0=entry_column[counter];
         if(j0 == id)
         {
            for(jj=0; jj<DOF; jj++)
            {
               if(i0 == j0&&ii==jj) continue; 
               j = (jj*DOF+ii)*size_entry_column+counter;
               b[jj*rows+i0] -= entry[j]*b_given; 
               entry[j] = 0.; 
               // Room for symmetry case
            }           
         } 
         //
         counter++;
       }         
     }
  }
  b[idiag] = vdiag*b_given;
}

/*\!
********************************************************************
   M^{-1}*A
   
          a_ij  i=j
   M = { 
          0     i!=j 
   Programm:  
   10/2007 WW
********************************************************************/
void CSparseMatrix::Precond_Jacobi(double *vec_s, double *vec_r)
{
  long i, j, idof;
  double diag = 0.;
  //
  if(DOF>1)
  {
    // Although this piece of code can deal with the case
    // of DOF = 1, we also prepare a special piece of code for
    // the case of DOF = 1 just for efficiency
    for(i=0; i<rows; i++)
    {
       for(idof=0; idof<DOF; idof++)
       {
          diag = entry[(idof*DOF+idof)*size_entry_column+diag_entry[i]];
          if(fabs(diag)<DBL_MIN)
//        if(fabs(diag)<DBL_EPSILON)
            diag = 1.0;
          //  cout<<"Diagonal entry is zero. Abort simulation!!  " <<endl;
          j = idof*rows+i;
          vec_r[j] = vec_s[j]/diag;
       }         
    }
    // 
  }
  else  // DOF = 1
  {
    for(i=0; i<rows; i++)
    {
       diag = entry[diag_entry[i]];
       //if(fabs(diag)<DBL_EPSILON)
       if(fabs(diag)<DBL_MIN)
          diag = 1.0;
       //   cout<<"Diagonal entry is zero. Abort simulation!!  " <<endl;
       //
       vec_r[i] = vec_s[i]/diag;
    }
  }
}



}// Namespace

