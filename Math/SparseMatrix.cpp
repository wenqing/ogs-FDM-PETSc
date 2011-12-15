/*========================================================================
 GeoSys - class Matrix (Definition)
          class vec  
 Task:       Sparse matrix definition. 
 Function:   See the definition below
 programming:
  22/08/2004  WW  
==========================================================================*/

/// Matrix
#include "SparseMatrix.h"
#include <iomanip>
#include <float.h>
#include <cmath>
//

#include "SparseTable.h"
#include "misc.h"

namespace Math_Group
{

using namespace std;
/*\!
********************************************************************
   Constructor of sparse matrix
   Arguments:  
      sparse_table: Sparse graph
      dof:  Degree of freedom given by PDE        
   08/2007 WW
   10/2007 WW
********************************************************************/
SparseMatrix::SparseMatrix(const SparseTable &sparse_table, const int dof):DOF(dof)
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
SparseMatrix::SparseMatrix(const SparseMatrix& m)
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
SparseMatrix::~SparseMatrix()
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
double& SparseMatrix::operator() (const long i, const long j) const
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
	  k = AuxFunctions::binarySearch(entry_column, jr, num_column_entries[ir], num_column_entries[ir+1]); 
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
void SparseMatrix::operator = (const double a)
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
void SparseMatrix::operator *= (const double a)
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
void SparseMatrix::operator += (const double a)
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
void SparseMatrix::operator = (const SparseMatrix& m)
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
void SparseMatrix::operator += (const SparseMatrix& m)
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
void SparseMatrix::operator -= (const SparseMatrix& m)
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
void SparseMatrix::Write(ostream &os)
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
void SparseMatrix::Write_BIN(ostream &os)
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
void SparseMatrix::multiVec(double *vec_s, double *vec_r)
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

       if(symmetry)
       {
          /// ptr is num_column_entries
          for (ii = 0; ii < rows; ii++)
          {
             for (j = num_column_entries[ii]; j < num_column_entries[ii+1]; j++)
             {          
                jj=entry_column[j];
                vec_r[ii] += entry[j]*vec_s[jj];           
                if(ii!=jj)
                  vec_r[jj] += entry[j]*vec_s[ii];
             }         
          }

       }  
       else
       {
          /// ptr is num_column_entries
          for (ii = 0; ii < rows; ii++)
          {
             for (j = num_column_entries[ii]; j < num_column_entries[ii+1]; j++)
             {          
                jj=entry_column[j];
                vec_r[ii] += entry[j]*vec_s[jj];
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
void SparseMatrix::Trans_MultiVec(double *vec_s, double *vec_r)
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
   set
        A(ii,ii) = x_i, 
        A(ii, j) = 0., j!=ii
        A(i, ii) = 0., i!=ii
        b_i -= A(i,k)b_k  // b_k is given
   Programm:  
   10/2007 WW
   03/2011 WW  CRS storage
********************************************************************/
void SparseMatrix::Diagonize(const long idiag, const double b_given, double *b)
{
  //
  double vdiag = 0.;
  long j, k, ii, jj, j0;
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
	 /*
     long i;	  
     for (i = 0; i < rows; i++)
     {
        j = AuxFunctions::binarySearch(entry_column, id, num_column_entries[i], num_column_entries[i+1]); 
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
	 */
  }
  else if(storage_type == JDS)
  {
     long row_in_parse_table, counter;

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

#ifdef colDEBUG
     long i0;
     // /// Clean column id    
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
#endif
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
void SparseMatrix::Precond_Jacobi(double *vec_s, double *vec_r)
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

