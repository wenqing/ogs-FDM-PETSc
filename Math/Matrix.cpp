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
#include "Matrix.h"
#include <iomanip>
#include <float.h>
#include <cmath>
//

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



}// Namespace

