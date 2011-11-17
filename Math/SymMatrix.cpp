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
#include "SymMatrix.h"
#include <iomanip>
#include <float.h>
#include <cmath>
//

#include "Matrix.h"


namespace Math_Group
{

// Symmetrical matrix 
SymMatrix::SymMatrix(const int dim):Matrix(0)
{
    Sym = true; 
    nrows = ncols = dim;    
    size = (int)nrows*(nrows+1)/2;
    data = new double[size];
    nrows0 = ncols0 = dim;
    for(int i=0; i<size; i++) data[i] = 0.0;
}

SymMatrix::SymMatrix():Matrix(0)
{
   Sym = true;
   nrows = 0;
   ncols = 0;
   nrows0 = 0;
   ncols0 = 0;
   size = 0;
   data = 0;
}
SymMatrix::SymMatrix(const SymMatrix& m):Matrix(0)
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

void SymMatrix::resize(const int dim)
{
 
   if(size>0)
   {
      delete [] data;
      data = NULL;
   }
     
    Sym = true; 
    nrows = ncols = dim;    
    size = (int)nrows*(nrows+1)/2;
    data = new double[size];
    nrows0 = ncols0 = dim;
    for(int i=0; i<size; i++) data[i] = 0.0;
}


//----------------------------------------------
#ifdef OverLoadNEW_DELETE
void* SymMatrix::operator new(size_t sz) {
  //printf("operator new: %d Bytes\n", sz);
  void* m = malloc(sz);
  if(!m) puts("out of memory");
  return m;
}
#endif


void SymMatrix::operator = (const double a)
{
    for(int i=0; i<size; i++) data[i] = a;
}
void SymMatrix::operator *= (const double a)
{
    for(int i=0; i<size; i++) data[i] *= a;
}
void SymMatrix::operator += (const double a)
{
    for(int i=0; i<size; i++) data[i] += a;
}


//
void SymMatrix::operator = (const SymMatrix& m)
{
#ifdef gDEBUG    
    if(nrows!=m.Rows()||ncols!=m.Cols())
    {
        cout<<"\n The sizes of the two matrices are not matched"<<endl;
        abort();
    }
#endif
    int id=0;
    for(int i=0; i<nrows; i++)
    {
       for(int j=0; j<ncols; j++)
       {
          if(j>i) continue;
          id = (int)i*(i+1)/2+j; // temporary  
          data[id] = m(i,j);
       }
    }
}

//
void SymMatrix::operator += (const SymMatrix& m)
{
#ifdef gDEBUG    
   if(nrows!=m.Rows())
    {
        cout<<"\n The sizes of the two matrices are not matched"<<endl;
        abort();
    }
#endif
    int id=0;
    for(int i=0; i<nrows; i++)
    {
       for(int j=0; j<ncols; j++)
       {
          if(j>i) continue;
          id = (int)i*(i+1)/2+j; // temporary  
          data[id] += m(i,j);
       }
    }
}

//
void SymMatrix::operator -= (const SymMatrix& m)
{
#ifdef gDEBUG    
    if(nrows!=m.Rows()) //Assertion, will be removed
    {
        cout<<"\n The sizes of the two matrices are not matched"<<endl;
        abort();
    }
#endif
    int id=0;
    for(int i=0; i<nrows; i++)
    {
       for(int j=0; j<ncols; j++)
       {
          if(j>i) continue;
          id = (int)i*(i+1)/2+j; // temporary  
          data[id] -= m(i,j);
       }
    }

}
//
double& SymMatrix::operator() (const int i, const int j) const
{
 #ifdef gDEBUG    
    if(i>=nrows||j>=nrows)
    {
        cout<<"\n Index exceeds the size of the matrix"<<endl;
        abort();
    }
 #endif
    
    int id=0;
    if(i>=j)
       id = (int)i*(i+1)/2+j; // temporary
    else
       id = (int)j*(j+1)/2+i; // temporary
    return data[id]; 
} 

void  SymMatrix::LimitSize(const int dim)
{
 #ifdef gDEBUG    
    if(dim>nrows0)
    {
        cout<<"\n Given size exceeds the original size of the matrix"<<endl;
        abort();
    }
 #endif
    nrows = ncols = dim;
    size = (int)nrows*(nrows+1)/2;
}




}// Namespace

