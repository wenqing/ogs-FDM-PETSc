/*========================================================================
 GeoSys - class SymMatrix   
 Task:       Symetrical matrix object for full matrices.  If the size of matrix is
             small, this class can do efficient matrix operation. 
 Function:   See the declaration below
 Design and programm WW
==========================================================================*/
#ifndef matrix_class_INC

#define matrix_class_INC

#include<iostream>
#include<fstream>



namespace Math_Group{

  using namespace std;
//

class Matrix;

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
}
//==========================================================================

#endif
