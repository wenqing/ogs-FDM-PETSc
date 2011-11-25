/*========================================================================
 GeoSys - class Matrix, (Declaration)   
 Task:       Matrix object for full matrices.  If the size of matrix is
             small, this class can do efficient matrix operation. 
 Function:   See the declaration below
 Design and programm WW
==========================================================================*/
#ifndef Matrix_INC

#define Matrix_INC

#include<iostream>
#include<fstream>


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

     void getTranspose(Matrix& m);

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

typedef Matrix Vec;

}
//==========================================================================

#endif
