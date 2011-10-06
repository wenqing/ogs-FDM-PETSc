/**************************************************************************
Task: Sparse matrix and linear equation solver
      Design and program by WW
Programing:
10/2007 WW/
**************************************************************************/
#ifndef Linear_EQS_INC
#define Linear_EQS_INC
#include<vector>
#include<cmath>
//
#include "Matrix.h"
#include "SparseMatrix.h"

namespace _FDM{class Numerics; class FiniteDifference;}

//

namespace Math_Group
{
/// Using PDE related classes
using _FDM::Numerics;
using _FDM::FiniteDifference;
//
class Linear_EQS
{
  public:
    Linear_EQS(const SparseTable &sparse_table, 
               const long dof, bool messg=true);
#if defined(USE_MPI)
    Linear_EQS(const long size);    
#endif
    ~Linear_EQS();
    // Configure numerics
    void ConfigNumerics(Numerics *num, const long n=0);
    // Preconditioner;
    void Precond(double *vec_s, double *vec_r); 
    void TransPrecond(double *vec_s, double *vec_r); 
//02.2011. WW #if defined(USE_MPI)
    void Precond_Jacobi(const double *vec_s, double *vec_r);
//#endif
    //

    void ComputePreconditioner();
    void ComputePreconditioner_Jacobi();
    //   
    // Solver
    int Solver();
    int CG();
    int BiCG();
    int BiCGStab();
    int Gauss() {return -1;}
    int QMRCGStab() {return -1;}
    int CGNR() {return -1;}
    int CGS();
    int Richardson() {return -1;}
    int JOR() {return -1;}
    int SOR() {return -1;}
    int AMG1R5() {return -1;}
    int UMF() {return -1;}
    //   
    int GMRES();
    //
    void Initialize();
    void Clean(); 
    //
    // Access to the members
    void SetDOF(const int dof_n) {A->SetDOF(dof_n);} // For different processes with different DOF of OPDE. _new. 10/2008. WW
    void SetKnownX_i(const long i, const double x_i);
    double X(const long i) const {return x[i];} 
    double RHS(const long i) const {return b[i];} 
    double NormX();
    double NormRHS() { return bNorm; } 

    long Size() const {if(A) return A->Dim(); else return size_global;}
    // Write
    void Write(ostream &os=cout);    
    void Write_BIN(ostream &os);    
  private:
    SparseMatrix *A;
    SparseMatrix *M;  // Preconditioner;
    double *b;
    double *x;
    double *pivots;  // For ILU preconditioner. _new by WW 08.10.2008
    double *z_vec;   // For ILU preconditioner. _new by WW 08.10.2008

    /// GMRES. 30.06.2010. WW  
    /// GMRES H matrix     
    mutable Matrix H;
    int m_gmres; /// number of columns of H matrix
    inline void Update(double *x, int k, Matrix &h, double *s);
    inline void Get_Plane_Rotation(double &dx, double &dy, double &cs, double &sn);
    inline void Set_Plane_Rotation(double &dx, double &dy, double &cs, double &sn);
    //

    // 
    string solver_name;
    string precond_name;
    // Buffer
    vector<double *> f_buffer;
    // Controls
    int precond_type;
    int solver_type;
    bool message;
    int iter, max_iter;
    double tol, bNorm, error; 
    long size_global;
    // Operators
    double dot (const double *xx,  const double *yy); 
    inline double Norm(const double *xx)  { return sqrt(dot(xx, xx)); }
    inline bool CheckNormRHS(const double normb_new);
    //
    void Message();

    /// Using PDE related classes
    friend class _FDM::FiniteDifference;
      
//
};   
}
using Math_Group::Linear_EQS;

#endif
