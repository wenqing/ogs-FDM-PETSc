/**************************************************************************
Task: Linear equation
Programing:
11/2007 WW/
**************************************************************************/

//#include"matrix_class.h"
#include"LinearEQS.h"

//
#include <iostream>
#include <iomanip>
#include <limits>
#include <cfloat>
#include <time.h>

#include"numerics.h"

//
namespace Math_Group
{
 using namespace std;
/**************************************************************************
Task: Linear equation::Constructor
Programing:
10/2007 WW/
**************************************************************************/
Linear_EQS::Linear_EQS(const SparseTable &sparse_table, 
                       const long dof, bool messg):M(NULL),
                       pivots(NULL), z_vec(NULL), message(messg)
{
  long i, size;
  A = new SparseMatrix(sparse_table, dof);
  size = A->Dim();
  size_global = 0;

  x = new double[size];
  b = new double[size];
  //
  for(i=0; i<size; i++)
  {
    b[i] = 0.;
  }
  iter = 0;
  bNorm = 1.0;
  error = 1.0e10;  

}
/**************************************************************************
Task: Linear equation::Destructor
Programing:
10/2007 WW/
**************************************************************************/
Linear_EQS::~Linear_EQS()
{
  if(A)
  { 
    delete A;
  }
  if(x) delete [] x;
  if(b) delete [] b;
  if(pivots) delete [] pivots;
  if(z_vec) delete [] z_vec;
  //
  A = NULL;
  x = NULL;
  b = NULL;
  pivots = NULL;
  z_vec = NULL;
  //

  ///GMRES. 30.06.2010. WW
  if(solver_type == 13)
     H.ReleaseMemory();
}
/**************************************************************************
Task: Linear equation::
Programing:
10/2007 WW/
**************************************************************************/
void Linear_EQS::ConfigNumerics(Numerics *num, const long n)
{
  //
  int i, nbuffer = 0; // Number of temperary float arrays
  long size;
  size = A->Dim();
  
  /// Just to avoid compilation warning
  i = n; 
 
  /// PDE related
  precond_type = num->GetPrecType();
  solver_type =  num->GetType();
  max_iter = num->GetMax_Iteration();
  tol = num->GetTolerance();

  switch(solver_type)
  {
    case 1:
      solver_name = "Gauss"; 
      break;
    case 2:
      solver_name = "BiCGSTab"; 
      nbuffer = 8;
      break;
    case 3:
      solver_name = "BiCG";
      nbuffer = 8;
      break;
    case 4:
      solver_name = "QMRCGStab"; 
      break;
    case 5:
      solver_name = "CG"; 
      nbuffer = 3;
      break;
    case 6:
      solver_name = "CGNR"; 
      break;
    case 7:
      solver_name = "CGS"; 
      nbuffer = 9;
      break;
    case 8:
      solver_name = "Richardson"; 
      break;
    case 9:
      solver_name = "JOR";
      break;
    case 10:
      solver_name = "SOR";
      break;
    case 11:
      solver_name = "AMG1R5";
      break;
    case 12:
      solver_name = "UMF"; 
    case 13:   // 06.2010. WW
      solver_name = "GMRES"; 
      m_gmres = num->GetSub_Dim();
      for(i=0; i<4; i++) 
      {
         double *new_array = new double[m_gmres+1];
         f_buffer.push_back(new_array);
      }  
      H.resize(m_gmres+1, m_gmres+1);
      nbuffer = m_gmres+4;
      break;
  }

  for(i=0; i<nbuffer; i++) 
  {
    double *new_array = new double[size];
    f_buffer.push_back(new_array);
  }  

  //---------------------------------------------
  switch(precond_type)
  {
    case 1:
      precond_name = "Jacobi";
      break;
    case 100:
      precond_name = "ILU not available. Use Jacobi"; 
      precond_type = 1;
      break;
    default:
      precond_name = "No preconditioner"; break;
  }
  //  
  //
  ///TODO max_iter = m_num->ls_max_iterations;
  ///TODO tol = m_num->ls_error_tolerance;
  //
}
/**************************************************************************
Task: Linear equation::Alocate memory for solver
Programing:
11/2007 WW/
**************************************************************************/
void Linear_EQS::Initialize()
{
  long i, size_A;

  (*A) = 0.;
  size_A = A->Dim();
  for(i=0; i<size_A; i++)
    b[i] = 0.;
  error = 1.0e10;                     
}
/**************************************************************************
Task: Linear equation::Alocate memory for solver
Programing:
11/2007 WW/
**************************************************************************/
void Linear_EQS::Clean()
{

  for(int i=0; i<(int)f_buffer.size(); i++)
  {
    if(f_buffer[i]) delete [] f_buffer[i];
    f_buffer[i] = NULL;
  }
  f_buffer.clear();
  if(pivots) delete [] pivots; 
  if(z_vec) delete [] z_vec; 
  if(M) delete M;
  M = NULL; 
  pivots = NULL;
  z_vec = NULL;
}
/**************************************************************************
Task: Linear equation::Write
Programing:
11/2007 WW/
**************************************************************************/
void Linear_EQS::Write(ostream &os)   
{
  long i, size_A;
  A->Write(os);
  size_A = A->Dim();
 
  //
  os<<" b ( RHS): " <<endl;  
  os.width(14);
  os.precision(8); 
  //
  for(i=0; i<size_A; i++)
    os<<setw(10)<<i<<" "<<setw(15)<<b[i]<<endl;
}
//**************************************************************************
/*!
    \brief Write the equation into a binary file
     
     Programing:
     03/2011 WW
*/
//**************************************************************************
void Linear_EQS::Write_BIN(ostream &os)   
{
 
  if((A->GetStorageType() != CRS )||(!A))
    return;
   
  A->Write_BIN(os);  
  os.write((char*) b, A->Dim()*sizeof(double)); 
}

/**************************************************************************
Task: Linear equation::Solver
Programing:
11/2007 WW/
**************************************************************************/
int Linear_EQS::Solver()
{
  //
  double eqs_time = -clock(); 

  iter = 0;
  ComputePreconditioner();
  switch(solver_type)
  {
    case 1:
      iter =  Gauss();
      break;
    case 2:
      iter = BiCGStab();
      break;
    case 3:
      iter =  BiCG();
      break;
    case 4:
      iter = QMRCGStab();
      break;
    case 5:
      iter =  CG();
      break; 
    case 6:
      iter = CGNR();
      break;
    case 7:
      iter = CGS();
      break;
    case 8:
      iter = Richardson();
      break;
    case 9:
      iter = JOR();
      break;
    case 10:
      iter = SOR();
      break;
    case 11:
      iter = AMG1R5();
      break;
    case 12:
      iter = UMF();
      break;
    case 13:
       iter =  GMRES();
      break;
  }

  eqs_time += clock ( );
  cout<<"Clock time in linear solver: "<< eqs_time/(CLOCKS_PER_SEC)<<"\n";

  return iter; 
}

// Preconditioners
/**************************************************************************
Task: Preconditioners
Programing:
10/2007 WW
**************************************************************************/
void Linear_EQS::ComputePreconditioner()
{
  switch(precond_type)
  {
    case 1:
      return;
    case 100:
      return;
    default:
      return;
  }
}
/**************************************************************************
Task: Linear equation::SetKnownXi
      Configure equation system when one entry of the vector of 
      unknown is given 
Programing:
10/2007 WW/
**************************************************************************/
void Linear_EQS::SetKnownX_i(const long i, const double x_i)
{
   A->Diagonize(i, x_i, b);  
}
/**************************************************************************
Task: Linear equation::Preconditioner
Programing:
08/2007 WW/
**************************************************************************/
void Linear_EQS::Precond(double *vec_s, double *vec_r)
{ 
  bool pre = true;
  switch(precond_type)
  {
    case 1:
      A->Precond_Jacobi(vec_s, vec_r);
      break;
    case 100:
      A->Precond_Jacobi(vec_s, vec_r);
      break;
    default:
      pre = false; //A->Precond_ILU(vec_s, vec_r);
      break; 
  }
  if(!pre) 
  {
     for(long i=0; i<A->Dim(); i++)
       vec_r[i] = vec_s[i];
  } 
}
/**************************************************************************
Task: Linear equation:: M^T x
  Transpose of preconditioner times a vector
Programing:
08/2008 WW/
**************************************************************************/
void Linear_EQS::TransPrecond(double *vec_s, double *vec_r)
{
   Precond(vec_s, vec_r);
}
/*\!
********************************************************************
   Dot production of two vectors
   Programm:  
   10/2007 WW
   12/2007 WW  Parallel
********************************************************************/
double Linear_EQS::dot (const double *xx,  const double *yy)
{
   long i;
   double val = 0.; 
   long size = 0;
   if(A)
     size = A->Dim(); 
   else  /// For JFNK. 01.09.2010. WW
     size = size_global; 
   for(i=0; i<size; i++)
     val += xx[i]*yy[i];
   return val;
}
/*\!
********************************************************************
   Dot production of two vectors
   Programm:  
   01/2008 WW  
********************************************************************/
double  Linear_EQS::NormX() 
{ 
  return sqrt(dot(x, x )); 
} 
//
/*\!
********************************************************************
   ConvergeTest
   Programm:  
   09/2007 WW
********************************************************************/
void Linear_EQS::Message()
{

  if (!message) return;
  cout.width(10);
  cout.precision(3);
  cout.setf(ios::scientific);
  //
  //system("color 0B");
  cout<<"\n================================================\n";         
  cout<<"Linear solver "<<solver_name<<" with "<<precond_name<<":\n";         
  //cout<<"\n------------------------------------------------ \n";         
  cout<<"Iterations |";         
  cout<<"Max Iters  |";         
  cout<<"Norm of b  |"; 
  cout<<"Error      |\n"; 
  //cout<<"\n------------------------------------------------ \n";         
  cout<<setw(11)<<iter<<"|"<<setw(11)<<max_iter<< "|"
      <<setw(11)<<bNorm<< "|"<<setw(11)<<error<<"|\n";          
  cout<<"================================================\n";         
  if (iter==max_iter) cout << " Maximum iteration reached !!! \n";
    cout.flush();    
//
}
/*\!
********************************************************************
   Check if the norm of b is samll enough for convengence. 
   normb_new is given to bNorm;
   Programm:  
   09/2007 WW
********************************************************************/
inline bool Linear_EQS::CheckNormRHS(const double normb_new)
{
  if(bNorm>DBL_EPSILON)
  {
    if((normb_new/bNorm)<tol)
    {
       error = normb_new/bNorm;
       bNorm = normb_new;
       Message();
       return true; 
    }    
  }
  bNorm = normb_new;
  if(bNorm<DBL_MIN)
  {
     error = 0.;
     Message();
     return true; 
  }
  return false;
}

/**************************************************************************
Task: Linear equation::CG
Programing:
11/2007 WW/
**************************************************************************/
int Linear_EQS::CG()
{
  long i, size;
  double rrM1;
  double *p, *r, *s;
  //
  size = A->Dim();
  p = f_buffer[0];
  r = f_buffer[1];
  s = f_buffer[2];
  //
  double bNorm_new = Norm(b);
  // Check if the norm of b is samll enough for convengence
  if(CheckNormRHS(bNorm_new))
    return 0;
  //
  // r0 = b-Ax 
  A->multiVec(x,s); 
  for(i=0; i<size; i++)
    r[i] = b[i]-s[i];	
  //
  // Preconditioning: M^{-1}r
  Precond(r, s);
  for(i=0; i<size; i++)
    p[i] = s[i];	
  // Check the convergence
  if ((error = Norm(r)/bNorm) < tol)
  {
    Message();
    return 1; 
  }
  //
  double rr = dot(r , s);
  //
  for (iter=1; iter<=max_iter; ++iter)
  {
    A->multiVec(p, s);
    const double alpha = rr / dot(p,s);
    // Update 
    for(i=0; i<size; i++)
    {
      x[i] += alpha*p[i];
      r[i] -= alpha*s[i];
    }
    if ((error = Norm(r)/bNorm) < tol)
    {
      Message();
      return iter <= max_iter; 
    }
    //
    Precond(r, s);
    //
    rrM1 = rr;
    rr   = dot(s, r);
    const double beta = rr / rrM1;
    for(i=0; i<size; i++)
      p[i] = s[i] + beta*p[i];
  }
  //
  Message(); 
  return iter <= max_iter;
}
/**************************************************************************
Task: Linear equation::BiCG
Programing:
08/2008 WW/
**************************************************************************/
int Linear_EQS::BiCG()
{
  long i, size;
  double rho1, rho2 = 0., alpha, beta;
  double *z, *zt, *p, *pt, *q, *qt, *r, *rt;
  //
  size = A->Dim();
  z = f_buffer[0];
  zt = f_buffer[1];
  p = f_buffer[2];
  pt = f_buffer[3];
  q = f_buffer[4];
  qt = f_buffer[5];
  r = f_buffer[6];
  rt = f_buffer[7];

#ifdef CG_test
 ofstream Dum("BiCG.txt", ios::out);
 Dum.width(20);
 Dum.precision(15);
 Dum.setf(ios::scientific);
#endif
  //
  double bNorm_new = Norm(b);
  // Check if the norm of b is samll enough for convengence
  if(CheckNormRHS(bNorm_new))
    return 0;
  //
  // r0 = b-Ax 
  A->multiVec(x,rt); 
  for(i=0; i<size; i++)
  {
    r[i] = b[i]-rt[i];	
    rt[i] = r[i];	
  }
  //
  // Check the convergence
  if ((error = Norm(r)/bNorm) < tol)
  {
    Message();
    return 1; 
  }

#ifdef CG_test
  Dum<<" Norm(r) "<< Norm(r)<<" bNorm "<<bNorm<<endl;
#endif


  //
  //
  for (iter=1; iter<=max_iter; ++iter)
  {
    Precond(r, z);   
    TransPrecond(rt, zt);
    rho1 = dot(z, rt);

#ifdef CG_test
  Dum<<" rho1 "<< rho1<<endl;
#endif

    //
    if (fabs(rho1) < DBL_MIN)
    {
       Message();
       return iter <= max_iter; 
    }
    //
    if(iter == 1)
    {
       for(i=0; i<size; i++)
       {
          p[i] = z[i];	
          pt[i] = zt[i];	
       }
    }
    else
    {
       beta = rho1/rho2;
       for(i=0; i<size; i++)
       {
          p[i] = z[i] + beta*p[i];	
          pt[i] = zt[i] + beta*pt[i];	
       }
    }
    // 
    A->multiVec(p, q);
    A->Trans_MultiVec(pt, qt);
    alpha = rho1/dot(pt,q);

#ifdef CG_test
  Dum<<" alpha "<< alpha<<endl;
 
#endif

    // 
    for(i=0; i<size; i++)
    {
       x[i] += alpha*p[i];	
       r[i] -= alpha*q[i];	
       rt[i] -= alpha*qt[i];	
    }
    //
    rho2 = rho1;
    if ((error = Norm(r)/bNorm) < tol)
    {
      Message();
      return iter <= max_iter; 
    }

#ifdef CG_test
  Dum<<" error = Norm(r)/bNorm "<< error<<endl;
#endif


    //
  }
  //
  Message(); 
#ifdef CG_test
  
  exit(0);
#endif

  return iter <= max_iter;
}

/*************************************************************************
GeoSys-Function:
Task: BiCGStab solver 
Programming: 
10/2007 WW  
**************************************************************************/
int Linear_EQS::BiCGStab() 
{
  long i, size;
  double rho_0, rho_1, alpha, beta, omega, tt=0., norm_r = 0.;
  double *r0, *r, *s, *s_h, *t, *v, *p, *p_h;
  // 
  if(A)
    size = A->Dim();
  else   //JFNK. 26.11.2010
    size = size_global;

  r0 = f_buffer[0];
  r = f_buffer[1];
  s = f_buffer[2];
  s_h = f_buffer[3];
  t = f_buffer[4];
  v = f_buffer[5];
  p = f_buffer[6];
  p_h = f_buffer[7];
  //
  rho_0 = alpha = omega = 1.0;
  //
  double bNorm_new = Norm(b);
  // Check if the norm of b is samll enough for convengence
  if(CheckNormRHS(bNorm_new))
    return 0;
  //

  //Norm of M r
  A->multiVec(x,s);   // s as buffer



  for(i=0; i<size; i++)
     r0[i] = b[i]-s[i];	  // r = b-Ax

  for(i=0; i<size; i++)
  {
     r[i] = r0[i];	
     v[i] = 0.;	
     p[i] = 0.;	
  }


  if ((error = Norm(r)/bNorm) < tol)
  {
    Message();
    return 0; 
  }
  //  
  for (iter = 1; iter <= max_iter; iter++)
  {
    rho_1 = dot(r0, r);
    if (fabs(rho_1) < DBL_MIN) // DBL_EPSILON
    {
      Message();
      return 0;
    }
    if (iter == 1)
    { 

      for(i=0; i<size; i++)
        p[i] = r[i];
    }
    else
    {
      beta = (rho_1/rho_0) * (alpha/omega);

      for(i=0; i<size; i++)
        p[i] = r[i] + beta * (p[i] - omega * v[i]);
    }
    // Preconditioner
    Precond(p, p_h);  
    A->multiVec(p_h, v);
    //
    alpha = rho_1 / dot(r0, v);
    //

    for(i=0; i<size; i++)
      s[i] = r[i] - alpha * v[i];
    if ((error = Norm(s)/bNorm) < tol) 
    {

      for(i=0; i<size; i++)
        x[i] += alpha * p_h[i];     
      Message();
      return iter <= max_iter;
    }
    //  M^{-1}s, 
    Precond(s, s_h); 
    // A* M^{-1}s
    A->multiVec(s_h, t);
    //
    tt = dot(t,t);
    if(tt>DBL_MIN)
      omega = dot(t,s) / tt;
    else
      omega = 1.0;
    // Update solution
    for(i=0; i<size; i++)
    {
      x[i] += alpha * p_h[i] + omega * s_h[i];
      r[i] = s[i] - omega * t[i];
    }


    rho_0 = rho_1;
    //  
    norm_r = Norm(r);
    if ((error = norm_r / bNorm) < tol)
    {
      Message();
      return iter <= max_iter;
    }
    if (fabs(omega) < DBL_MIN)
    {
      error = norm_r / bNorm;
      Message(); 
      return iter <= max_iter;
    }
  }
  //
  Message(); 
  //
  return iter <= max_iter;
}

/*************************************************************************
GeoSys-Function:
Task: CGS solver 
Programming: 
11/2007 WW  
**************************************************************************/
int Linear_EQS::CGS() 
{
  long i, size;
  double rho_1, rho_2, alpha, beta;
  double *r0, *r, *p, *p_h, *q, *q_h, *v, *u, *u_h;
  // 
  size = A->Dim();
  r0 = f_buffer[0];
  r = f_buffer[1];
  p = f_buffer[2];
  p_h = f_buffer[3];
  q = f_buffer[4];
  q_h = f_buffer[5];
  v = f_buffer[6];
  u = f_buffer[7];
  u_h = f_buffer[8];
  //
  rho_1 = rho_2 = 1.0;
  //
  double bNorm_new = Norm(b);
  // Check if the norm of b is samll enough for convengence
  if(CheckNormRHS(bNorm_new))
    return 0;
  //
  A->multiVec(x,v);   // v as buffer
  for(i=0; i<size; i++)
  {
    r0[i] = b[i]-v[i];	  // r = b-Ax
    r[i] = r0[i];	
    v[i] = 0.;	
  }
  if ((error = Norm(r)/bNorm) < tol)
  {
    Message();
    return 0; 
  }
  //  
  for (iter = 1; iter <= max_iter; iter++)
  {
    rho_1 = dot(r0, r);
    if (fabs(rho_1) < DBL_MIN) //  DBL_EPSILON
    {
      Message();
      return 0;
    }
    if (iter == 1)
    { 
      for(i=0; i<size; i++)
        p[i] = u[i] = r[i];
    }
    else
    {
      beta = rho_1/rho_2;
      for(i=0; i<size; i++)
      {
        u[i] = r[i] + beta * q[i];
        p[i] = u[i] + beta * (q[i] + beta * p[i]);
      } 
    }
    // Preconditioner
    Precond(p, p_h);  
    // A M^{-1}p-->v
    A->multiVec(p_h, v);
    //
    alpha = rho_1 / dot(r0, v);
    //
    for(i=0; i<size; i++)
    {
       q[i] = u[i] - alpha * v[i];
       q_h[i] = u[i] + q[i];
    }
     // Preconditioner
    Precond(q_h, u_h);  
    for(i=0; i<size; i++)
      x[i] += alpha * u_h[i];
    //
    A->multiVec(u_h, q_h);
    //
    for(i=0; i<size; i++)
      r[i] -= alpha * q_h[i];
    rho_2 = rho_1;
    if ((error = Norm(r) / bNorm) < tol)
    {
      Message();
      return iter <= max_iter;
    }
  }
  //
  Message(); 
  //
  return iter <= max_iter;
}

#define aGMRES
#ifdef aGMRES

//-----------------------------------------------------------------
/*!
     GMRES solver.
     
     by WW. 06.2010
*/
//-----------------------------------------------------------------
/// For GMRES
inline void Linear_EQS::Get_Plane_Rotation(double &dx, double &dy, double &cs, double &sn)
{  
   if (dy == 0.0) 
   {
       cs = 1.0;
       sn = 0.0;
   } 
   else if (fabs(dy) > fabs(dx)) 
   {
       double temp = dx / dy;
       sn = 1.0 / sqrt( 1.0 + temp*temp );
       cs = temp * sn;
   } 
   else 
   {
       double temp = dy / dx;
       cs = 1.0 / sqrt( 1.0 + temp*temp );
       sn = temp * cs;
   }
}

/// For GMRES. 
inline void Linear_EQS::Set_Plane_Rotation(double &dx, double &dy, double &cs, double &sn)
{
    double temp  =  cs * dx + sn * dy;
    dy = -sn * dx + cs * dy;
    dx = temp;
}

/// Update solution in GMRES
inline void Linear_EQS::Update(double *x, int k, Matrix &h, double *s)
{
  long i;
  long m, j;
  long   size = 0;
  int v_idx0 = 7;
  
  if(A)
    size = A->Dim();
  else
    size = size_global;

  double  *v_j;

  m = m_gmres; 

  double *y = f_buffer[3];
  for (j = 0; j <m+1; j++)
    y[j] = s[j];

  // Back solve 
  for (i = k; i >= 0; i--) {
    y[i] /= h(i,i);
    for ( j = i - 1; j >= 0; j--)
      y[j] -= h(j,i) * y[i];
  }

  for (j = 0; j <= k; j++)
  {
     v_j =  f_buffer[v_idx0+j];
     for(i=0; i<size; i++)
       x[i] += v_j[i] * y[j];
  }  
}
/// GMRES solver. WW
int Linear_EQS::GMRES()
{
 
  long i, k, l, m, size; 
  double normb, beta;
  double *s, *cs, *sn,  *v, *w, *r, *t, *v_k;

  normb = Norm(b);
  // Check if the norm of b is samll enough for convengence
  if(CheckNormRHS(normb))
    return 0;
  
  //
  m = m_gmres; 
  if(A)
    size = A->Dim();
  else
    size = size_global;
  s = f_buffer[0];
  cs = f_buffer[1];
  sn = f_buffer[2];
  w = f_buffer[4];
  r = f_buffer[5];
  t = f_buffer[6]; // Buffer array

  int v_idx0 = 7;

  // Norm of Mb
  Precond(b, r);
  // Here Mb-->r  
  normb = Norm(r);

  //Norm of M r
  A->multiVec(x,w);   // Ax-->w
  for(l=0; l<size; l++)
    r[l] = b[l]-w[l];	  // r = b-Ax.

  Precond(r, w);         // Mr-->w
  beta = Norm(w);
  
  if (normb < DBL_MIN)
     normb = 1;
  
  //if ((error = Norm(r) / normb) <= tol) 
  if ((error = beta / normb) <= tol) 
  {
      Message();
      return 0; 
  }

  iter = 1;
  while (iter <= max_iter) 
  {
     v =  f_buffer[v_idx0];
     for(l=0; l<size; l++)
        v[l] = r[l] / beta;    //  r/beta
     for(l=0; l<m+1; l++)
        s[l] = 0.0;
     s[0] = beta;

     for (i = 0; i < m && iter <= max_iter; i++, iter++)
     {
        v =  f_buffer[v_idx0+i];
        A->multiVec(v, t);
        Precond(t, w); 

        for (k = 0; k <= i; k++) 
        {
            v_k = f_buffer[v_idx0+k];
            H(k, i) = dot(w, v_k);
       
            for(l=0; l<size; l++)
               w[l] -= H(k, i) * v_k[l];
        }
        H(i+1, i) = Norm(w);
        v_k = f_buffer[v_idx0+i+1];
        for(l=0; l<size; l++)      
           v_k[l] = w[l] / H(i+1, i); 

        for (k = 0; k < i; k++)
          Set_Plane_Rotation(H(k,i), H(k+1,i), cs[k], sn[k]);
      
        Get_Plane_Rotation(H(i,i), H(i+1,i), cs[i], sn[i]);
        Set_Plane_Rotation(H(i,i), H(i+1,i), cs[i], sn[i]);
        Set_Plane_Rotation(s[i], s[i+1], cs[i], sn[i]);
      
        if ((error = fabs(s[i+1]) / normb) < tol)
        {
            Update(x, i, H, s);
            Message();
            return iter <= max_iter;
        }
     }
  
     Update(x, i - 1, H, s);
     A->multiVec(x, t);
  
     for(l=0; l<size; l++)
        w[l] = b[l]-t[l];	  // r = b-Ax.
     Precond(w, r);    // M*r 

     beta = Norm(r);
     if ((error = beta / normb) < tol)
     {
        Message();
        return iter <= max_iter;
     }
  }
  
  Message();
  return iter <= max_iter;
}
//-----------------------------------------------------------------
#endif //GMRES



/*************************************************************************
GeoSys-Function:
Task: Parallel preconditioner, inverse
Programming: 
12/2007 WW  
02/2011 WW
**************************************************************************/
void Linear_EQS::Precond_Jacobi(const double *vec_s, double *vec_r)
{
   double val;
   long size;
   size = A->Dim();
   
   for(long i=0; i<size; i++)
   {
      val = pivots[i];
      //  <DBL_EPSILON
      if(fabs(val)<DBL_MIN)
        val = 1.0;
      vec_r[i] = vec_s[i]/val;
   }

}

//------------------------------------------------------------------------
} // namespace

