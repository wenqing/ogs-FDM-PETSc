/*!
   \file File defines class Numerics
  
   04.2011. WW
*/
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>

#include "numerics.h"

using namespace std;
namespace _FDM
{
//--------------- class Numerics ------------------------------
   Numerics::Numerics(ifstream &ins)
   {
     string aline;
     std::stringstream ss;

     for(int i=0; i<4; i++)
     {
        getline(ins, aline); 
        aline = string_To_lower(aline);
        if(aline.find("linear")!=string::npos) 
        {
           ss.str(aline);
           /// skip key
           ss>>aline;
           ss>>aline;
           name = aline;
           if(aline.find("Guass")!=string::npos)
             type = 1; 
           else if(aline.find("BiCGSTab")!=string::npos)
             type = 2; 
           else if(aline.find("BiCG")!=string::npos)
             type = 3; 
           else if(aline.find("CG")!=string::npos)
             type = 5; 
           else if(aline.find("CGS")!=string::npos)
             type = 7; 
           else if(aline.find("GMRES")!=string::npos)
           {
             type = 13; 
             ss>>sub_dim;
           }
           else
             type = 5; 
           ss.clear();
        } 
        if(aline.find("preconditioner")!=string::npos) 
        {
           ss.str(aline);
           /// skip key
           ss>>aline;
           ss>>aline;
           prec_name = aline;
           if(aline.find("jacobi")!=string::npos)
             prec_type = 1; 
           else if(aline.find("ILU")!=string::npos)
             prec_type = 100; 
           else 
             prec_type = -1; 
           ss.clear();
        }
        if(aline.find("tolerance")!=string::npos) 
        {
           ss.str(aline);
           /// skip key
           ss>>aline;
           ss>>tol;
           ss.clear(); 
        }
        if(aline.find("max_iteration")!=string::npos) 
        {
           ss.str(aline);
           /// skip key
           ss>>aline;
           ss>>max_ite;
           ss.clear(); 
        }
     }

   }

      /// Output boundary condition 
   void Numerics::Write(ostream &os)
   {
      os<<"--- Solver"<<endl;  
      os<<"\t linear: "<<name;
      if(type == 13) 
        os<<" "<<sub_dim;   
      os    <<endl;
    
      os<<"\t preconditioer: "<<prec_name<<endl;
      os<<"\t tolerance: "<<tol<<endl;
      os<<"\t max_iteration: "<<tol<<endl<<endl;

   }
}

/*!
     Computer area of a triangle
*/
double ComputeDetTri(const double *x1, const double *x2,
                                const double *x3)
{
    static double u[3], v[3], z[3];
    
    u[0] = x3[0] - x1[0];	
    u[1] = x3[1] - x1[1];
    u[2] = x3[2] - x1[2];

    v[0] = x2[0] - x1[0];	
    v[1] = x2[1] - x1[1];
    v[2] = x2[2] - x1[2];
 
    z[0] = u[1]*v[2] - u[2]*v[1];
    z[1] = u[2]*v[0] - u[0]*v[2];
    z[2] = u[0]*v[1] - u[1]*v[0];

    return 0.5*sqrt(z[0]*z[0]+z[1]*z[1]+z[2]*z[2] );   
} 
