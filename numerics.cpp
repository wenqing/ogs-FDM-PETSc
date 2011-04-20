/*!
   \file File defines class Numerics
  
   04.2011. WW
*/
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "numerics.h"

using namespace std;
namespace _FDM
{
//--------------- class Numerics ------------------------------
   Numerics::Numerics(ifstream &ins)
   {
     long ID = 0;
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