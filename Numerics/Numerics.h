#ifndef numerics_INC
#define numerics_INC
#include<iostream>

#include "misc.h"

namespace _FDM
{
   class Numerics
   {
       public:
         Numerics(std::ifstream &ins);
         ~Numerics() {}

		 void Write(std::ostream &os = std::cout);

         int GetType() const {return type;}
         int GetPrecType() const {return prec_type;}
         int GetSub_Dim() const {return sub_dim;}
         int GetMax_Iteration() const {return max_ite;}
         double GetTolerance() const {return tol;}

       private:

         int type;
         int prec_type;
         int max_ite;
         int sub_dim; /// For GMRES
         double tol; 
          

         //For output;
         std::string name;
         std::string prec_name; 
   };
}

extern double ComputeDetTri(const double *x1, const double *x2,
                                const double *x3);
#endif