#ifndef numerics_INC
#define numerics_INC
#include<iostream>
#include "misc.h"

namespace _FDM
{
   class Numerics
   {
       public:
         Numerics(ifstream &ins);
         ~Numerics() {}

         void Write(ostream &os = cout);

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
         string name;
         string prec_name; 
   };
}

#endif