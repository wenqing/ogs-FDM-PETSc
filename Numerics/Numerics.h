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

      int getType() const
      {
         return type;
      }
      int getPrecType() const
      {
         return prec_type;
      }
      int getSub_Dim() const
      {
         return sub_dim;
      }
      int getMax_Iteration() const
      {
         return max_ite;
      }
      double getTolerance() const
      {
         return tol;
      }

      const char *getSolverName() const
      {
         return name.c_str();
      }
      const char *getPreConditionerName() const
      {
         return prec_name.c_str();
      }

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
#endif

