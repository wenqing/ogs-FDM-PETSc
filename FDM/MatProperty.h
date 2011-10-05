#ifndef Mat_Property_INC
#define Mat_Property_INC
#include<iostream>
#include<fstream>

#include"misc.h"
using namespace std;
namespace _FDM
{
   class FiniteDifference;

   class Mat_Property
   {
       public:
         Mat_Property(ifstream &ins);
         ~Mat_Property() {}

         void Write(ostream &os = cout);

       private:
         real conductivity;
         real storage;  

         friend class FiniteDifference;     
   };
   
}
#endif