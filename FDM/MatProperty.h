#ifndef Mat_Property_INC
#define Mat_Property_INC
#include<iostream>
#include<fstream>

#include"misc.h"

namespace _FDM
{
   class FiniteDifference;

   class Mat_Property
   {
       public:
        Mat_Property(std::ifstream &ins);
         ~Mat_Property() {}

		 void Write(std::ostream &os = std::cout);

       private:
         real conductivity;
         real storage;  

         friend class FiniteDifference;     
   };
   
}
#endif