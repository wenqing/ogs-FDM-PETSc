#ifndef mat_INC
#define mat_INC
#include<iostream>
#include<fstream>

#include"misc.h"
using namespace std;
namespace _FDM
{
   class Mat_Property
   {
       public:
         Mat_Property(ifstream &ins);
         ~Mat_Property() {}

         void Write(ostream &os = cout);

       private:
         real conductivity;
         real storage;      
   };
   
}
#endif