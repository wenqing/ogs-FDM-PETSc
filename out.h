/*!
  \brief
   Class to handle output. 

   05.2011 WW
*/
#ifndef out_INC
#define out_INC
#include<iostream>
#include<fstream>

#include"misc.h"
using namespace std;
namespace _FDM
{
   class FiniteDifference;

   class Output
   {
       public:
         Output(ifstream &ins);
         ~Output();

         void Write(ostream &os = cout);

       private:
         vector<float> at_times;
         int steps;

         string fname;
         ofstream *os; 
         /// class Geo_Root;  

         friend class FiniteDifference;     
   };
   
}
#endif