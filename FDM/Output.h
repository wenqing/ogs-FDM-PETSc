/*!
  \brief
   Class to handle output. 

   05.2011 WW
*/
#ifndef Output_INC
#define Output_INC
#include<iostream>
#include<fstream>

#include"misc.h"
using namespace std;
namespace _FDM
{
   class FiniteDifference;
   class Geo_Entity;

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
        
         /// geomtry;
         Geo_Entity *geo_entity;  

         friend class FiniteDifference;     
   };
   
}
#endif