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

namespace Geometry_Group {class Geo_Entity; class Geometry;}
namespace _FDM
{
   class FiniteDifference;

   class Output
   {
       public:
         Output(std::string f_path, std::string f_name, std::ifstream &ins, Geometry_Group::Geometry *geometry);
         ~Output();

		 void Write(std::ostream &os = std::cout);

       private:
         std::vector<float> at_times;
         int steps;

         std::string fname;
         std::ofstream *os; 
        
         std::string file_name; 
         std::string file_path; 

         /// geomtry;
         Geometry_Group::Geo_Entity *geo_entity;  

         friend class FiniteDifference;     
   };
   
}
#endif
