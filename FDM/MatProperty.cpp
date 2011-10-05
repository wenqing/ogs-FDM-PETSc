/*!
\file mat.cpp
     
    Contains the definition of class Mat_Property;

    WW. 03.2011  
*/

#include "MatProperty.h"

#include <vector>
#include <sstream>

#include "misc.h"
/*!
   \class class _material

    Handle meterial properties of groundwater flow
     
    by WW. 03.2011
*/
namespace _FDM
{
  Mat_Property::Mat_Property(ifstream &ins)
  {
     vector<string> key(2);
     vector<float> keyval(2);
     
     key[0] = "conductivity:";
     key[1] = "storage:";

     Read_Block(ins, key, keyval);
      
     conductivity = keyval[0];
     storage = keyval[1];      
 }

/*!
  \fn  Mat_Property::Write(ostream &os = cou)
   Write the material parameters

*/ 
  void Mat_Property::Write(ostream &os)
  {
      os<<"--- Material"<<endl;
      os<<"\t conductivity: \t"<<conductivity<<endl;
      os<<"\t storage: \t"<<storage<<endl<<endl; 
  }

}