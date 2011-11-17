/*!
  \file Point.cpp

   Definition of member fucntions of class Point

   14.04.2011. WW   
*/

#include "Point.h"


#include <vector>
#include <sstream>
#include <limits>
#include <cfloat>

#include "Numerics.h"
#include "misc.h"

#include "GeoEntity.h"


namespace Geometry_Group
{
   using namespace std;
   /*!
      \fn Point::Write(ostream &os)
        
       output the coordinate of a point
   */
   void Point::Write(ostream &os)
   {
      os<<index<<" "<< coordinates[0]<<" "<<coordinates[1]<<"  0."<<endl;      
   }
   void Point::Write_VTK(ostream &os)
   {
      os<< coordinates[0]<<" "<<coordinates[1]<<"  0."<<endl;      
   }

}


