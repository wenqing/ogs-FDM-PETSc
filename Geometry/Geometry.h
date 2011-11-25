/*!
  \file Polyline.h
   
   Declaration of fundtions, which are used to handel the geometrical data
   
   05.09.2011. WW
*/
#ifndef Geometry_INC
#define Geometry_INC

#include<iostream>
#include<fstream>
#include<vector>

namespace Geometry_Group
{
class Point;
class Polyline;

class Geometry
{
   public:
      Geometry() {}
	  ~Geometry() {GeoReleaseMemory(); }


      void GeoRead(std::string file_name); 
      void GeoReleaseMemory(); 
      Polyline *getPolylineByName(std::string name); 
      Point *getPointByID(long ID); 
      void WriteGeoData(std::ostream &os = std::cout);  
  
   private:

      /// Contains all points of the domain
      std::vector<Point*> points;
      std::vector<Polyline*> polylines;
};
}
#endif

