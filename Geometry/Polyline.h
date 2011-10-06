/*!
  \file Polyline.h
   
   Declaration of class  Polyline, which is used
   to assign  boundary conditions
   
   13.04.2011. WW
*/
#ifndef Polyline_INC
#define Polyline_INC

#include<iostream>
#include<fstream>
#include<cmath>

#include "misc.h"
#include "GeoEntity.h"



namespace Geometry_Group
{
   
   //using _FDM::FiniteDifference;
   //using _FDM::ConditionDataBC; 
   class Point;
   class Geometry;

   /*!
      \class Polyline
          
      Define a polyine that consists of points
   */
    class Polyline: public Geo_Entity 
    {
        public:
          Polyline(std::ifstream &ins, std::string ply_name, Geometry *geo);
          ~Polyline();

        std::string Name() const {return name;}

        void Write(std::ostream &os = std::cout);
        
        bool PointInDomain(double x, double y);
        real MinDisttanceTo_a_Point(const Point *pnt);

        private:
          std::vector<Point*> points;
          std::string name;
          friend class FiniteDifference;
          friend class ConditionDataBC;           
    };
   
}



#endif