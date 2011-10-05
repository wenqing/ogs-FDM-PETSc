/*!
  \file Polyline.h
   
   Declaration of class  Polyline, which is used
   to assign  boundary conditions
   
   13.04.2011. WW
*/
#ifndef Polyline_INC
#define Polyline_INC

#include "misc.h"
#include<iostream>
#include<fstream>
#include<cmath>

namespace _FDM
{
   
   class FiniteDifference;
   class ConditionData; 

   class Point;
   class Geo_Entity;
   /*!
      \class Polyline
          
      Define a polyine that consists of points
   */
    class Polyline: public Geo_Entity 
    {
        public:
          Polyline(ifstream &ins, string ply_name);
          ~Polyline();

        string Name() const {return name;}

        void Write(ostream &os = cout);
        
        bool PointInDomain(double x, double y);
        real MinDisttanceTo_a_Point(const Point *pnt);

        private:
          vector<Point*> points;
          string name;
          friend class FiniteDifference;
          friend class ConditionData;           
    };
   
}

/// Contains all points of the domain
extern vector<_FDM::Point*> points;
extern vector<_FDM::Polyline*> polylines;
extern void GeoRead(); 
extern void GeoReleaseMemory(); 
extern _FDM::Polyline *GetPolylineByName(string name); 
extern _FDM::Point *GetPointByID(long ID); 
extern void WriteGeoData(ostream &os = cout);



#endif