/*!
  \file geo.h
   
   Declaration of class Point and Polyline, which are used
   to assign  boundary conditions
   
   13.04.2011. WW
*/
#ifndef geo_INC
#define geo_INC

#include "misc.h"
#include<iostream>
#include<fstream>

namespace _FDM
{
   
   class FiniteDifference;
   class BoundayCondition; 
   /*!
     \class  Point
         
   */
   enum Geo_Point_Type {none, neumann, dirichl};
   class Point
   {
      public:
        Point(long id, real x, real y) 
        {
           index = id;
           coordinates = new real[2];
           coordinates[0] = x;
           coordinates[1] = y;
           point_type = none;
        }
        
        ~Point()
        {
           delete [] coordinates;            
        }  

        real X() const {return coordinates[0];}
        real Y() const {return coordinates[1];}
 
        void Write(ostream &os = cout);
        long Index() const {return index;}
      private:
        long index;
        long grid_i;
        long grid_j;

        real *coordinates; 
        Geo_Point_Type point_type;
        

        vector<Point*> neighbors; 

        friend class Polyline;
        friend class FiniteDifference;
        friend class BoundayCondition; 

   };

   /*!
      \class Polyline
          
      Define a polyine that consists of points
   */
    class Polyline
    {
        public:
          Polyline(ifstream &ins, string ply_name);
          ~Polyline()
            {
               points.clear();
            }

        string Name() const {return name;}

        void Write(ostream &os = cout);
        
        bool PointInDomain(double x, double y);

        private:
          vector<Point*> points;
          string name;
          friend class Polyline;
          friend class FiniteDifference;
          friend class BoundayCondition;           
    };
   
}

/// Contains all points of the domain
extern vector<_FDM::Point*> points;
extern vector<_FDM::Polyline*> polylines;
extern void GeoRead(); 
extern _FDM::Polyline *GetPolylineByName(string name); 
extern _FDM::Point *GetPointByID(long ID); 
extern void WriteGeoData(ostream &os = cout);



#endif