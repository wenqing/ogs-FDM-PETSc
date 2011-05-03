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
#include<cmath>

namespace _FDM
{
   
   class FiniteDifference;
   class ConditionData; 
   /*!
     \class  Point
         
   */
   enum Point_Type {none, intern, border, nm_11, nm_12, nm_13, nm_14, nm_21, nm_22, nm_23, nm_24};
   enum BC_Type{Neumann, Dirichlet};
   enum NeighborCell_Type {NE, NW, SE, SW};
   enum NeighborPoint_Type {C, E, N, W, S};
  
   /// class Geo_Root;
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
        real GetDistanceTo(const Point *a_p)
           {
              real *a_coord = a_p->coordinates; 
              return sqrt((coordinates[0]-a_coord[0])
                        *(coordinates[0]-a_coord[0])
                      +  (coordinates[1]-a_coord[1])
                        *(coordinates[1]-a_coord[1]));
            }
 
        void Write(ostream &os = cout);
        void Write_VTK(ostream &os = cout);
        long Index() const {return index;}
        long GetNeighborIndex(const int ii) const {return neighbor_points[ii]; } 
        long GetNumNeighborPoints() const {return (long)neighbor_points.size(); } 
        

      private:
        long index;
        long grid_i;
        long grid_j;

        real *coordinates; 
        /// Save Dirichlet or Neumann value here
        real value;
        Point_Type point_type;
        BC_Type bc_type;
        

        vector<long> neighbor_points; 
        vector<NeighborPoint_Type> np_position; 
        vector<NeighborCell_Type> neighbor_cell_type; 

        friend class Polyline;
        friend class FiniteDifference;
        friend class ConditionData; 

   };

   /*!
      \class Polyline
          
      Define a polyine that consists of points
   */
    class Polyline
    {
        public:
          Polyline(ifstream &ins, string ply_name);
          ~Polyline();

        string Name() const {return name;}

        void Write(ostream &os = cout);
        
        bool PointInDomain(double x, double y);

        private:
          vector<Point*> points;
          string name;
          friend class Polyline;
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