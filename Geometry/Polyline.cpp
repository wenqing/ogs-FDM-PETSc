/*!
  \file Polyline.cpp

   Definition of member fucntions of class Polyline, and the function 
   used to read geometrical data

   14.04.2011. WW   
*/

#include "Polyline.h"

#include <vector>
#include <sstream>
#include <limits>
#include <cfloat>

#include "Numerics.h"
#include "misc.h"

#include "Point.h"
#include "GeoEntity.h"
#include "Geometry.h"


namespace Geometry_Group
{

   using namespace std;
   using namespace AuxFunctions;
   /*!
      \fn constructor of class Polyline
      
   */
   Polyline::Polyline(ifstream &ins, string ply_name, Geometry *geo): Geo_Entity(ply), name(ply_name)
   {
      long id; 
      string aline;
      std::stringstream ss;
  
      for(;;)
      {
         getline(ins, aline); 
         if(aline.find("...")!=string::npos)
           break; 
         
         ss.str(aline);
         ss>> id ;
         ss.clear();

         points.push_back(geo->getPointByID(id)); 
      }
         
   }
    
   /// Destructor
   Polyline::~Polyline()
   {
      points.clear();
   }

   /*!
      \fn Point::Write(ostream &os)
        
       output the coordinate of a polyline
   */
   void Polyline::Write(ostream &os)
   {
      os<<"--- polyline "<<name<<endl;
      for(int i=0; i<(int)points.size(); i++)
        os<<points[i]->Index()<<endl;  
      os<<"...\n"<<endl;    
   }

   /*!
     \fn  PointInDomain(double x, double y)

     Determine whether a point is in the domain 

     04.2011. WW
   */ 
   bool Polyline::PointInDomain(double x, double y)
   {
      Point *pnt_i, *pnt_j; 
      real xi, yi, xj, yj, dx, dy, ku, kv;
      int size = (int)points.size();
      int   i, j;
      j = size-1 ;
      bool  inside = false;

	  
      for (i=0; i<size; i++) 
      {
         pnt_i = points[i];
         pnt_j = points[j];
         xi =  pnt_i->X();
         yi =  pnt_i->Y();
         xj =  pnt_j->X();
         yj =  pnt_j->Y();

         dx = x-xi;
         dy = y-yi;
         if(sqrt(dx*dx+dy*dy)<DBL_EPSILON)
		 {
		    inside = true;
		    break; 	 
		 }

         if(fabs(xj-xi)>DBL_MIN)               
         {
             ku = (yj-yi)/(xj-xi);
             /// p, pi, pj on the same line  
             if(fabs(x-xi ) > DBL_MIN)
			 {
                kv = (y-yi)/(x-xi);
				/// p_pi and pj_pi are parallel
                if(fabs(ku-kv) < DBL_MIN)
				{
                   kv = (x-xi)/(xj-xi);
                   if(kv>=0.0&&kv<=1.0) // p is inbetween pi and pj
				   {
                      inside = true;
                      break;					  
				   }
				}

			 }
		 }
		 else // Vertical line: xi = xj
		 {
             if(fabs(x-xi)<DBL_MIN)    
			 {
                kv = (y-yi)/(yj-yi);
                if(kv>=0.0&&kv<=1.0) // p is inbetween pi and pj
			    {
                    inside = true;
                    break;					  
				}
			 }
		 }

         if (( (yi<y) && (yj>=y) )||( (yj<y) && (yi>=y)) ) 
         {
            if (xi+((y-yi)/(yj-yi))*(xj-xi)<x)
            {
              inside = !inside; 
            }
         }
         j=i;
      }
	  
       return inside;  
   }
   /*!
     \fn real MinDisttanceTo_a_Point(Point *pnt);

      Calculate the distance to a point 
      
      05.2011 WW   
   */ 
   real Polyline::MinDisttanceTo_a_Point(const Point *pnt)
   {
      size_t i, k;
      real a[3], b[3], c[3], dist; 
      real min_dist = DBL_MAX;
     
      for(i=0; i<points.size()-1; i++) 
      {
         for(k=0; k<3; k++) 
         {
            a[k] = pnt->coordinates[k];
            b[k] = points[i]->coordinates[k];
            c[k] = points[i+1]->coordinates[k];
         }
         dist = 2.0*ComputeDetTri(a,b,c)/points[i]->getDistanceTo(points[i+1]);
         if(dist<min_dist)
           min_dist = dist;
      }
      return min_dist; 
   } 
}


