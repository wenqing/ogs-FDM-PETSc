/*!
  \file Polyline.cpp

   Definition of member fucntions of class Polyline, and the function 
   used to read geometrical data

   14.04.2011. WW   
*/

#include "Geometry.h"

#include <sstream>
#include <limits>
#include <cfloat>
#include <cstdlib>

#include "Numerics.h"
#include "misc.h"

#include "Point.h"
#include "Polyline.h"
#include "GeoEntity.h"


namespace Geometry_Group
{

using namespace std;
using namespace AuxFunctions;
/*!
  \fn  ReadPolyline()
  Read geometrical data

  14.04.2011. WW
*/
void Geometry::GeoRead(string file_name)
{   
   long id; 
   float xy[2];
   string aline;
   std::stringstream ss;

   string geo_fname = file_name+".geo"; 
   ifstream ins(geo_fname.c_str());

   if(!ins.good()) 
   {
      cout<<"Could not find file "<<geo_fname<<". Stop now!"<<endl;
      exit(1);
   } 

   cout<<">> Read geometry data."<<endl;
   while(!ins.eof())
   {
      getline(ins, aline); 
      if(CheckComment(aline))
          continue;
	  aline = AuxFunctions::string_To_lower(aline);
      if(aline.find("point")!=string::npos)  
      {   
         for(;;)
         {
            getline(ins, aline); 
            if(aline.find("...")!=string::npos)
              break; 
         
            ss.str(aline);
            ss>> id >> xy[0] >> xy[1];
            ss.clear();

            points.push_back(new Geometry_Group::Point(id, xy[0], xy[1])); 
         }
      }

      if(aline.find("polyline")!=string::npos) 
      {
         ss.str(aline);
         // Skip "---" and "polyline"
         ss>>aline>>aline>>aline;
         ss.clear();
         polylines.push_back(new Geometry_Group::Polyline(ins, aline, this)); 

      }
       

   }

} 

/*!
    \fn Polyline *getPolylineByName(string name); 
   
     Find a polyline by name 
      
     04.2011. WW 
*/
Polyline *Geometry::getPolylineByName(string name)
{
   int i;
  
   for(i=0; i<(int)polylines.size(); i++)
   {
      if(polylines[i]->Name().find(name)!=string::npos)
      {
         return polylines[i]; 
      }
   }          
   return NULL;
}
/*!
    \fn FDM::Point *getPointByID(long ID); 
   
     Find a polyline by name 
      
     04.2011. WW 
*/

Point *Geometry::getPointByID(long ID)
{
   int i;
  
   for(i=0; i<(int)points.size(); i++)
   {
      if(points[i]->Index() == ID)
      {
         return points[i]; 
      }
   }          
   return NULL;
}
//-------------------------------------------
/*!
   \fn  WriteGeoData(ostream &os = cout);
   
    Output Geometrical Data
*/
void Geometry::WriteGeoData(ostream &os)
{
   size_t i;

   os<<"--- point"<<endl;
   for(i=0; i<points.size(); i++)
     points[i]->Write(os);
   os<<"...\n"<<endl;

  
   for(i=0; i<polylines.size(); i++)
      polylines[i]->Write(os);

}

void Geometry::GeoReleaseMemory()
{
   DeleteVector(points);
   DeleteVector(polylines);
}

}//Namespace

