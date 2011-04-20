/*!
   \file File defines class BoundayCondition
  
   04.2011. WW
*/
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "geo.h"
#include "bc.h"

using namespace std;
namespace _FDM
{
   //--------------- class  BoundayCondition ---------------------
   BoundayCondition::BoundayCondition(ifstream &ins)
   {
     long ID = 0;
     string aline;
     std::stringstream ss;

     point = NULL;
     ply = NULL;
     for(int i=0; i<2; i++)
     {
        getline(ins, aline); 
        aline = string_To_lower(aline);
        if(aline.find("geometry")!=string::npos) 
        {
           ss.str(aline);
           /// skip key
           ss>>aline;
           ss>>aline;
           if(aline.find("polyline")!=string::npos)
           {
              /// polyline name
              ss>>aline; 
              ply = GetPolylineByName(aline);             
           }  
           else if(aline.find("point")!=string::npos)
           {
              /// point ID
              ss>>ID; 
              point = GetPointByID(ID);             

           }
           ss.clear();
        } 
        if(aline.find("value")!=string::npos) 
        {
           ss.str(aline);
           /// skip key
           ss>>aline;
           ss>>value;         
           ss.clear();
        }
     }

   }
 
   /// Output boundary condition 
   void BoundayCondition::Write(ostream &os)
   {
      os<<"\t geometry: ";
      if(ply)
         os<<"polyline "<<ply->Name()<<endl;
      else if(point)
         os<<"point "<<point->Index()<<endl;
      os<<"\t value: "<<value<<endl<<endl; 
   }

   /// As the funtion name
   void BoundayCondition::SetGeoEntityType(string type_name)
   {
       Geo_Point_Type geo_pnt_type;
       if(type_name.find("neumann")!=string::npos)
          geo_pnt_type = neumann;
       if(type_name.find("dirichlet")!=string::npos)
          geo_pnt_type = dirichl;
      
       if(point)
          point->point_type = geo_pnt_type;
       if(ply)
       {
          for(int i=0; i<(int)ply->points.size(); i++)
            ply->points[i]->point_type = geo_pnt_type; 
       }   
         
   }

}