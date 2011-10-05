/*!
   \file File defines class ConditionData
  
   04.2011. WW
*/
#include <iomanip>
#include <limits>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "geo.h"
#include "bc.h"

using namespace std;
namespace _FDM
{
   //--------------- class  ConditionData ---------------------
   ConditionData::ConditionData(ifstream &ins)
   {        

     long ID = 0;
     string aline;
     std::stringstream ss;

     point = NULL;
     ply = NULL;
     for(int i=0; i<2; i++)
     {
        getline(ins, aline); 
        if(CheckComment(aline))
          continue;


        aline = string_To_lower(aline);
        if(aline.find("geometry")!=string::npos) 
        {
           ss.str(aline);
           if(aline.find("polyline")!=string::npos)
           {
              /// skip key
              ss>>aline;
              ss>>aline;
              /// polyline name
              ss>>aline; 
              ply = GetPolylineByName(aline);             
           }  
           else if(aline.find("point")!=string::npos)
           {
              /// skip key
              ss>>aline;
              ss>>aline;
              /// point ID
              ss>>ID; 
              point = GetPointByID(ID);             

           }
           //else if(aline.find("domain")!=string::npos)
           else
           {
              cout<<"Geometry type is nit specified. Please check .dat file"<<endl;
              exit(1);
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
 
   //--------------------------------------------------
   /*!
   \fn  Point* ConditionData::GetClosedPoint(Point *pnt)
      
      For a given point pnt, find a closed point of geometry
      entity
      
      04.2011. WW
   */
   Point* ConditionData::GetClosedPoint(const Point *pnt, const double tol)
   {
      Point *bc_pnt;
 
      real dist; 
      real min_dist = DBL_MAX;
      bc_pnt = NULL;

      if(ply)
      {
         min_dist = ply->MinDisttanceTo_a_Point(pnt);
         bc_pnt = ply->points[0];
      }

      else if(point)
      {
         dist = point->GetDistanceTo(pnt);
         if(dist<min_dist)
         {
            min_dist = dist;
            bc_pnt = point; 
         }

      }

      if(min_dist>tol)
        return NULL;
 
      return bc_pnt;

   }

   /// Output boundary condition 
   void ConditionData::Write(ostream &os)
   {
      os<<"\t geometry: ";
      if(ply)
         os<<"polyline "<<ply->Name()<<endl;
      else if(point)
         os<<"point "<<point->Index()<<endl;
      os<<"\t value: "<<value<<endl<<endl; 
   }

   /// As the funtion name
   void ConditionData::SetGeoEntityType(string type_name)
   {
       BC_Type pnt_bc_type;
       if(type_name.find("neumann")!=string::npos)
          pnt_bc_type = Neumann;  //Defaults
       else if(type_name.find("dirichlet")!=string::npos)
          pnt_bc_type = Dirichlet;
       else if(type_name.find("source")!=string::npos)
          pnt_bc_type = Source_term;
      
       if(point)
          point->bc_type = pnt_bc_type;
       if(ply)
       {
          for(int i=0; i<(int)ply->points.size(); i++)
            ply->points[i]->bc_type = pnt_bc_type; 
       }   
         
   }

}