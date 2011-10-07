/*!
  \brief
   Definitions class Outout. 

   05.2011 WW
*/

#include "Output.h"

#include <sstream>

#include "misc.h"
#include "Point.h"
#include "Polyline.h"
#include "Geometry.h"

namespace _FDM
{
   using namespace std;
   using namespace Geometry_Group;
   using namespace AuxFunctions;
	

   Output::Output(string f_path, string f_name, ifstream &ins,  Geometry_Group::Geometry *geometry)
   {
     string aline;
     std::stringstream ss;

     file_path = f_path;
     file_name = f_name; 	

     int ID;
     float ti = 0.; 
     steps = -1;
     bool done = false; 
     bool geo = false; 

     geo_entity = NULL;

     for(;;)
     {
        getline(ins, aline);  

        if(CheckComment(aline))
           continue;

        if(aline.find("...")!=string::npos)
           break;
		aline = AuxFunctions::string_To_lower(aline);
        if(aline.find("geometry")!=string::npos) 
        {
           ss.str(aline);
           /// skip key
           ss>>aline;
           ss>>aline;
           ss.clear();
		   aline = AuxFunctions::string_To_lower(aline);
           geo = false; 

           if(aline.find("domain")!=string::npos)
           {
              geo = true; 
              fname = file_name +"_domain";  
           }
           else if(aline.find("polyline")!=string::npos)
           {
              geo = true; 

              ss.str(aline);
              /// skip key
              ss>>aline;
              ss>>aline;  
              ss.clear();
              Polyline *ply = geometry->getPolylineByName(aline);       

              fname = file_name +"_ply" + ply->Name();
              geo_entity = ply;  
           }
           else if(aline.find("point")!=string::npos)
           {
              geo = true; 
              ss.str(aline);
              /// skip key
              ss>>aline;
              /// point ID
              ss>>ID; 
              ss.clear();

              geo_entity = geometry->getPointByID(ID); 
              ss<<ID; 
              fname = file_name +"_point" + ss.str();
              ss.clear();
           }

           if(geo)
           {
              /// If Data of all steps can be dumped to a single file, use: 
              //      os = new ofstream(fname.c_str(), ios::trunc);
              os = new ofstream(); //(fname.c_str(), ios::trunc);
              os->clear();
              os->close();
              continue; 

           }  
        } 
        if(aline.find("times")!=string::npos) 
        {
           for(;;)
           {
              getline(ins, aline); 
              if(aline.find("...")!=string::npos)
              {
                 done = true;
                 break;
              }
              ss.str(aline);
              ss>>ti;
              at_times.push_back(ti);
              ss.clear();               
           } 
        }
        if(done)
          break;
        if(aline.find("steps")!=string::npos) 
        {
           ss.str(aline);
           /// skip key
           ss>>aline;
           ss>>steps;
        }

     }
  

     
   }

   Output::~Output()
   {
      os->close();
      delete os;
   }
   void Output::Write(std::ostream &os)
   {
      os<<"--- Output"<<endl;  
      os<<"\t geometry: "<<"domain"<<endl;
      if(steps == -1)  
      {  
         os<<"\t times: "<<endl; 
         for(int i=0; i<(int)at_times.size(); i++)
            os<<"\t"<<at_times[i]<<endl;
         os<<endl;
      } 
      else
         os<<"\t step: "<<steps<<endl; 
      os<<"\t ..."<<endl;
   } 
}

