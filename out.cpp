/*!
  \brief
   Definitions class Outout. 

   05.2011 WW
*/
#include <sstream>
#include "out.h"
#include "misc.h"
using namespace std;
namespace _FDM
{
   Output::Output(ifstream &ins)
   {
     string aline;
     std::stringstream ss;
     float ti = 0.; 
     steps = -1;
     bool done = false; 

     for(;;)
     {
        getline(ins, aline); 
        if(aline.find("...")!=string::npos)
           break;
        aline = string_To_lower(aline);
        if(aline.find("geometry")!=string::npos) 
        {
           ss.str(aline);
           /// skip key
           ss>>aline;
           ss>>aline;
           ss.clear();
           aline = string_To_lower(aline);
           if(aline.find("domain")!=string::npos)
           {
              fname = file_name +"_domain";  
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
