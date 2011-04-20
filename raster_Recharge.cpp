/*!
\file Definition of class Raster_Recharge
*/
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include "raster_Recharge.h"

/// Constrcutor
Raster_Recharge::Raster_Recharge(string fpath, istream &ins)
{
  string aline;
  double step;
  std::stringstream ss;
  string key, uname;  

  file_path =  fpath;
 
  getline(ins, aline);     
  ss.str(aline);
  ss>> key >> uname;
  ss.clear();
   
  getline(ins, aline);     
  ss.str(aline);
  ss>> key >> ratio;
  ss.clear();

     
  step = 0.;

  while(!ins.eof())
  {       
     getline(ins, aline);     
     ss.str(aline);
     ss>> key;
     ss.clear();       
     step += 1.0; 
     precip_times.push_back(step);
     precip_files.push_back(key);  
  }

  /// GIS_shape_head[0]:  ncols
  /// GIS_shape_head[1]:  nrows
  /// GIS_shape_head[2]:  xllcorner
  /// GIS_shape_head[3]:  yllcorner
  /// GIS_shape_head[4]:  cellsize
  /// GIS_shape_head[5]:  NONDATA_value
  GIS_shape_head = new double[6]; 
  recharge_cell_value = NULL;    
}

//--------------------------------------------------------------------------------------------------
void Raster_Recharge::Read_Raster(const double current_time)
{
  int i, k, size; 
  double tim_0, tim_1; 
  long l;

  std::string aline;
  std::stringstream ss;

  string fileA;

  size = (int)precip_times.size();
  k = -1;
  for(i=0; i<size-1; i++)
  {
     tim_0 = precip_times[i];
     tim_1 = precip_times[i+1];
     if((current_time>tim_0||fabs(current_time-tim_0)<DBL_MIN)&&
       (current_time<tim_1))
     {
        k = i;
        break;
     }             
  }
  if(k==-1)
    fileA = precip_files[0];
  else
    fileA = precip_files[k];

    
  fileA = file_path + fileA;
  std::ifstream ins(fileA.c_str());
  if(!ins.good())
  {
      std::cout << "Can not find file " << std::endl;
      return ;
  }

  getline(ins, aline);
  ss.str(aline);
  ss>> aline >> ncols;
  ss.clear();

  getline(ins, aline);
  ss.str(aline);
  ss>> aline >> nrows;
  ss.clear();

  getline(ins, aline);
  ss.str(aline);
  ss>> aline >> x0;
  ss.clear();

  getline(ins, aline);
  ss.str(aline);
  ss>> aline >> y0;
  ss.clear();

  getline(ins, aline);
  ss.str(aline);
  ss>> aline >> csize;
  ss.clear();

  getline(ins, aline);
  ss.str(aline);
  ss>> aline >> ndata_v;
  ss.clear();

  if(!recharge_cell_value)
    recharge_cell_value = new double[nrows * ncols];

  for(l=0; l<nrows * ncols; l++)
    ins>>recharge_cell_value[l];
  ins.close();
}


double Raster_Recharge::Assign_Grid_Point(const double x, const double y)
{
   long nx, ny;
   double val;

   nx = (long)((x-x0)/csize);
   ny = (long)((y-y0)/csize);
   ny = nrows-ny;
   if(ny<0) ny = 0;
   if(ny>nrows) ny = nrows;

   if(nx*csize+x0>=x)  nx -= 1;
   if(ny*csize+y0>=y)  ny -= 1;
   if(nx>=ncols-1) nx = ncols-2;
   if(ny>=nrows-1) ny = nrows-2;
   if(nx<0) nx = 0;
   if(ny<0) ny = 0;

   val = recharge_cell_value[ncols*ny+nx];
   if(fabs(val-ndata_v)<DBL_MIN)
      val = 0.;

   return val*ratio;
}
