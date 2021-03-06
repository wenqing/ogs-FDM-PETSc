/*! \class Raster_Recharge_H
    \brief A class for recharge data

      This class handle a series raster files, which save
    the recharge data. Each raster file contains the rechange
    03.2011. WW
*/
#ifndef RasterRecharge_H
#define RasterRecharge_H

#include<vector>
#include<iostream>
#include<string>

class RasterRecharge
{
   public:
      RasterRecharge(std::string f_path, std::string f_name);
      ~RasterRecharge()
      {
         delete [] GIS_shape_head;
         if(recharge_cell_value)
            delete [] recharge_cell_value;
      }
      /*!
         \fn Read(const double current_time);
         Open a raster file.
         \param  current_time :current time
      */
      void Read_Raster(const double current_time);
      /*!
         \fn Assign_Grid_Point();
         Calculate the point infiltration.
         \param  x            :x coordinate
         \param  y            :y coordinate
      */
      double Assign_Grid_Point(const double x, const double y);
   private:

      double ratio;
      double *GIS_shape_head;
      double *recharge_cell_value;

      /// Import shape file.
      long ncols, nrows;
      /// (x_0, y_0): coordinate of the left down corner
      double x0, y0, csize, ndata_v;


      std::vector<double> precip_times;
      std::vector<std::string> precip_files;

      std::string file_name;
      std::string file_path;

};
#endif

