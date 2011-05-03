/*! \class Raster_Recharge_H
    \brief A class for recharge data

      This class handle a series raster files, which save
    the recharge data. Each raster file contains the rechange   
    03.2011. WW
*/
#ifndef Raster_Recharge_H
#define Raster_Recharge_H 

#include<vector>
#include<iostream>
#include<string>
using namespace std;
class Raster_Recharge
{
   public:
    Raster_Recharge();
    ~Raster_Recharge()
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


    vector<double> precip_times;
    vector<string> precip_files;

};
#endif