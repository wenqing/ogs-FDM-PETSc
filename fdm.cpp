/*!
\file fdm.cpp
   
   The file contains the defitions of member functions of 
 class FiniteDifference

  WW 04.2011

*/

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "fdm.h"
#include "mat.h"
#include "geo.h"
#include "bc.h"
#include "numerics.h"
#include "equation_class.h"

using namespace std;
namespace _FDM
{
//--------------- class  FiniteDifference ---------------------
   /// Constructor 
   FiniteDifference:: FiniteDifference()
   {
      vector<string> key;
      vector<float> keyval;
      
      string dat_fname = file_name+".dat"; 
      ifstream ins(dat_fname.c_str());


      boundary = GetPolylineByName("boundary");
      if(!boundary)
      {
         cout<<"Boundary of the domain is not defined. "<<endl;
         cout<<"The boundary must be defined as polyline with name 'boundary' "<<endl;
         exit(1);
      }      

      
      if(!ins.good()) 
      {
         cout<<"Could not find file "<<dat_fname<<". Stop now!"<<endl;
         exit(1);
      } 

      string aline;
      mat = NULL;
      while(!ins.eof())
      {
         getline(ins, aline);

         //if(aline.find("...")!=string::npos)
         //   break;

         if(aline.find("#")!=string::npos)
           continue;

         aline = string_To_lower(aline);
         if(aline.find("material")!=string::npos)
         {
            mat = new Mat_Property(ins);
            continue;
         } 
         if(aline.find("solver")!=string::npos)
         {
            num = new Numerics(ins);
            continue;
         } 
        
         if(aline.find("grid")!=string::npos)
         {
            key.resize(5);
            keyval.resize(5);
     
            key[0] = "ncols:";
            key[1] = "nrows:";
            key[2] = "xllcorner:";
            key[3] = "yllcorner:";
            key[4] = "cellsize:";

            Read_Block(ins, key, keyval);
           
            ncols = (long)keyval[0]; 
            nrows = (long)keyval[1]; 
            xll0 = keyval[2]; 
            yll0 = keyval[3]; 
            cell_size = keyval[4];
            continue;                        
         }

         if(aline.find("time")!=string::npos)
         {

            if(aline.find("minut")!=string::npos)
            {
               tim_fac=60.0;
               time_unit = "minut";
            }
            else if(aline.find("hour")!=string::npos)
            {  
               tim_fac=3600.0;
               time_unit = "hour";
            }
            else if(aline.find("day")!=string::npos)
            {
               tim_fac=86400.0;
               time_unit = "day";
            } 
            else if(aline.find("month")!=string::npos)
            { 
               tim_fac=2592000.0;
               time_unit = "month";
            }
            else if(aline.find("year")!=string::npos) 
            {  
               tim_fac=31536000;
               time_unit = "year";                
            }

            key.resize(3);
            keyval.resize(3);
     
            key[0] = "start_time:";
            key[1] = "end_time:";
            key[2] = "step_size:";

            Read_Block(ins, key, keyval);
           
            T0 = keyval[0]; 
            T1 = keyval[1]; 
            dt = keyval[2]; 
            continue;
         }

         if(aline.find("neumann")!=string::npos)
         {
            BC_Neumann.push_back(new BoundayCondition(ins));
            BC_Neumann[BC_Neumann.size()-1]->SetGeoEntityType("neumann");
         } 
         if(aline.find("dirichlet")!=string::npos)
         {
            BC_Dirichlet.push_back(new BoundayCondition(ins));
            BC_Dirichlet[BC_Dirichlet.size()-1]->SetGeoEntityType("dirichlet");
         }
                  
      }
  
      CaterogorizeGridPoint();
      WriteGrid_VTK();
     
   }
   

//----------------------------------------------------------
   /// Destructor
   FiniteDifference::~FiniteDifference()   
   {
      delete mat;
   } 

//----------------------------------------------------------
   /*!
      \fn FiniteDifference::Write(ostream &os = cout)

       Output all paramters of the FDM
       
       04.2011 WW
       
   */
   void FiniteDifference::Write(ostream &os)
   {
       mat->Write(os);
       num->Write(os);   

       /// Data for geometry and grid
       os<<"--- Grid"<<endl;
       os<<"\t ncols:    \t"<< nrows  <<endl;
       os<<"\t nrows:    \t"<< ncols  <<endl;
       os<<"\t xllcorner:\t"<< xll0   <<endl;
       os<<"\t yllcorner:\t"<< yll0   <<endl;
       os<<"\t cellsize: \t"<< cell_size <<endl<<endl;

       /// Time step
       os<<"--- Time ("<<time_unit<<")"<<endl;
       os<<"\t start_time:\t"<<T0<<endl;
       os<<"\t end_time:\t"<<T1<<endl;
       os<<"\t step_size:\t"<<dt<<endl<<endl;

       /// Boundary condition
       int i;
       for(i=0; i<(int)BC_Neumann.size(); i++)
       {
          os<<"--- Neumann BC"<<endl;
          BC_Neumann[i]->Write(os);            
       }
       for(i=0; i<(int)BC_Dirichlet.size(); i++)
       {
          os<<"--- Dirichlet BC"<<endl;
          BC_Dirichlet[i]->Write(os);            
       }
 
       os<<"..."<<endl;
    }

    /*!
      \fn  AssembleEQS()
        
      Assemble equation system of FDM   

    */
    void FiniteDifference::AssembleEQS()
    {
         
    }
  
    /*!
      \fn  CaterogorizeGridPoint();
        
      To find the arribution of points, whethere they are interior,
      assgined with boundary conditions based on the theory by 
      B.J. Noye and R.J. Arnold, Accurate finite difference approxmations 
      for the Neumann conditions on a curved boundary, Appl. Math. Modelling, 
      1990, Vol. 14, pp2-13.   

    */
    void FiniteDifference::CaterogorizeGridPoint()
    {
        long i, j;
        double x0, x1, y0, y1;
        long size = (ncols+1)*(nrows+1);

        pnt_eqs_index.resize(size);
        for(i=0; i<size; i++)
          pnt_eqs_index[i] = -1;

        /// Loop over grid cells
        for(i=0; i<nrows; i++)
        {
           y0 = yll0 + cell_size*i;
           y1 = yll0 + cell_size*(i+1);
           for(j=0; j<ncols; j++)
           {
              x0 = xll0 + cell_size*j;
              x1 = xll0 + cell_size*(j+1);

              if(  boundary->PointInDomain(x0,y0)
                 &&boundary->PointInDomain(x1,y0)
                 &&boundary->PointInDomain(x1,y1)
                 &&boundary->PointInDomain(x0,y1))
              {  
                 if(pnt_eqs_index[i*(ncols+1)+j] == -1)
                 {
                    Point *new_grid_pnt = new Point((long)grid_point_in_use.size(), x0, y0);
                    pnt_eqs_index[i*(ncols+1)+j] = new_grid_pnt->Index();
                    new_grid_pnt->grid_i = i;
                    new_grid_pnt->grid_j = j; 
                    grid_point_in_use.push_back(new_grid_pnt);
                 }                 
                 if(pnt_eqs_index[(i+1)*(ncols+1)+j] == -1)
                 {
                    Point *new_grid_pnt = new Point((long)grid_point_in_use.size(), x0, y1);
                    pnt_eqs_index[(i+1)*(ncols+1)+j] = new_grid_pnt->Index();
                    new_grid_pnt->grid_i = i+1;
                    new_grid_pnt->grid_j = j; 
                    grid_point_in_use.push_back(new_grid_pnt);
                 }                 
                 if(pnt_eqs_index[(i+1)*(ncols+1)+j+1] == -1)
                 {
                    Point *new_grid_pnt = new Point((long)grid_point_in_use.size(), x1, y1);
                    pnt_eqs_index[(i+1)*(ncols+1)+j+1] = new_grid_pnt->Index();
                    new_grid_pnt->grid_i = i+1;
                    new_grid_pnt->grid_j = j+1; 
                    grid_point_in_use.push_back(new_grid_pnt);
                 }                 
                 if(pnt_eqs_index[i*(ncols+1)+j+1] == -1)
                 {
                    Point *new_grid_pnt = new Point((long)grid_point_in_use.size(), x1, y0);
                    pnt_eqs_index[i*(ncols+1)+j+1] = new_grid_pnt->Index();
                    new_grid_pnt->grid_i = i+1;
                    new_grid_pnt->grid_j = j+1; 
                    grid_point_in_use.push_back(new_grid_pnt);
                 }                 
                  
              }
           }

        }   
     
        //


    } 
    //---------------------------------
    /*
       \fn WriteGrid_VTK()
       
       Output Grid and boundary to VTK
 
       04.2011. WW 
       
    */ 
    void FiniteDifference::WriteGrid_VTK()
    {
       long i, j, n_points;
       double x, y;
       long size = (ncols+1)*(nrows+1);
       string dat_fname = file_name+"_grid.vtk"; 
       ofstream os(dat_fname.c_str(), ios::trunc);

       setw(12);
       os.precision(12);

       os<<"# vtk DataFile Version 4.0\nGrid of Oobs\nASCII\n"<<endl;
       os<<"DATASET UNSTRUCTURED_GRID"<<endl;
       os<<"POINTS "<<(nrows+1)*(ncols+1)<<" double"<<endl;

       /// Loop over grid points 
       for(i=0; i<nrows+1; i++)
       {
          y = yll0 + cell_size*i;
          for(j=0; j<ncols+1; j++)
          {
             x = xll0 + cell_size*j;             
             os<<x<<"  "<<y<<"  0."<<endl;
           
         }
       }
   
       os<<"\nCELLS "<<nrows*ncols<<" "<<nrows*ncols*5<<endl;
       /// Loop over grid cells
       for(i=0; i<nrows; i++)
       {
          for(j=0; j<ncols; j++)
          {
             os<<"4 ";
             os<<i*(ncols+1)+j<<" "<<i*(ncols+1)+j+1<<" "
               <<(i+1)*(ncols+1)+j+1<<" "<<(i+1)*(ncols+1)+j<<endl;            
          }

       }
       // CELL types
       os << "CELL_TYPES " << nrows*ncols << endl; 
       for(i=0; i<nrows*ncols; i++)
         os<<"9 "<<endl;

      
       os<<"POINT_DATA "<<size<<endl;
       os<<"SCALARS Points_in_Domain float 1\nLOOKUP_TABLE default"<<endl;
       for(i=0; i<size; i++)
         os<<(pnt_eqs_index[i]>-1 ? 1:0)<<endl;


       os<<"CELL_DATA "<<nrows*ncols<<endl;
       os<<"SCALARS Material_Group int 1\nLOOKUP_TABLE default"<<endl;
       int mat_id = 0;
       for(i=0; i<nrows; i++)
       {
           for(j=0; j<ncols; j++)
           {

              mat_id = 0;
              if(  pnt_eqs_index[i*(ncols+1)+j] > -1
                 &&pnt_eqs_index[(i+1)*(ncols+1)+j] > -1
                 &&pnt_eqs_index[(i+1)*(ncols+1)+j+1] > -1
                 &&pnt_eqs_index[i*(ncols+1)+j+1] > -1)
              mat_id = 1;
              os<<mat_id<<endl;

           }
       } 
       os.close();

       dat_fname = file_name+"_boundary.vtk"; 
       os.open(dat_fname.c_str(), ios::trunc);


       os<<"# vtk DataFile Version 4.0\nGrid of Oobs\nASCII\n"<<endl;
       os<<"DATASET UNSTRUCTURED_GRID"<<endl;
       os<<"POINTS "<<boundary->points.size()<<" double"<<endl;

       /// Loop over grid points 
       n_points = (long)boundary->points.size(); 
       for(i=0; i<n_points; i++)
         os<<boundary->points[i]->X()<<"  "<<boundary->points[i]->Y()<<"  0."<<endl;
       os<<endl;   
  
       os<<"\nCELLS "<<n_points<<" "<<n_points*3<<endl;
       /// Boundary 
       j = n_points-1; 
       for(i=0; i<n_points; i++)
       {
          os<<"2 "<<j<<" "<<i<<endl;
          j = i;
          
       }
      

       // CELL types
       os << "CELL_TYPES " << n_points << endl; 
       for(i=0; i<n_points; i++)
         os<<"3 "<<endl;
  

       os<<"POINT_DATA "<<n_points<<endl;
       os<<"SCALARS BC_type int 1\nLOOKUP_TABLE default"<<endl;
       for(i=0; i<n_points; i++)
         os<<boundary->points[i]->point_type<<endl;


       os.close();
    }

}


