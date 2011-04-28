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
#include <limits>

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

      /// Generate a linear solver
      sp = new SparseTable(this);
      eqs = new Linear_EQS(*sp, 1);      

      //sp->Write();

      WriteGrid_VTK();
     
   }
   

//----------------------------------------------------------
   /// Destructor
   FiniteDifference::~FiniteDifference()   
   {
      delete mat;
      delete sp;
      delete eqs;
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
        int k, l;
        long i, j;
        double x0, x1, y0, y1;

        /// Sort neighbor point index
        long idx_buff[5], buff1, buff0;
        NeighborPoint_Type nbp_type; 


        long size = (ncols+1)*(nrows+1);

        cell_status.resize(ncols*nrows);
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
              
              cell_status[i*ncols+j] = false;
              x0 = xll0 + cell_size*j;
              x1 = xll0 + cell_size*(j+1);

              if(  boundary->PointInDomain(x0,y0)
                 &&boundary->PointInDomain(x1,y0)
                 &&boundary->PointInDomain(x1,y1)
                 &&boundary->PointInDomain(x0,y1))
              {  
                 cell_status[i*ncols+j] = true;
                 if(pnt_eqs_index[i*(ncols+1)+j] == -1)
                 {
                    Point *new_grid_pnt = new Point((long)grid_point_in_use.size(), x0, y0);
                    new_grid_pnt->point_type = intern;
                    pnt_eqs_index[i*(ncols+1)+j] = new_grid_pnt->Index();
                    new_grid_pnt->grid_i = i;
                    new_grid_pnt->grid_j = j; 
                    grid_point_in_use.push_back(new_grid_pnt);
                 }                 
                 if(pnt_eqs_index[(i+1)*(ncols+1)+j] == -1)
                 {
                    Point *new_grid_pnt = new Point((long)grid_point_in_use.size(), x0, y1);
                    new_grid_pnt->point_type = intern;
                    pnt_eqs_index[(i+1)*(ncols+1)+j] = new_grid_pnt->Index();
                    new_grid_pnt->grid_i = i+1;
                    new_grid_pnt->grid_j = j; 
                    grid_point_in_use.push_back(new_grid_pnt);
                 }                 
                 if(pnt_eqs_index[(i+1)*(ncols+1)+j+1] == -1)
                 {
                    Point *new_grid_pnt = new Point((long)grid_point_in_use.size(), x1, y1);
                    new_grid_pnt->point_type = intern;
                    pnt_eqs_index[(i+1)*(ncols+1)+j+1] = new_grid_pnt->Index();
                    new_grid_pnt->grid_i = i+1;
                    new_grid_pnt->grid_j = j+1; 
                    grid_point_in_use.push_back(new_grid_pnt);
                 }                 
                 if(pnt_eqs_index[i*(ncols+1)+j+1] == -1)
                 {
                    Point *new_grid_pnt = new Point((long)grid_point_in_use.size(), x1, y0);
                    new_grid_pnt->point_type = intern;
                    pnt_eqs_index[i*(ncols+1)+j+1] = new_grid_pnt->Index();
                    new_grid_pnt->grid_i = i+1;
                    new_grid_pnt->grid_j = j+1; 
                    grid_point_in_use.push_back(new_grid_pnt);
                 }                 
                  
                 /// Add this cell as a neighbor 
                 grid_point_in_use[pnt_eqs_index[i*(ncols+1)+j]]->neighbor_cell_type.push_back(NE);
                 grid_point_in_use[pnt_eqs_index[i*(ncols+1)+j+1]]->neighbor_cell_type.push_back(NW);
                 grid_point_in_use[pnt_eqs_index[(i+1)*(ncols+1)+j+1]]->neighbor_cell_type.push_back(SW);
                 grid_point_in_use[pnt_eqs_index[(i+1)*(ncols+1)+j]]->neighbor_cell_type.push_back(SE);
              }
           }

        }   
     
         
        /// Configure the topology of the grid
        long jw, je, jn, js;
        for(i=0; i<nrows+1; i++)
        {
           for(j=0; j<ncols+1; j++)
           {
              /// If the point is not in use              
              if(pnt_eqs_index[i*(ncols+1)+j] == -1)
                 continue;
             
              Point *pnt = grid_point_in_use[pnt_eqs_index[i*(ncols+1)+j]];

              /// This point is the center
              pnt->neighbor_points.push_back(pnt->Index());
              pnt->np_position.push_back(C);

              je = i*(ncols+1)+j+1;
              jw = i*(ncols+1)+j-1;
              jn = (i+1)*(ncols+1)+j;
              js = (i-1)*(ncols+1)+j;

              
            
              if(pnt->neighbor_cell_type.size() ==1 )   
              {
                 switch(pnt->neighbor_cell_type[0])
                 {
                   case NE:
                     //   |
                     //   |
                     //   ----
                     pnt->point_type = nm_24;  
                     pnt->neighbor_points.push_back(grid_point_in_use[pnt_eqs_index[jn]]->Index());
                     pnt->np_position.push_back(N);
                     pnt->neighbor_points.push_back(grid_point_in_use[pnt_eqs_index[je]]->Index());
                     pnt->np_position.push_back(E);
                     break;
                   case NW:
                     //      |
                     //      |
                     //   ----
                     pnt->point_type = nm_23;  
                     pnt->neighbor_points.push_back(grid_point_in_use[pnt_eqs_index[jn]]->Index());
                     pnt->np_position.push_back(N);
                     pnt->neighbor_points.push_back(grid_point_in_use[pnt_eqs_index[jw]]->Index());
                     pnt->np_position.push_back(W);
                     break;
                   case SW:
                     //  ---
                     //     |
                     //     |
                     pnt->point_type = nm_22;  
                     pnt->neighbor_points.push_back(grid_point_in_use[pnt_eqs_index[js]]->Index());
                     pnt->np_position.push_back(S);
                     pnt->neighbor_points.push_back(grid_point_in_use[pnt_eqs_index[jw]]->Index());
                     pnt->np_position.push_back(W);
                     break;
                   case SE:
                     //   ---
                     //  |
                     //  |
                     pnt->neighbor_points.push_back(grid_point_in_use[pnt_eqs_index[je]]->Index());
                     pnt->np_position.push_back(E);
                     pnt->neighbor_points.push_back(grid_point_in_use[pnt_eqs_index[js]]->Index());
                     pnt->np_position.push_back(S);
                     pnt->point_type = nm_21;  
                     break;
                 
                 } 
              }
              else if(pnt->neighbor_cell_type.size() ==2 )
              {
                 /// 
                 ///-----.---.---.-----
                 ///-------------------
                 ///-------------------                  
                 if(   (pnt->neighbor_cell_type[0] == NE&&pnt->neighbor_cell_type[1] == NW)
                     ||(pnt->neighbor_cell_type[1] == NE&&pnt->neighbor_cell_type[0] == NW))
                 {
                     pnt->point_type = nm_14;  
                     
                     pnt->neighbor_points.push_back(grid_point_in_use[pnt_eqs_index[je]]->Index());
                     pnt->np_position.push_back(E);
                     pnt->neighbor_points.push_back(grid_point_in_use[pnt_eqs_index[jw]]->Index());
                     pnt->np_position.push_back(W);
                 }
                 ///       ||||
                 ///       .|||
                 ///       ||||
                 else if( (pnt->neighbor_cell_type[0] == SW&&pnt->neighbor_cell_type[1] == NW)
                        ||(pnt->neighbor_cell_type[1] == SW&&pnt->neighbor_cell_type[0] == NW))
                 {
                     pnt->point_type = nm_12;  

                     pnt->neighbor_points.push_back(grid_point_in_use[pnt_eqs_index[jn]]->Index());
                     pnt->np_position.push_back(N);
                     pnt->neighbor_points.push_back(grid_point_in_use[pnt_eqs_index[js]]->Index());
                     pnt->np_position.push_back(S);
                     
                 }
                 ///||||
                 ///|||.
                 ///||||
                 else if( (pnt->neighbor_cell_type[0] == SE&&pnt->neighbor_cell_type[1] == NE)
                        ||(pnt->neighbor_cell_type[1] == SE&&pnt->neighbor_cell_type[0] == NE))
                 {
                     pnt->point_type = nm_11;  
                     
                     pnt->neighbor_points.push_back(grid_point_in_use[pnt_eqs_index[jn]]->Index());
                     pnt->np_position.push_back(N);
                     pnt->neighbor_points.push_back(grid_point_in_use[pnt_eqs_index[js]]->Index());
                     pnt->np_position.push_back(S);
                 }
                 ///-------------------
                 ///-------------------                  
                 ///-----.---.---.-----
                 /// 
                 else if( (pnt->neighbor_cell_type[0] == SW&&pnt->neighbor_cell_type[1] == SE)
                        ||(pnt->neighbor_cell_type[1] == SW&&pnt->neighbor_cell_type[0] == SE))
                 {
                     pnt->point_type = nm_13;  
                     
                     pnt->neighbor_points.push_back(grid_point_in_use[pnt_eqs_index[je]]->Index());
                     pnt->np_position.push_back(E);
                     pnt->neighbor_points.push_back(grid_point_in_use[pnt_eqs_index[jw]]->Index());
                     pnt->np_position.push_back(W);
                 }

              }

              // This is a point inside the domain or the point is at a corner but has
              // four neighbours.
              else if(pnt->neighbor_cell_type.size()>2) 
              {
                  pnt->neighbor_points.push_back(grid_point_in_use[pnt_eqs_index[jn]]->Index());
                  pnt->np_position.push_back(N);
                  pnt->neighbor_points.push_back(grid_point_in_use[pnt_eqs_index[je]]->Index());
                  pnt->np_position.push_back(E);
                  pnt->neighbor_points.push_back(grid_point_in_use[pnt_eqs_index[jw]]->Index());
                  pnt->np_position.push_back(W);
                  pnt->neighbor_points.push_back(grid_point_in_use[pnt_eqs_index[js]]->Index());
                  pnt->np_position.push_back(S);    

                  if(pnt->neighbor_cell_type.size()==3)
                    pnt->point_type = border;
                  else if(pnt->neighbor_cell_type.size()==4)
                    pnt->point_type = intern;
          
              }


              for(k=0; k<(int)pnt->neighbor_points.size(); k++) 
                idx_buff[k] = pnt->neighbor_points[k];
              for(k=0; k<(int)pnt->neighbor_points.size(); k++) 
              {
                 buff0 = idx_buff[k];  
                 buff1 = pnt->neighbor_points[k];
                 nbp_type =  pnt->np_position[k]; 
                 l = k;
                 while(l>0&&idx_buff[l-1]>buff0)  
                 {
                    idx_buff[l] = idx_buff[l-1];
                    pnt->neighbor_points[l] = pnt->neighbor_points[l-1];
                    pnt->np_position[l] = pnt->np_position[l-1]; 
                    l--;
                 }
                 idx_buff[l] = buff0;  
                 pnt->neighbor_points[l] = buff1;
                 pnt->np_position[l] = nbp_type; 
              }
              
              /// Assign boundary condition to this point, if it is close to 
              /// the geometry entity for the boundary conditions.
              if(!CheckDirichletBC(pnt))
              {
                 if(pnt->point_type != intern)
                   CheckNuemannBC(pnt); 
              }
             
           }
        }

    } 
    //---------------------------------
    /*!
       \fn inline bool CheckDirichletBC(Point *pnt)
       
       Determine whether a grid point on the boundary is
       a point assigned with the Dirichlet boundary by geometry
       entity, which mush be with threhold of the grid size. 
 
       04.2011. WW 
       
    */ 

    inline bool FiniteDifference::CheckDirichletBC(Point *pnt)
    {
        int i;
        Point *bc_pnt;

        bc_pnt = NULL;
        for(i=0; i<(int)BC_Dirichlet.size(); i++)
        {
           bc_pnt = BC_Dirichlet[i]->GetClosedPoint(pnt, cell_size);
           if(bc_pnt)
             bc_pnt->value = BC_Dirichlet[i]->value;
        }
         
        if(!bc_pnt)
          return false;

              
        pnt->bc_type = Dirichlet;
        pnt->value = bc_pnt->value;


        return true;
       
         
    }

    //---------------------------------
    /*!
       \fn inline void CheckNeumannBC(Point *pnt)
       
       Determine whether a grid point on the boundary is
       a point assigned with the Neumann boundary by geometry
       entity, which mush be with threhold of the grid size. 
 
       04.2011. WW 
       
    */ 

    inline void FiniteDifference::CheckNuemannBC(Point *pnt)
    {
        int i;
        Point *bc_pnt;
        bc_pnt = NULL;

        for(i=0; i<(int)BC_Neumann.size(); i++)
        {
           bc_pnt = BC_Neumann[i]->GetClosedPoint(pnt, cell_size);
           if(bc_pnt)
             bc_pnt->value = BC_Neumann[i]->value;
        }
         
        if(!bc_pnt)
          return;

        pnt->bc_type = Neumann;
        pnt->value = bc_pnt->value;
         
    }
    //---------------------------------
    /*!
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
              if(cell_status[i*ncols+j])
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


