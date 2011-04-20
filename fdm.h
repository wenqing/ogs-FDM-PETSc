#ifndef fdm_INC
#define fdm_INC

#include<iostream>

#include "misc.h"

/*!
 \class  FiniteDifference
 
   A class of the finite difference caculation of 
   groundwater flow

   WW 04.2011
   
*/
namespace Math_Group {class Linear_EQS;}

namespace _FDM
{
   class Mat_Property;
   class Point;
   class Polyline; 
   class BoundayCondition;
   class Numerics;

   using Math_Group::Linear_EQS;

//--------------- class  FiniteDifference ---------------------
   enum Grid_Point_Type {out, intern, nm_11, nm_12, nm_13, nm_14, nm_21, nm_22, nm_23, nm_24, dirichlet};
   class FiniteDifference 
   {
      public:
        FiniteDifference();
        ~FiniteDifference();   

        void AssembleEQS(); 

        void Write(ostream &os = cout);
        void WriteGrid_VTK();

      private:

        /// Unknowns
        double *u0;
        double *u1;

        /// EQS
        Linear_EQS *eqs;


        /// Data for geometry and grid
        long nrows;
        long ncols;
        /// Coodinates of the low left corner of the grid 
        real xll0;
        real yll0;        
        /// Cell size
        real cell_size; 

        ///
        vector<Point*>  grid_point_in_use;
        
 
        /// Time step
        real dt;
        /// Start time 
        real T0;
        /// End time 
        real T1;
        /// Time factor
        real tim_fac;
        ///
        string time_unit; 

        /// Material
        Mat_Property *mat;

        /// Numerics
        Numerics *num; 

        /// Boundary Condition
        vector<BoundayCondition*> BC_Neumann; 
        vector<BoundayCondition*> BC_Dirichlet; 

        /// Boundary of Domain
        Polyline *boundary;
  
        /// Save equation indicies of grid points 
        vector<long> pnt_eqs_index;
        void CaterogorizeGridPoint();

        friend class BoundayCondition; 

   };
 
}
#endif