#ifndef FiniteDifference_INC
#define FiniteDifference_INC

#include<iostream>

#include "misc.h"
#include "geo.h"

/*!
 \class  FiniteDifference
 
   A class of the finite difference caculation of 
   groundwater flow

   WW 04.2011
   
*/
namespace Math_Group {class Linear_EQS; class SparseTable;}

class RasterRecharge;
namespace _FDM
{
   class Mat_Property;
   class Point;
   class Polyline; 
   class ConditionData;
   class Numerics;
   class Output;

   using Math_Group::Linear_EQS;
   using Math_Group::SparseTable;

//--------------- class  FiniteDifference ---------------------
   class FiniteDifference 
   {
      public:
		  FiniteDifference() {}
        ~FiniteDifference();   

        void Initialize();
        void TimeSteping(); 
        void AssembleEQS(); 

        void Output_Results(const float c_tim, const int i_step);
        void Write(ostream &os = cout);
        void WriteGrid_VTK();

      private:

        /// Unknowns
        double *u0;
        double *u1;

        /// EQS
        SparseTable *sp; 
        Linear_EQS *eqs;


        /// Data for geometry and grid
        long nrows;
        long ncols;
        /// Coodinates of the low left corner of the grid 
        float xll0;
        float yll0;        
        /// Cell size
        float cell_size; 

 
        /// Time step
        float dt;
        /// Start time 
        float T0;
        /// End time 
        float T1;
        /// Time factor
        real tim_fac;
        ///
        string time_unit; 

        /// Material
        Mat_Property *mat;

        /// Numerics
        Numerics *num; 

        /// Boundary Condition
        vector<ConditionData*> BC_Neumann; 
        vector<ConditionData*> BC_Dirichlet; 
        vector<ConditionData*> Source_Sink; 
        vector<long> BC_Dirichlet_points; 
        ConditionData *ic;
         
        /// Source term
        RasterRecharge *rrecharge;

        
        /// Output
        vector<Output*> outp;

        /// Boundary of Domain
        Polyline *boundary;

        /// Output control (so far only domain is considered)
        /// If other geometry
  
        /// Save equation indicies of grid points 
        vector<long> pnt_eqs_index;
        ///
        vector<Point*>  grid_point_in_use;
        
        /// Cell mark. 
        vector<bool> cell_status;
        ///
        long num_cell_in_use;
  
        void CaterogorizeGridPoint();




        /// Set Dirichlet boundary condition
        bool CheckDirichletBC(Point *pnt);
        void CheckNuemannBC(Point *pnt);
        void CheckSourceSink(Point *pnt);
        void SetBC_at_PointOnLine(long i, Point *pnt, NeighborPoint_Type nbt);
        void SetBC_at_Point_atCCorner(long i, Point *pnt);
        
        void Output_Domain_VTK(ostream &os);

        friend class ConditionData; 
        friend class Math_Group::SparseTable; 

   };
 
}
#endif