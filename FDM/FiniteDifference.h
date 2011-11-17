#ifndef FiniteDifference_INC
#define FiniteDifference_INC

#include<iostream>

#ifdef USE_PETSC
#include "PETScLinearSolver.h"
#endif

#include "Point.h"

/*!
 \class  FiniteDifference
 
   A class of the finite difference caculation of 
   groundwater flow

   WW 04.2011
   
*/
namespace Geometry_Group{ class Geometry;}

namespace Math_Group {class Linear_EQS; class SparseTable;}

class RasterRecharge;
namespace _FDM
{
   class Mat_Property;
   class ConditionDataBC;
   class Numerics;
   class Output;

   using Math_Group::Linear_EQS;
   using Math_Group::SparseTable;

   using Geometry_Group::Point;
   using Geometry_Group::Polyline;


//--------------- class  FiniteDifference ---------------------
   class FiniteDifference 
   {
      public:
        FiniteDifference(std::string f_path, std::string f_name)
          : eqs(NULL),geo_grid(NULL),mat(NULL), num(NULL), ic(NULL),
            rrecharge(NULL),  boundary(NULL),
            file_name(f_name), file_path(f_path) 		  
		  {  }
        ~FiniteDifference();   
	
        void Initialize();
        
        void TimeSteping(); 
        void AssembleEQS(); 

        void Output_Results(const float c_tim, const int i_step);
        void Write(std::ostream &os = std::cout);
        
        void WriteGrid_VTK();

#ifdef USE_PETSC
        void set_MPI_rank_size(const int rank, const int size)
        {
           size_MPI = size;
           rank_MPI = rank;  
        } 
#endif

      private:

        /// Unknowns
        double *u0;
        double *u1;

#ifdef USE_PETSC
        int size_MPI;
        int rank_MPI;
        PETScLinearSolver *eqs;
#else
        /// EQS
        SparseTable *sp; 
        Linear_EQS *eqs;
#endif

        /// Data for geometry and grid
        long nrows;
        long ncols;
        /// Coodinates of the low left corner of the grid 
        float xll0;
        float yll0;        
        /// Cell size
        float cell_size; 

        Geometry_Group::Geometry *geo_grid; 
 
        /// Time step
        float dt;
        /// Start time 
        float T0;
        /// End time 
        float T1;
        /// Time factor
        real tim_fac;
        ///
        std::string time_unit; 

        /// Material
        Mat_Property *mat;

        /// Numerics
        Numerics *num; 

        /// Boundary Condition
        std::vector<ConditionDataBC*> BC_Neumann; 
        std::vector<ConditionDataBC*> BC_Dirichlet; 
        std::vector<ConditionDataBC*> Source_Sink; 
        std::vector<long> BC_Dirichlet_points; 
        ConditionDataBC *ic;
         
        /// Source term
        RasterRecharge *rrecharge;

        
        /// Output
        std::vector<Output*> outp;

        /// Boundary of Domain
        Polyline *boundary;

        /// Output control (so far only domain is considered)
        /// If other geometry
  
        /// Save equation indicies of grid points 
        std::vector<long> pnt_eqs_index;
        ///
        std::vector<Point*>  grid_point_in_use;
        
        /// Cell mark. 
        std::vector<bool> cell_status;
        ///
        long num_cell_in_use;
  
        void CaterogorizeGridPoint();


        std::string file_name; 
        std::string file_path; 

        /// set Dirichlet boundary condition
        bool CheckDirichletBC(Point *pnt);
        void CheckNuemannBC(Point *pnt);
        void CheckSourceSink(Point *pnt);
		void setBC_at_PointOnLine(long i, Point *pnt, Geometry_Group::NeighborPoint_Type nbt);
        void setBC_at_Point_atCCorner(long i, Point *pnt);
        
        void Output_Domain_VTK(std::ostream &os);

        friend class ConditionDataBC; 
        friend class Math_Group::SparseTable; 

   };
 
}
#endif

