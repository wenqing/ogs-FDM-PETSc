#ifndef bc_INC
#define bc_INC
#include<iostream>
#include "misc.h"

namespace _FDM
{
   class Mat_Property;
   class Point;
   class Polyline;
    
   /*!
      \class BoundayCondition
       
      Manage Boundary conditions  
   */
   class BoundayCondition
   {
       public:
         BoundayCondition(ifstream &ins);
         ~BoundayCondition() {}

         void Write(ostream &os = cout);

         void SetGeoEntityType(string type_name);

         Point* GetClosedPoint(const Point *pnt, const double tol);
       private:

         //BC_Type bc_type;
         Point *point;
         Polyline *ply; 
         real value;   
         
         friend class FiniteDifference;   
   };

}

#endif