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
      \class ConditionData
       
      Manage Initial Boundary conditions  
   */
   class ConditionData
   {
       public:
         ConditionData(ifstream &ins);
         ~ConditionData() {}

         void Write(ostream &os = cout);

         void SetGeoEntityType(string type_name);

         Point* GetClosedPoint(const Point *pnt, const double tol);
       private:


         Point *point;
         Polyline *ply; 
         real value;   
         
         friend class FiniteDifference;   
   };

}

#endif