#ifndef ConditionDataBC_INC
#define ConditionDataBC_INC
#include<iostream>
#include "misc.h"

namespace _FDM
{
   class Mat_Property;
   class Point;
   class Polyline;
    
   /*!
      \class ConditionDataBC
       
      Manage Initial Boundary conditions  
   */
   class ConditionDataBC
   {
       public:
         ConditionDataBC(ifstream &ins);
         ~ConditionDataBC() {}

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