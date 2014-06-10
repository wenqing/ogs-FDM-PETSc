#ifndef ConditionDataBC_INC
#define ConditionDataBC_INC
#include<iostream>
#include "misc.h"


namespace Geometry_Group
{
class Point;
class Polyline;
class Geometry;
}
namespace _FDM
{
class Mat_Property;

using Geometry_Group::Geometry;

/*!
   \class ConditionDataBC

   Manage Initial Boundary conditions
*/
class ConditionDataBC
{
   public:
      ConditionDataBC(std::ifstream &ins, Geometry *geo);
      ~ConditionDataBC() {}

      void Write(std::ostream &os = std::cout);

      void setGeoEntityType(std::string type_name);

      Geometry_Group::Point* getClosedPoint(const Geometry_Group::Point *pnt, const double tol);
   private:


      Geometry_Group::Point *point;
      Geometry_Group::Polyline *ply;
      real value;

      friend class FiniteDifference;
};

}

#endif
