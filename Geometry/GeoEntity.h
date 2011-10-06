/*!
  \file Geo_Entity.h
   
   Declaration of class Geo_Entity
   
   13.04.2011. WW
*/
#ifndef Geo_Entity_INC
#define Geo_Entity_INC

#include "misc.h"
#include<iostream>
#include<fstream>
#include<cmath>

namespace _FDM {class FiniteDifference; class ConditionDataBC;}

namespace Geometry_Group
{
   
   using _FDM::FiniteDifference;
   using _FDM::ConditionDataBC; 
   /*!
     \class  Point
         
   */
   enum Point_Type {none, intern, border, nm_11, nm_12, nm_13, nm_14, nm_21, nm_22, nm_23, nm_24};
   enum BC_Type{Neumann, Dirichlet, Source_term};
   enum NeighborCell_Type {NE, NW, SE, SW};
   enum NeighborPoint_Type {C, E, N, W, S};
   enum Geo_Type {point, ply};
  
   class Geo_Entity
   {
      public:
        Geo_Entity(Geo_Type gtyp) { gtype = gtyp; }
        ~Geo_Entity() {} 
 
      private: 
        Geo_Type gtype;
   };
   
}


#endif

