// -*- C++ -*-
//----------------------------------------------------------------------
//
//                    William S. Klug, Luigi Perotti
//                University of California Los Angeles
//                   (C) 2004-2007 All Rights Reserved
//
//----------------------------------------------------------------------
//
// $Log$
//
//----------------------------------------------------------------------

/*! 
  \file ShapeBrick9.h

  \brief Class for linear in the plane /quadratic in the thickness 9-Node Triangular Brick Finite Element shape functions.

*/

#if !defined(__ShapeBrick9_h__)
#define __ShapeBrick9_h__

#include "Shape.h"

namespace voom
{
  //! Linear in the plane/ quadratic in the thickness 9-node triangular brick isoparametric shape functions.
  /*! These 9-node brick shape functions are to be used for 3-D
    problems.  The nodal numbering arrangement
    and the geometry of the standard domain are determined by the
    method <tt>nodalCoordinates()</tt>.
  */
  class ShapeBrick9 : public Shape<3> {
    
  public:

    //! constructor
    ShapeBrick9(const CoordinateArray & s) {
      const unsigned int n=9;
      _functions.resize(n);
      _derivatives.resize(n);
      _positions.resize(n);

      _positions[0] = 0.0, 0.0, 0.0;
      _positions[1] = 1.0, 0.0, 0.0;
      _positions[2] = 0.0, 1.0, 0.0;
      _positions[3] = 0.0, 0.0, 0.5;
      _positions[4] = 1.0, 0.0, 0.5;
      _positions[5] = 0.0, 1.0, 0.5;
      _positions[6] = 0.0, 0.0, 1.0;
      _positions[7] = 1.0, 0.0, 1.0;
      _positions[8] = 0.0, 1.0, 1.0;

      compute(s);
    }

    //! compute the shape functions and derivatives
    void compute(const CoordinateArray & s);

    //! Return parametric coordinates of nodes.
    /*! Nodal arrangement   	\verbatim
      Node 1 = 0.0, 0.0, 0.0;
      Node 2 = 1.0, 0.0, 0.0;
      Node 3 = 0.0, 1.0, 0.0;
      Node 4 = 0.0, 0.0, 0.5;
      Node 5 = 1.0, 0.0, 0.5;
      Node 6 = 0.0, 1.0, 0.5;
      Node 7 = 0.0, 0.0, 1.0;
      Node 8 = 1.0, 0.0, 1.0;
      Node 9 = 0.0, 1.0, 1.0;  \endverbatim
    
    */
    PositionContainer nodalCoordinates() {
       return _positions;
       }
  };
}


#endif //#define __ShapeBrick9_h__
