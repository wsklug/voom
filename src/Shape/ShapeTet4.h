// -*- C++ -*-
//----------------------------------------------------------------------
//
//                 William S. Klug, Melissa M. Gibbons
//                University of California Los Angeles
//                 (C) 2004-2006 All Rights Reserved
//
//----------------------------------------------------------------------
//
// $Log$
//
//----------------------------------------------------------------------

/*! 
  \file ShapeTet4.h

  \brief Class for Linear 4-Node Tetrahedral Finite Element shape functions.

*/

#if !defined(__ShapeTet4_h__)
#define __ShapeTet4_h__

#include "Shape.h"

namespace voom
{
  //! Linear 4-node tetrahedral isoparametric shape functions.
  /*! These 4-node tetrahedral shape functions are to be used for 3-D
    problems.  Three of the four barycentric coordinates are
    used as curvilinear coordinates.  The nodal numbering arrangement
    and the geometry of the standard domain are determined by the
    method <tt>nodalCoordinates()</tt>.
  */
  class ShapeTet4 : public Shape<3,3> {
    
  public:

    //! Constructor.  Calls compute function.
    ShapeTet4(const CoordinateArray & s, const PositionContainer & x) 
      : Shape<3,3>(x) {
      std::size_t nNodes = x.size();
      if(nNodes!=4) {
	std::cout << "ShapeTet4::ShapeTet4() needs 4 nodal position vectors; "
		  << nNodes << " provided."<<std::endl;
	exit(0);
      }
      compute(s);
    }

    //! compute the shape functions and derivatives
    void compute(const CoordinateArray & s);

    //! Return parametric coordinates of nodes.
    /*! Nodal arrangement   	\verbatim
      Node 1 = 1.0, 0.0, 0.0;
      Node 2 = 0.0, 1.0, 0.0;
      Node 3 = 0.0, 0.0, 1.0;
      Node 4 = 0.0, 0.0, 0.0;    \endverbatim
    
    */
    PositionContainer nodalCoordinates() {
      PositionContainer sa(4);
      sa[0] = 1.0, 0.0, 0.0;
      sa[1] = 0.0, 1.0, 0.0;
      sa[2] = 0.0, 0.0, 1.0;
      sa[3] = 0.0, 0.0, 0.0;
      return sa;
    }
  };
}


#endif //#define __ShapeTet4_h__
