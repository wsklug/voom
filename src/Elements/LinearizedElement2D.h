// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                 (C) 2004-2005 All Rights Reserved
//
//----------------------------------------------------------------------
//
// $Log$
// Revision 1.3  2005/05/23 17:27:53  klug
// Added index map; renamed strainEnergy as energy.
//
// Revision 1.2  2005/04/11 16:26:59  klug
// Renamed LinearizedElement.cc to LinearizedElement.icc
//
// Revision 1.1  2005/04/11 15:57:21  klug
// Initial checkin.
//
//----------------------------------------------------------------------

/*! 
  \file LinearizedElement2D.h

  \brief Class for a finite element.

*/

#if !defined(__LinearizedElement2D_h__)
#define __LinearizedElement2D_h__

#include <blitz/array.h>
#include <vector>
 
#include "Element.h"

#if defined(_B)
#undef _B
#endif

namespace voom
{

  //! Templated class for a 2-D small-strain finite element.
  /*!
    Material objects used with this class should be derived from LinearizedMaterial.
    Node objects should be either DeformationNode<2> or derived from it.
  */
  template< class Node_t,
	    class Quadrature_t,
	    class Material_t,
	    class Shape_t >
  class LinearizedElement2D 
    : public Element< Node_t, Quadrature_t, Material_t, Shape_t >
  {
    
  public:

    typedef typename Node_t::PositionVector PositionVector;

    typedef Element<Node_t,Quadrature_t,Material_t,Shape_t> Base;

    typedef typename Base::NodeContainer 		NodeContainer;
    typedef typename Base::NodeIterator 		NodeIterator;
    typedef typename Base::ConstNodeIterator 		ConstNodeIterator;
    typedef typename Base::QuadPointIterator 		QuadPointIterator;
    typedef typename Base::ConstQuadPointIterator 	ConstQuadPointIterator;
    typedef typename Base::ElementVector 		ElementVector;
    typedef typename Base::ElementMatrix 		ElementMatrix;
    typedef typename Base::DofIndexMap 			DofIndexMap;

    //! 3-D array for strain-displacement matrices.
    /*!  These are evaluated at quadrature points, and multiply
      the vector of nodal displacements to yield strain "vector" in
      Voigt notation: 

      B(p,alpha,ai), p=0,...,\#quad-1, alpha=0,1,2, ai=2a+i=0,...,2*\#nodes-1.
      
    */
    typedef typename std::vector< blitz::Array<double,2> > BmatrixContainer; 

    //! default constructor
    LinearizedElement2D();

    //! construct from nodes
    LinearizedElement2D(const Quadrature_t & quad,
			const Material_t & mat,
			const NodeContainer & nodes) 
    { _initialize(quad,mat,nodes); }

    //! destructor
    ~LinearizedElement2D() {;}
    
    //! Do mechanics on element; compute energy, forces, and/or stiffness.
    void compute(bool f0, bool f1, bool f2, bool updateDOF=true);

  private:
  
    BmatrixContainer _B;

    void _initialize(const Quadrature_t & quad,
		     const Material_t & mat, 
		     const NodeContainer & nodes);
    
    //! Copy nodal field values to local array.
    void _copyNodalField();

  };  // end of class declaration
	
} // namespace voom


#include "LinearizedElement2D.icc"

#endif // __LinearizedElement2D_h__
