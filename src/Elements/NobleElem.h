// -*- C++ -*-
//----------------------------------------------------------------------
//
//                   HoHai Van and William S. Klug
//                University of California Los Angeles
//                   (C) 2006 All Rights Reserved
//
//----------------------------------------------------------------------
//
//
//----------------------------------------------------------------------

#if !defined(__NobleElem_h__)
#define __NobleElem_h__

#include <vector>
#include "QuadQuadrature.h"
#include "ShapeQ4.h"

namespace voom
{
  /*! This class is the element class which calculates the quadrature for 
      each element in the mesh.  The input will be the node matrix consisting
      of the x,y coordinates, voltage, and dv/dt at each of the nodes in the
      local element and the quadrature matrix containing the quadrature point
      values of voltage, m, n, h gating variables.
  */



  class NobleElem 

  {
    
  public:
      

    //! NodeMatrix_t is a 4 x 4 matrix containing the x coord, y coord, v, dv/dt of each node 	
    typedef tvmet::Matrix<double, 4, 5> NodeMatrix_t;

    //! QuadMatrix_t is a 4x4 matrix containing the voltage, m, n, h at the quad points
    typedef tvmet::Matrix<double, 4, 4> QuadMatrix_t;

    typedef QuadQuadrature Quadrature_t;
    typedef ShapeQ4 Shape_t;
        
    //! default constructor
    NobleElem(const Quadrature_t & quad, NodeMatrix_t & nodes, QuadMatrix_t & quads, const double i_stim);

    //! virtual destructor
    virtual ~NobleElem() {;}

    //! Calculate the stiffness matrix, ionic currents.
    virtual void compute(const double i_stim, NodeMatrix_t & nodes, QuadMatrix_t & quads);

    struct QuadPointStruct {
      double weight;
      Shape_t::FunctionContainer shapeFunctions;
      Shape_t::DerivativeContainer shapeDerivatives;

      QuadPointStruct(double w, const Shape_t & s, NodeMatrix_t & nds);
      };

    typedef std::vector<QuadPointStruct> QuadPointContainer;
    typedef QuadPointContainer::iterator QuadPointIterator;
    typedef QuadPointContainer::const_iterator ConstQuadPointIterator;
    
    //! Access the container of quadrature points
    const QuadPointContainer & quadraturePoints() {return _quadPoints;}

  protected:
    QuadPointContainer _quadPoints;


  };
	
} // namespace voom

#endif // __NobleElem_h__
