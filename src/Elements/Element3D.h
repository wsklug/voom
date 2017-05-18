// -*- C++ -*-
//----------------------------------------------------------------------
//
//                           Luigi Perotti
//                University of California Los Angeles
//                   (C) 2006 All Rights Reserved
//
//----------------------------------------------------------------------

/*! 
  \file BrickElement.h

  \brief 3D FE

*/

#if !defined(__Element3D_h__)
#define __Element3D_h__

#include <vector>
#include <cstdio>
#include <ctime>

#include "Element.h"
#include "Node.h"
#include "Quadrature.h"
#include "Shape.h"
#include "Material.h"
#include "VoomMath.h"

using namespace std;

namespace voom
{

  class Element3D : public Element
  {

  public:

    typedef DeformationNode<3> Node; // nickname for mechanics nodes
  
    typedef std::vector< Node* > NodeContainer;
    typedef std::vector< Node* >::iterator NodeIterator;
    typedef std::vector< Node* >::const_iterator ConstNodeIterator;

    //! virtual destructor
    virtual ~Element3D() {;}

    //! Constructor

    Element3D( const NodeContainer & Nodes,
	       Material * Mat,
	       Quadrature<3> & Quad,
	       Shape<3> & Sh,
	       double k = -1.0);

    struct QuadPointStruct 
    {
      double	 weight;
      Material   *material;
      Shape<3>::FunctionContainer shapeFunctions;
      Shape<3>::DerivativeContainer shapeDerivatives;
      QuadPointStruct(double w, Material * Mat, Shape<3> & Sh, const NodeContainer & Nodes);
    };

    typedef std::vector<QuadPointStruct> QuadPointContainer;
    typedef QuadPointContainer::iterator QuadPointIterator;
    typedef QuadPointContainer::const_iterator ConstQuadPointIterator;

    const QuadPointContainer & quadPoints() const {return _quadPoints;}
    const NodeContainer & nodes() const {return _nodes;}

    //! Do mechanics on element; compute energy, forces, and/or stiffness.
    virtual void compute(bool f0, bool f1, bool f2);
    vector<pair<Vector3D, vector<double > > > invariants(int &);

    void reset();
		
    //! access
    double volume() const { return _volume; }
    double strainEnergy() const { return _strainEnergy; }

    //
    //   data
    //
  private:

    NodeContainer _nodes;
    QuadPointContainer _quadPoints;
    double _strainEnergy;
    double _volume;
    double _k;
    Quadrature<3> * _quad;
    Shape<3> * _sh;

  };  // end of class
} // namespace voom



#endif // __Element3D_h__
