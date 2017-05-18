// -*- C++ -*-
/*! 
  \file FiniteContraction.h

*/

#if !defined(__FiniteContraction_h__)
#define __FiniteContraction_h__

#include <blitz/array.h>
#include <tvmet/Matrix.h>
#include <tvmet/Vector.h>
#include <vector>
#include <cstdio>
#include <ctime>

#include "Element.h"
#include "Node.h"
#include "Math.h"

namespace voom
{

  /*!  Concrete class for a .........
  */
  template<class DefQuadrature_t, class Material_t, class DefShape_t, 
  		   class VoltQuadrature_t, class EpMaterial_t, class VoltShape_t>
  class FiniteContraction 
    : public Element
  {
    
  public:
    // typedefs
    typedef tvmet::Vector<double, 3> Vector3D;
    typedef tvmet::Matrix<double, 3, 3> Matrix3D;

    typedef DeformationNode<3> Node_t;
    typedef ScalarFieldNode<3> Vnode_t;
    


    typedef typename std::vector<Node_t*> NodeContainer;
    typedef typename NodeContainer::iterator NodeIterator;
    typedef typename NodeContainer::const_iterator ConstNodeIterator;
    
    typedef typename std::vector<Vnode_t*> VnodeContainer;
    typedef typename VnodeContainer::iterator VnodeIterator;
    typedef typename VnodeContainer::const_iterator ConstVnodeIterator;

    //! virtual destructor
    virtual ~FiniteContraction() {;}

    //! Constructor
    FiniteContraction( const DefQuadrature_t & defQuad,
		       const Material_t & mat,
		       const NodeContainer & defNodes,
		       const VoltQuadrature_t & vQuad,
		       const EpMaterial_t & epMat,
		       const VnodeContainer & voltNodes
		      );

  public:
   
    //! Do mechanics on element; compute energy, forces, and/or stiffness.
    virtual void compute(bool f0, bool f1, bool f2);
    
    virtual void computeEP(bool f0, bool f1, bool f2);



    double strainEnergy() const { return _strainEnergy; }

    Matrix3D cauchyStress() const { return _cauchy; }

    Matrix3D pkStress() const { return _PKStress; }

    double vonMisesStress() const { return _vonMisesStress; }

    //! Access the container of deformation nodes
    const NodeContainer & defNodes() { return _dNodes; }
    
    //! Access the container of deformation nodes
    const VnodeContainer & voltageNodes() { return _vNodes; }
        
    //! compute positions
    Vector3D computePosition(const double s1, const double s2);

    //! check new positions
    void checkPositions();

    //void setElementNum(int ei) { _elNum = ei; }

    struct DefQuadPointStruct 
    {
      double	 weight;
      typename DefShape_t::FunctionContainer defShapeFunctions;
      typename DefShape_t::DerivativeContainer defShapeDerivatives;
      Material_t material;
      double voltage;
      DefQuadPointStruct(double w, const Material_t & m, const DefShape_t & s, const NodeContainer & defnds, double volt);
    };
    
    struct VoltQuadPointStruct 
        {
          double	 weight;
          typename VoltShape_t::FunctionContainer voltShapeFunctions;
          typename VoltShape_t::DerivativeContainer voltShapeDerivatives;
          EpMaterial_t material;
          VoltQuadPointStruct(double w, const EpMaterial_t & epmat, const DefShape_t & d_s, 
        		  const VoltShape_t & v_s, const NodeContainer & defnds);
        };

    typedef typename std::vector<DefQuadPointStruct> DefQuadPointContainer;
    typedef typename DefQuadPointContainer::iterator DefQuadPointIterator;
    typedef typename DefQuadPointContainer::const_iterator ConstDefQuadPointIterator;
    
    //! Access the container of quadrature points
    const DefQuadPointContainer & defQuadraturePoints() const { return _defQuadPoints; }
    DefQuadPointContainer & defQuadraturePoints() { return _defQuadPoints; }
    
    
    typedef typename std::vector<VoltQuadPointStruct> VoltQuadPointContainer;
    typedef typename VoltQuadPointContainer::iterator VoltQuadPointIterator;
    typedef typename VoltQuadPointContainer::const_iterator ConstVoltQuadPointIterator;
    
    const VoltQuadPointContainer & voltQuadraturePoints() const { return _voltQuadPoints; }
    VoltQuadPointContainer & voltQuadraturePoints() { return _voltQuadPoints; }
    
    //
    //   data
    //
  private:
    //! forces on the element
    blitz::Array< Vector3D, 1> _internalForce;

    //! strain energy
    double _strainEnergy;

    // First Piola-Kirchhoff stress
    Matrix3D _PKStress;

    // Cauchy stress
    Matrix3D _cauchy;

    // von Mises stress
    double _vonMisesStress;

    //int _elNum;

    NodeContainer _dNodes;  
    VnodeContainer _vNodes; 
    DefQuadPointContainer _defQuadPoints;
    VoltQuadPointContainer _voltQuadPoints;

  };  // end of class
} // namespace voom

#include "FiniteContraction.icc"

#endif // __FiniteContraction_h__
