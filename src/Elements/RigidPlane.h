// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2010 All Rights Reserved
//
//----------------------------------------------------------------------
//

#if !defined(__RigidPlane_h__)
#define __RigidPlane_h__

#include "Node.h"
#include "Element.h"
#include "VoomMath.h"

namespace voom
{
  
  //! Augmented Lagrange fricitonal contact with a rigid plate
  /*! Enforces contact constraint using the Augmented Lagrange
      approach. Coulomb Friction is included following the algorithm
      listed in chapter 6 of P. Wriggers, "Computational Contact
      Mechanics", 2nd ed., Springer, 2006.
  */
  template< int dim>
  class RigidPlane : public Element
  {
  public:

    // typedefs


    typedef tvmet::Vector<double, dim> VectorND;

    //! Constructor
    /*! Inputs: 

        nodes = a container of pointers to nodes for which contact should
                be enforced;

	k = penalty coefficient (i.e., spring constant)

	Z = the Z position of the plate (in normal direction);

	n = Unit normal to the plane.

        friction = the (Coulomb) friction coefficient (default=0); 
    */
    RigidPlane( const std::vector< DeformationNode<dim> * > & nodes, 
		double k, double Z, 
		const VectorND & n, double friction=0.0 );
 
    //! return the penalty coefficient
    double penaltyCoefficient() const {return _k;}

    //! set the penalty coefficient
    void setPenaltyCoefficient( double k ) { _k = k; }

    //! count up how many nodes are being constrained
    int active() const;

    //! check for active contacts, and update multipliers
    void updateContact();

    //! change the Z position
    void setZ(double Z) { _Z = Z; }

    //! return total force in Z-direction
    double FZ() const {return _FZ;}

    //! return maximum penetration distance
    double penetration() const {return _penetration;}

    //! Compute contact energy and forces
    virtual void compute( bool f0, bool f1, bool f2 );
     
  private:

    //! penalty coefficient
    double _k;
    //! Z position
    double _Z;
    //! total force in Z-direction
    double _FZ;
    //! friction coefficient
    double _mu;
    //! maximum nodal penetration depth
    double _penetration;

    //! normal to the plane
    VectorND _n;

    //! nodes being constrained
    std::vector< DeformationNode<dim> * > _defNodes;
    //! boolean flags indicating nodal penetrations
    std::vector<bool> _active;
    //! contact forces
    blitz::Array<VectorND,1> _forces;
    //! pressure/normal multipliers
    blitz::Array<double,1>   _pressureMultipliers;
    //! friciton/tangent multipliers
    blitz::Array<VectorND,1> _frictionMultipliers;

    //! stick point for frictional contact
    blitz::Array<VectorND,1> _x_stick;

  };

} // namespace voom

#include "RigidPlane.icc"

#endif // __RigidPlane_h__
