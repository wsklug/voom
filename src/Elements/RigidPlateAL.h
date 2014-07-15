// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2010 All Rights Reserved
//
//----------------------------------------------------------------------
//

#if !defined(__RigidPlateAL_h__)
#define __RigidPlateAL_h__

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
  class RigidPlateAL : public Element
  {
  public:

    // typedefs

    //! 3-D deformation node
    typedef DeformationNode<3> DefNode;
    //! Container for node pointers
    typedef std::vector< DefNode* > DefNodeContainer;
    //! Iterator for node pointers
    typedef DefNodeContainer::iterator DefNodeIterator;
    //! Const iterator for node pointers
    typedef DefNodeContainer::const_iterator ConstDefNodeIterator;

    //! Constructor
    /*! Inputs: 

        nodes = a container of pointers to nodes for which contact should
                be enforced;

	k = penalty coefficient (i.e., spring constant)

	Z = the Z position of the plate;

	up = direction flag: is the plate pointing up or down; 

        friction = the (Coulomb) friction coefficient (default=0); 
    */
    RigidPlateAL( const DefNodeContainer & nodes, double k, double Z, 
		  bool up=true, double friction=0.0 );
 
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

    //! is plate pointing up or down?
    bool _up;

    //! nodes being constrained
    DefNodeContainer _defNodes;
    //! boolean flags indicating nodal penetrations
    std::vector<bool> _active;
    //! contact forces
    blitz::Array<Vector3D,1> _forces;
    //! pressure/normal multipliers
    blitz::Array<double,1>   _pressureMultipliers;
    //! friciton/tangent multipliers
    blitz::Array<Vector3D,1> _frictionMultipliers;

    //! stick point for frictional contact
    blitz::Array<Vector3D,1> _x_stick;

    //! normal vector
    Vector3D _normal;

  };

} // namespace voom

#endif // __RigidPlateAL_h__
