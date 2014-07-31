// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2010 All Rights Reserved
//
//----------------------------------------------------------------------
//

#if !defined(__RigidHemisphereAL_h__)
#define __RigidHemisphereAL_h__

#include "Node.h"
#include "Element.h"
#include "VoomMath.h"

namespace voom
{
  
  //! Augmented Lagrange fricitonal contact with a rigid hemisphere
  /*! Enforces contact constraint using the Augmented Lagrange
      approach. Coulomb Friction is included following the algorithm
      listed in chapter 6 of P. Wriggers, "Computational Contact
      Mechanics", 2nd ed., Springer, 2006.
  */
  class RigidHemisphereAL : public Element
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

        1) a container of pointers to nodes for which contact should
           be enforced;

	2) penalty coefficient (i.e., spring constant)

	3) radius of the hemisphere; 

	4) 3-D position of the center of
           the hemisphere; 

        5) friction coefficient (default=0); 
    */
    RigidHemisphereAL( const DefNodeContainer & nodes, 
		       double k, double R, Vector3D xc, double friction=0.0 );

    //! return the penalty coefficient
    double penaltyCoefficient() const {return _k;}

    //! set the penalty coefficient
    void setPenaltyCoefficient( double k ) { _k = k; }

    //! count up how many nodes are being constrained
    int active() const;

    //! return whether a node is penetrating or not
    bool active_node(int i) { 
     if(i>_active.size()){
       std::cout<<"RgidiHemisphereAL::active_node::Error Node number out of range"<<std::endl;
       exit(0);}
     return _active[i];
    }

    //! check for active contacts, and update multipliers
    void updateContact();

    //! return total force in Z-direction
    double FZ() const {return _FZ;}

    //! return maximum penetration distance
    double penetration() const {return _penetration;}

    //! Compute contact energy and forces
    virtual void compute( bool f0, bool f1, bool f2 );
     
  private:

    //! penalty coefficient
    double _k;
    //! hemisphere radius
    double _R;
    //! total force in Z-direction
    double _FZ;
    //! friction coefficient
    double _mu;
    //! maximum nodal penetration depth
    double _penetration;
    //! curvature center of hemisphere
    Vector3D _xc;
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
    //! normal vectors (at contact point)
    blitz::Array<Vector3D,1> _normals;

  };

} // namespace voom

#endif // __RigidHemisphereAL_h__
