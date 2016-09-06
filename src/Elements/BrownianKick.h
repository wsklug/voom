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
// Created 2015/11/25 03:06:00  amit
// 
//
//
//----------------------------------------------------------------------

#if !defined(__BrownianKick_h__)
#define __BrownianKick_h__

#include <blitz/array.h>
#include <vector>
#include <random/normal.h>
#include <random/discrete-uniform.h>
#include "Node.h"
#include "Element.h"

namespace voom
{

  //! BrownianKick element provides a layer of random Brownian forces to a Body.
  /*! This class provides a simple way to implement Brownian
      dynamics. The idea is that we can implement viscous
      regularization and add a BrownianKick element to a Body. The
      resultant energy resembles the energy of Brownian dynamics.
   */
 
  class BrownianKick  : public Element
  {
    
  public:

    typedef DeformationNode<3> Node; // nickname for mechanics nodes
    typedef std::vector< Node* > NodeContainer;    

    //! Constructor
    BrownianKick( const NodeContainer &nodes, double Cd, double D, double dt);

    //! Destructor
    ~BrownianKick();
    
    //! Do mechanics on element; compute energy, forces, and/or stiffness.
    void compute(bool f0, bool f1, bool f2);

    //! Set the Brownian random displacements
    void updateParallelKick();

    //! Set the Brownian random displacements
    void updateSerialKick();

    //! Set the diffusion coefficient
    void setDiffusionCoeff( double D ){ _D = D; }

    //! Set drag Coefficient
    void setDragCoefficient( double Cd ){ _Cd = Cd; }

    //! Set time-step
    void setTimeStep( double dt ){ _dt = dt; }

    //! Access the container of nodes
    const NodeContainer& nodes() const { return _nodes; }    
    
  private:
    NodeContainer _nodes;
    double _Cd; //Drag coeff or viscosity coeff
    double _D; //Diffusion coefficient
    double _dt; //time step
    std::vector<Vector3D> _delta_xB; //3-D Brownian displacement

    int _nodeCount;

    ranlib::NormalUnit<double> _rng;
    ranlib::DiscreteUniform<int>* _dis;
  };

} // namespace voom

#endif // __BrownianKick_h__
