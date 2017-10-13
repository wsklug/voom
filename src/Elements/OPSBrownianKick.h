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

#if !defined(__OPSBrownianKick_h__)
#define __OPSBrownianKick_h__

#include <blitz/array.h>
#include <vector>
#include <random/normal.h>
#include <random/discrete-uniform.h>
#include "Node.h"
#include "Element.h"

namespace voom
{

  //! OPSBrownianKick element provides a layer of random Brownian forces to a Body.
  /*! This class provides a simple way to implement Brownian
      dynamics. The idea is that we can implement viscous
      regularization and add a OPSBrownianKick element to a Body. The
      resultant energy resembles the energy of Brownian dynamics.
   */
 
  class OPSBrownianKick  : public Element
  {
    
  public:

    typedef std::vector< OPSNode* > NodeContainer;

    //! Constructor
        OPSBrownianKick(const NodeContainer &opsNodes, double v);

	//! Another constructor
    OPSBrownianKick( const NodeContainer &nodes, double v, double cut);

    //! Destructor
    ~OPSBrownianKick();
    
    void brownianStep();

    //! Do mechanics on element; compute energy, forces, and/or stiffness.
    void compute(bool f0, bool f1, bool f2);

    //! Set the Brownian random displacements
    void updateParallelKick();
    
    //!Set the Brownian random displacements projected to surface of a sphere
    void updateProjectedKick();
    
    //!Set rigid rotation kicks
    void updateRotationKick();

	//!Set truncated projected kicks
	void truncatedProjectedKick();

	//!Set 2D kicks
	void update2DKick();

	//!Set 1D kicks
	void update1DKick();
    
    //! Set the Brownian random displacements
    void updateSerialKick();

    //! Set the variance of the random variable
    void setCoefficient(double v){ _coefficient = 1.41421356237*v; }

    //! Access the container of nodes
    const NodeContainer& nodes() const { return _nodes; }    

	//!Get Kick Stats
	std::vector<double> getKickStats();
    
  private:
    NodeContainer _nodes;
    double _coefficient;
	double _cutOff;
    std::vector<Vector3D> _randomKick; //3-D Brownian displacement
    std::vector<Vector3D> _prevX;

    int _nodeCount;

    ranlib::NormalUnit<double> _rng;
    ranlib::DiscreteUniform<int>* _dis;
  };

} // namespace voom

#endif // __OPSBrownianKick_h__
