// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          William S. Klug
//                University of California Los Angeles
//                   (C) 2005 All Rights Reserved
//
//----------------------------------------------------------------------
//

#if !defined(__BrownianDynamics3D_h__)
#define __BrownianDynamics3D_h__

#include<iostream>
#include<iomanip>
#include<cstring>
#include<string>
#include<blitz/array.h>
#include<random/normal.h>
#include<vector>
#include "Solver.h"

namespace voom
{

  /*!  A concrete class for a Brownian dynamics solver for a Finite
       Element model.
  */

  class BrownianDynamics3D
  {
    
  public:
    
    typedef blitz::Array<double,1> Array;
    
    typedef std::vector< Body* > 			BodyContainer;
    typedef std::vector< Body* >::iterator 		BodyIterator;
    typedef std::vector< Body* >::const_iterator 	ConstBodyIterator;

    typedef BrownianNode<3> Node_t;
    
    typedef std::vector< Node_t* > 			NodeContainer;
    typedef std::vector< Node_t* >::iterator 		NodeIterator;
    typedef std::vector< Node_t* >::const_iterator 	ConstNodeIterator;

    typedef std::vector< Constraint* > ConstraintContainer;
    typedef ConstraintContainer::iterator ConstraintIterator;
    typedef ConstraintContainer::const_iterator ConstConstraintIterator;
    

    BrownianDynamics3D(NodeContainer & n,
		     int printStride,
		     bool debug=false) 
      :  _nodes(n), _printStride(printStride), _debug(debug) 
    { 
      int dof = 0;
      for(ConstNodeIterator n=_nodes.begin(); n!=_nodes.end(); n++) 
	dof += (*n)->dof();

      std::cout << std::endl
		<< "Constructing BrownianDynamics object with " 
		<< _nodes.size() << " nodes and " << dof << " dof."
		<< std::endl
		<< std::endl;

    }

    //! destructor
    virtual ~BrownianDynamics3D() {}

    //! overloading pure virtual function solve()
    int run(int nSteps, double dt);

    void pushBackConstraint( Constraint * c ) { _constraints.push_back( c ); }

    void pushBackBody( Body * bd ) { _bodies.push_back( bd ); }

    //! check consistency
    bool checkConsistency(bool verbose = false);

    void computeAndAssemble( bool f0, bool f1, bool f2 );
     
    double energy() const {return _E;}

  private:	

    double _E;
    int _printStride;

    bool _debug;

    //! container of bodies which contribute to the total model energy
    BodyContainer _bodies;

    ConstraintContainer _constraints;

    //! container of the nodes that represent all of the dof in the model
    NodeContainer _nodes;

  };
  
}; // namespace voom

#endif // __BrownianDynamics_h__
