// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          William S. Klug
//                University of California Los Angeles
//                   (C) 2005 All Rights Reserved
//
//----------------------------------------------------------------------
//

#if !defined(__BrownianMotorDynamics_h__)
#define __BrownianMotorDynamics_h__

#include<iostream>
#include<iomanip>
#include<cstring>
#include<string>
#include<blitz/array.h>
#include<random/normal.h>
#include<vector>
#include "Solver.h"
#include "Motor.h"

namespace voom
{

  /*!  A concrete class for a Brownian dynamics solver for a Finite
       Element model with molecular motors.
  */

  class BrownianMotorDynamics
  {
    
  public:
    
    typedef blitz::Array<double,1> Array;
    
    typedef std::vector< Body* > 			BodyContainer;
    typedef std::vector< Body* >::iterator 		BodyIterator;
    typedef std::vector< Body* >::const_iterator 	ConstBodyIterator;

    typedef BrownianNode<2> Node_t;
    
    typedef std::vector< Node_t* > 			NodeContainer;
    typedef std::vector< Node_t* >::iterator 		NodeIterator;
    typedef std::vector< Node_t* >::const_iterator 	ConstNodeIterator;

    typedef std::vector< Constraint* > ConstraintContainer;
    typedef ConstraintContainer::iterator ConstraintIterator;
    typedef ConstraintContainer::const_iterator ConstConstraintIterator;

    typedef Motor<2> MolMot;
    typedef std::vector< MolMot* > MotorContainer;
    typedef MotorContainer::iterator MotorIterator;
    typedef MotorContainer::const_iterator ConstMotorIterator;
    

    BrownianMotorDynamics(NodeContainer & n, MotorContainer & m, int printStride, bool debug=false) 
      :  _nodes(n), _motors(m), _printStride(printStride), _debug(debug) 
    { 
      int dof = 0;
      for(ConstNodeIterator n=_nodes.begin(); n!=_nodes.end(); n++) 
	dof += (*n)->dof();

      std::cout << std::endl
		<< "Constructing BrownianMotorDynamics object with " 
		<< _nodes.size() << " nodes, " << dof << "dof, and " << _motors.size() << "motors."
		<< std::endl
		<< std::endl;

    }

    //! destructor
    virtual ~BrownianMotorDynamics() {}

    //! overloading pure virtual function solve()
    int run(int nSteps, double dt);

    int doMotorHalfStep(double dt);

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

    MotorContainer _motors;

  };
  
}; // namespace voom

#endif // __BrownianMotorDynamics_h__
