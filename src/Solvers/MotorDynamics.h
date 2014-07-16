// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          William S. Klug
//                University of California Los Angeles
//                   (C) 2005 All Rights Reserved
//
//----------------------------------------------------------------------
//

#if !defined(__MotorDynamics_h__)
#define __MotorDynamics_h__

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

  class MotorDynamics
  {
    
  public:
    
    typedef blitz::Array<double,1> Array;
    
    typedef std::vector< Body* > BodyContainer;
    typedef std::vector< Body* >::iterator BodyIterator;
    typedef std::vector< Body* >::const_iterator ConstBodyIterator;

    typedef BrownianNode<2> Node_t;
    
    typedef std::vector< Node_t* > NodeContainer;
    typedef std::vector< Node_t* >::iterator NodeIterator;
    typedef std::vector< Node_t* >::const_iterator ConstNodeIterator;

    typedef Motor<N> MolMot;
    typedef std::vector< MolMot* > MotorContainer;
    typedef typename MotorContainer::iterator MotorIterator;
    typedef typename MotorContainer::const_iterator ConstMotorIterator;
    

    MotorDynamics(MotorContainer & mc,
		     int printStride,
		     bool debug=false) 
      :  _motors(mc), _printStride(printStride), _debug(debug) 
    {
      std::cout << std::endl
		<< "Constructing MotorDynamics object with " 
		<< _motors.size() << " motors "
		<< std::endl;
    }

    //! destructor
    virtual ~MotorDynamics() {}

    //! overloading pure virtual function solve()
    int run(int nSteps, double dt);

    void pushBackBody( Body * bd ) { _bodies.push_back( bd ); }

    void computeAndAssemble( bool f0, bool f1, bool f2 );
     
    double energy() const {return _E;}

  private:	

    double _E;
    int _printStride;

    bool _debug;

    //! container of bodies which contribute to the total model energy
    BodyContainer _bodies;

    //! container of the nodes that represent all of the dof in the model
    MotorContainer _motors;

  };
  
}; // namespace voom

#endif // __MotorDynamics_h__
