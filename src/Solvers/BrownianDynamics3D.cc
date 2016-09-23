// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          William S. Klug
//                University of California Los Angeles
//                   (C) 2005 All Rights Reserved
//
//----------------------------------------------------------------------
//

#include "BrownianDynamics3D.h"
#include "VoomMath.h"

namespace voom {

  int BrownianDynamics3D::run( int nSteps, double dt) {
    
    //
    // main loop
    //
    double t=0;
    for(int step=0; step < nSteps; step++ ) {
    
      // compute only Brownian forces at initial positions
      computeAndAssemble( false, false, true );
      
      // do a half time-step
      for( NodeIterator n=_nodes.begin(); n!=_nodes.end(); n++ ) {	

	// get force and drag
	const Node_t::Point & f = (*n)->force();
	const Node_t::Matrix & D = (*n)->drag();

	// mobility = inverse of drag
	Node_t::Matrix M(0.0);
	invert(D,M);
	
	// compute displacement
	Node_t::Point dx;
	dx = -0.5*dt*M*f;
	
	// save initial position as reference
	(*n)->resetPosition();

	// add displacement to node
	for(int i=0; i<(*n)->dof(); i++) (*n)->addPoint(i,dx(i));
	
      }

      // Compute deterministic forces at new half-step positions, but
      // don't update Brownian forces
      computeAndAssemble( true, true, false );

      // Compute velocities and do a full time-step from the initial positions
      for( NodeIterator n=_nodes.begin(); n!=_nodes.end(); n++ ) {	
	
	// get force and drag
	const Node_t::Point & f = (*n)->force();
	const Node_t::Matrix & D = (*n)->drag();

	// mobility = inverse of drag
        Node_t::Matrix M(0.0);
	invert(D,M);
	
	// compute displacement
	Node_t::Point dx;
	dx = -dt*M*f;
	
	// displace from initial position
	(*n)->setPoint( (*n)->position() );
	for(int i=0; i<(*n)->dof(); i++) (*n)->addPoint(i,dx(i));
	
      }
     
      t += dt;

      if ( _printStride > 0 && step % _printStride == 0) {
	// comptue energy and force and print stuff out
	computeAndAssemble( true, true, false );
	std::cout << "BrownianDynamics: step "<< step
		  << std::setprecision( 16 ) 
		  << " | energy = "<<_E<< std::endl;
	
	char s[20];
	for(int i=0; i<_bodies.size(); i++) {
	  std::cout << "Printing body " << i << std::endl;
	  sprintf(s,"body%d-step%d",i, step);
	  _bodies[i]->print(s);	  
	}
	std::cout << "Printed bodies." << std::endl;
      }
      
    }
    
    computeAndAssemble(true,true,false);
    std::cout << "End Brownian dynamics run."<<std::endl
 	      << "         Energy = "<<_E<<std::endl;
    return 1;
    
  }


  //! Do mechanics, assemble, send data to solver
  void BrownianDynamics3D::computeAndAssemble(bool f0, bool f1, bool f2) 
  {
    if( _constraints.size() == 0 && _bodies.size() == 0 ) {
      std::cout << std::endl
		<< "BrownianDynamics::computeAndAssemble has no bodies or constraints to compute." 
		<< std::endl
		<< std::endl;
    }

    // zero out all forces and stiffness in nodes before computing bodies
    for(NodeIterator n=_nodes.begin(); n!=_nodes.end(); n++) {

      for(int i=0; i<(*n)->dof(); i++) (*n)->setForce(i,0.0);

      (*n)->setMobility( Node_t::Matrix(0.0) );
      (*n)->setDrag( Node_t::Matrix(0.0) );
    }

    // Predictor/corrector approach for constraint
    for(ConstraintIterator c=_constraints.begin(); c!=_constraints.end(); c++) {
      (*c)->predict();
    }

    // compute bodies
    for(BodyIterator b=_bodies.begin(); b!=_bodies.end(); b++) {      
      (*b)->compute( f0, f1, f2);
    }

    // Predictor/corrector approach for constraint
    for(ConstraintIterator c=_constraints.begin(); c!=_constraints.end(); c++) {
      (*c)->correct();
    }

    // assemble
    _E = 0.0;
    for(ConstBodyIterator b=_bodies.begin(); b!=_bodies.end(); b++) {
      if(f0) _E += (*b)->energy();
    }

#ifdef WITH_MPI
    if(f0) {
      double myf=_E;
      MPI_Allreduce(&myf, &(_E), 1, MPI_DOUBLE, 
		    MPI_SUM, MPI_COMM_WORLD);
    }
#endif 

  } // end BrownianDynamics::computeAndAssemble()

  //! check consistency of derivatives
  bool BrownianDynamics3D::checkConsistency(bool verbose) {

  }


} // namespace voom
