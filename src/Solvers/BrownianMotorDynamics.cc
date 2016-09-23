// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          William S. Klug
//                University of California Los Angeles
//                   (C) 2005 All Rights Reserved
//
//----------------------------------------------------------------------
//

#include "BrownianMotorDynamics.h"

namespace voom {

  int BrownianMotorDynamics::run( int nSteps, double dt) {
    double t = 0.0;
    //
    // main loop
    //
    for(int step=0; step < nSteps; step++ ) {
      for( NodeIterator n=_nodes.begin(); n!=_nodes.end(); n++ ) {
      // zero out all forces and stiffness in nodes //
	
	for(int i=0; i<(*n)->dof(); i++) (*n)->setForce(i,0.0);
	(*n)->setMobility( Node_t::Matrix(0.0) );
	(*n)->setDrag( Node_t::Matrix(0.0) );
      }
      if ( _printStride > 0 && step % _printStride == 0) {
	// comptue energy and force and print stuff out
	computeAndAssemble( true, true, false );
	std::cout << "BrownianMotorDynamics: step "<< step
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
      else {
	// compute only forces
	computeAndAssemble( false, true, false );
      }
      
    // get force and mobility
      for( NodeIterator n=_nodes.begin(); n!=_nodes.end(); n++ ) {
	const Node_t::Point & f = (*n)->force();
	const Node_t::Matrix & D = (*n)->drag();
	
	Tensor2D M(0.0);
	double det;
	det = D(0,0)*D(1,1)-D(0,1)*D(1,0);
	M(0,0) = D(1,1)/det;
	M(0,1) = -D(0,1)/det; 
	M(1,0) = -D(1,0)/det; 
	M(1,1) = D(0,0)/det; 
     
	// compute displacement
	Node_t::Point dx;
	dx = -dt*M*f;//*****Added a minus sign*****

	// add displacement to node
	for(int i=0; i<(*n)->dof(); i++) (*n)->addPoint(i,dx(i));

      }

      for( MotorIterator m=_motors.begin(); m!=_motors.end(); m++) {
	if((*m)->isAttached()) (*m)->stepMotor();
      }
      t += dt;
    }
    
//     std::cout << "End Brownian dynamics run."<<std::endl
// 	      << "         Energy = "<<_E<<std::endl;
    return 1;
 
  }

  int BrownianMotorDynamics::doMotorHalfStep(double dt) {
    for (MotorIterator m=_motors.begin(); m!=_motors.end(); m++) {
      (*m)->setTimeStep(dt/2.0);
      (*m)->stepMotor();
      (*m)->setTimeStep(dt);
    }
  }


  //! Do mechanics, assemble, send data to solver
  void BrownianMotorDynamics::computeAndAssemble(bool f0, bool f1, bool f2) 
  {
    if( _constraints.size() == 0 && _bodies.size() == 0 ) {
      std::cout << std::endl
		<< "BrownianMotorDynamics::computeAndAssemble has no bodies or constraints to compute." 
		<< std::endl
		<< std::endl;
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

  } // end BrownianMotorDynamics::computeAndAssemble()

  //! check consistency of derivatives
  bool BrownianMotorDynamics::checkConsistency(bool verbose) {

  }


} // namespace voom
