// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          William S. Klug
//                University of California Los Angeles
//                   (C) 2005 All Rights Reserved
//
//----------------------------------------------------------------------
//

#include "MotorDynamics.h"

namespace voom {

  int MotorDynamics::run( int nSteps, double dt) {

    //
    // main loop
    //
    double t=0;
    for(int step=0; step < nSteps; step++ ) {

      for( MotorIterator m=_motors.begin(); m!=_motors.end(); m++ ) {

	      if ( _printStride > 0 && step % _printStride == 0) {
	// compute energy and force and print stuff out
	computeAndAssemble( true, true, false );
	std::cout << "MotorDynamics: step "<< step
		  << std::setprecision( 16 ) 
		  << " | energy = "<<_E<< std::endl;

	char s[20];
	for(int i=0; i<_bodies.size(); i++) {
	  std::cout << "Printing body " << i << std::endl;
	  sprintf(s,"body%d-step%d",i, step);
	  _bodies[i]->print(s);	  
	}
	std::cout << "Printed bodies." << std::endl;
      } else {
	// compute only forces
	computeAndAssemble( false, true, false );
      }
      
    // get force and mobility
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

      t += dt;


      
    }
    
    std::cout << "End Motor dynamics run."<<std::endl
	      << "         Energy = "<<_E<<std::endl;
    return 1;
 
  }


  //! Do mechanics, assemble, send data to solver
  void MotorDynamics::computeAndAssemble(bool f0, bool f1, bool f2) 
  {
    if( _motors.size() == 0 && _bodies.size() == 0 ) {
      std::cout << std::endl
		<< "MotorDynamics::computeAndAssemble has no motors or constraints to compute." 
		<< std::endl
		<< std::endl;
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

  } // end MotorDynamics::computeAndAssemble()

  //! check consistency of derivatives
  bool MotorDynamics::checkConsistency(bool verbose) {

  }


} // namespace voom
