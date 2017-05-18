// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          William S. Klug
//                University of California Los Angeles
//                   (C) 2005 All Rights Reserved
//
//----------------------------------------------------------------------

#include<iostream>
#include<fstream>
#include<cstdio>
#include<string>

#include "ViscousRelaxation.h"

using namespace blitz;

namespace voom
{

//   //! update f
//   double ViscousRelaxation::_computeFunction( Vector_t & x ) {
//     for( Vector_t::const_iterator i=x.begin(); i!=x.end(); ++i ) {
//       if( std::abs( *i ) > 1.0e5 ) {
// 	std::cerr << *i << std::endl
// 		  << x << std::endl;
// 	_model->print("IncipientDeath");
// 	blitz::cycleArrays( x, _x );
// 	_model->putField( *this );
// 	_model->computeAndAssemble( *this, true, true, false );
// 	_model->print("Death");
// 	exit(0);
//       }
//     }
//     blitz::cycleArrays( x, _x );
//     double temp1=_f;
//     _model->putField( *this );
//     _model->computeAndAssemble( *this, true, false, false );
//     blitz::cycleArrays( x, _x );    
//     double temp2 =_f;
//     _f = temp1;
//     return temp2;
//   }

  //! update gradf
  void ViscousRelaxation::_computeGradient( Vector_t & x, Vector_t & grad ) {
    for( Vector_t::const_iterator i=x.begin(); i!=x.end(); ++i ) {
      if( std::abs( *i ) > 1.0e5 ) {
	std::cerr << *i << std::endl
		  << x << std::endl;
	_model->print("End");
	_model->computeAndAssemble( *this, true, true, false );
	exit(0);
      }
    }
    
    blitz::cycleArrays( x, _x );
    blitz::cycleArrays( grad, _gradf );
    _model->putField( *this );
    _model->computeAndAssemble( *this, false, true, false );
    blitz::cycleArrays( x, _x );
    blitz::cycleArrays( grad, _gradf );
    return;
  }

  int ViscousRelaxation::solve(Model * m)
  {
//     _model->compute(true,true,false);
//     _f = _model->getEnergy();
//     _model->getPositions( _x );
//     _model->getResidual( _gradf );

    _model = m;
    if( _size != _model->dof() || _x.size() != _model->dof() ) resize( _model->dof() );

    assert( _x.size() == _model->dof() );
    _model->getField( *this );


    //
    // compute initial residual norm 
    //
    _model->computeAndAssemble( *this, false, true, false );
    double initialNorm=sqrt(sum(sqr(_gradf)));
    double norm=initialNorm;
    double tolerance = std::max(_absTol, _tol*initialNorm);
    if(_debug) {
      cout << "VR: initial residual norm = "<< norm //<< endl
	   << "    tolerance = " << tolerance << endl;
    }
    //
    // main loop
    //
    for(unsigned iter=0,totalIter=0; totalIter<_maxIter; iter++,totalIter++ ) {
      
      _x -= _dt*_gradf;
      _model->putField( *this );
      _model->computeAndAssemble( *this, true, true, false );
 
      double normOld = norm;
      norm = sqrt(sum(sqr(_gradf)));
      if ( norm < tolerance /* || norm < 1.0e-6 HACK!!! */ ) {
	cout << "VR converged with residual norm = "<<norm
	     << " and energy = "<<_f
	     <<" after "<<totalIter<<" iterations."<<endl;
	_model->print("vr-converged");
	return 0;
      } else if (  totalIter > 0 && totalIter % _printStride == 0) {
	cout << "VR: iteration "<<totalIter 
	     << setprecision( 16 ) 
	     << " | residual norm = "<<norm
	     << " | energy = "<<_f<<endl;
	char s[20];
	sprintf(s,"%04d", totalIter);
	std::string fileName(s);
	_model->print(fileName);
	
      }
      
    }
    
    cout << "VR failed to converge after "<<_maxIter<<" iterations:"<<endl
	 << "\t initialNorm = "<<initialNorm<<"; norm = "<<norm
	 << "; energy = "<<_f<<endl;
    return 1;

  }

}; // namespace voom

	
