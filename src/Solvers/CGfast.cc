// -*- C++ -*-
//----------------------------------------------------------------------
//
//                   Feng Feng and William S. Klug
//                University of California Los Angeles
//                   (C) 2004 All Rights Reserved
//
//----------------------------------------------------------------------
//
// $Log$
// Revision 1.2  2005/08/22 22:28:28  klug
// Model::setField renamed putField
//
// Revision 1.1  2005/05/23 18:07:36  klug
// Initial checkin.
//
//----------------------------------------------------------------------
//
//    Reference:
//
//       Jonathan Richard Shewchuk, "An introduction to the conjugate
//       gradient method without agonizing pain"
//       Edition 5/4:  August 4, 1994
//       School of Computer Science, Carnegie Mellon University
//       Pittsburgh, PA 15213
//
//  http://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.ps
//
/////////////////////////////////////////////////////////////////////////

/*! 
  \file CGfast.cc

  \brief Implementation of a concrete class for nonlinear conjugate
  gradient solver for static equilibrium of a Finite Element model.

*/

#include<iostream>
#include<fstream>
#include<cstdio>
#include<string>

#include "CGfast.h"

using namespace blitz;

namespace voom
{

  //! update f
  void CGfast::_compute(bool f0, bool f1, bool f2) {

    if(_debug) {
      for( Vector_t::const_iterator i=_x.begin(); i!=_x.end(); ++i ) {
	if( std::abs( *i ) > 1.0e5 ) {
	  std::cerr << *i << std::endl
		    << _x << std::endl;
	  _model->print("IncipientDeath");
	  _model->putField( *this );
	  _model->computeAndAssemble( *this, f0, f1, false );
	  _model->print("Death");
	  exit(0);
	}
      }
    }
    _model->putField( *this );
    _model->computeAndAssemble( *this, f0, f1, f2 );
    _model->getField( *this );

    if(f2) {
      _h = 1.0 + abs(_h);
      //std::cout << _h << std::endl;
    }
    return;
  }


  int CGfast::solve(Model * m) {

    _model = m;
    if( _size != _model->dof() ) resize( _model->dof() );
    const blitz::TinyVector< int, 1 > shp = _x.shape();
    Vector_t gradOld(shp); 
    Vector_t searchDir(shp); 
    
    _model->getField( *this );

    // Initial Gradients
    _compute(false,true,false);
    gradOld = _g;

    // compute initial diagonal Hessian for preconditioning
    std::cout << "Computing diagonal preconditioner..." << std::endl;
    _compute(false,false,true);
    std::cout << " done." << std::endl;
    
    Vector_t s(shp);
    //precondition the steepest descent 
    s = -_g/_h;

    // set search direction to preconditioned steepest descent
    searchDir = s;
    
    double deltaNew = -sum(_g*searchDir);

    double delta0 = deltaNew;

    // compute initial residual norm 
    double initialNorm=sqrt(sum(sqr(gradOld)));
    double norm=initialNorm;
    double inftyNorm = max(abs(_g));
    double tolerance = _absTol;
//     double tolerance = std::max(_absTol, _tol*initialNorm);
    //tolerance = std::min( tolerance, initialNorm*1.0e-3 );
    cout << "CG: initial residual " << endl
	 << "       2-norm = "<< norm << endl
	 << "infinity-norm = "<< inftyNorm << endl
	 << "    tolerance = " << tolerance << endl;

    //
    // CG iteration
    //
    int computeCalls=0;
    for(int iter=0,totalIter=0; totalIter<_maxIter; iter++,totalIter++ ) {

      gradOld = _g;
      double normOld = norm;
      ///////////////////////////////////////////////////
      // secant line search
      double delta = sum(sqr(searchDir));
      double dfOld = sum(_g*searchDir);

      // check for descent (PR may not yield this)
      //
      if(dfOld>0) {
	if(_debug)cout << "CG: _g = "<<_g<<endl
		       << "    searchDir = "<<searchDir<<endl
		       << "    f'(0) = "<<dfOld<<endl;
	cout << endl << "    f'(0) > 0. Restarting CG."<<endl<<endl;
	s = -_g/_h;
	searchDir = s;
      }

      double alpha = _sigma;
      _x += alpha*searchDir;
      
      double eps2 = _tolLS*_tolLS;
      int LSiter=0; 
      do {
	_compute(false,true);
	computeCalls++;

	double df = sum(_g*searchDir);

// 	if( df*df < eps2*delta ) break; // LS converged with small derivative

	if(abs(df - dfOld) < _sigma) {
	  std::cout << "Breaking out of line search: df - dfOld = "
		    << df - dfOld 
		    << std::endl;
	  break;
	}
	alpha *= df/(dfOld-df);
	dfOld = df;
	_sigma += alpha; // sum total distance along searchDir
	_x += alpha*searchDir;

	LSiter++;
      } while( LSiter < _maxIterLS && 
// 	       dfOld*dfOld > eps2*delta );
	       alpha*alpha*delta > eps2 ); 
      
      // end line search
      ///////////////////////////////////////////////////

      // compute new gradient and function
      _compute(true,true);
      computeCalls++;

      double deltaOld = deltaNew;
      double deltaMid = -sum(_g*s);

      _compute(false,false,true);

      s = -_g/_h;
      deltaNew = -sum(_g*s);

      // compute norm and check for convergence
      norm = sqrt(sum(sqr(_g)));
      inftyNorm = max(abs(_g));
      // check for convergence
      if ( inftyNorm < tolerance /* || norm < 1.0e-6 HACK!!! */ ) {
	cout << "CG converged with residual 2-norm = "<<norm
	     << " infinity norm = " << inftyNorm 
	     <<" after "<<totalIter<<" iterations "
	     << computeCalls << " computeCalls." <<endl;
	_model->print("cg-converged");
	return 0;
      }	

      ///////////////////////////////////////////////////
      // update search direction
      double beta;
      switch (_cgMethod) {
      case SD:
	// steepest descent
	beta = 0.0;
	break;
      case FR:
	// Fletcher-Reeves
// 	beta = norm/normOld;
// 	beta *= beta;
	beta = deltaNew/deltaOld;
	break;
      case PR:
	// Polak-Ribiere
// 	beta = std::max(0.0, 
// 			sum(_g*(_g-gradOld))/(normOld*normOld) );
	beta = std::max(0.0, (deltaNew-deltaMid)/deltaOld);
	break;
      }
      if(_debug) 
	cout << "CG: beta = "<<beta<<endl
	     << "    searchDirOld="<<searchDir<<endl
	     << "    _g      ="<<_g<<endl;
      
      double nuMax = 0.2;
      double newDotOld = sum(gradOld*_g);
      double normSqrdNew = sum(sqr(_g));
      double nu = fabs( newDotOld )/normSqrdNew;
//    double normSqrdOld = sum(sqr(gradOld));
//    double nu = fabs( newDotOld )/normSqrdOld;
//    double nu = std::abs( newDotOld )/sqrt(normSqrdNew*normSqrdOld);
      if ( /*nu > nuMax ||*/ iter == _restartStride ) {
	iter = 0;
	beta = 0.0;
	searchDir = 0.0;
	cout << "Refreshing descent direction."
	     << endl;
      }
      searchDir *= beta;
      searchDir += s;
      if(_debug) cout << "    searchDirNew="<<searchDir<<endl;
      
      ///////////////////////////////////////////////////

      ///////////////////////////////////////////////////
      // print out info
      if ( totalIter % _printStride == 0 ) {
	cout << "CG: iteration "<<totalIter 
	     << " calls " << computeCalls
	     << setprecision( 8 ) 
	     << " | alpha = "<<alpha
	     << " | sigma = "<<_sigma
	     << " | nu = " << nu 
	     << " | residual = "<<inftyNorm
	     << " | energy = "<<_f<<endl;
// 	char s[20];
// 	sprintf(s,"%04d", totalIter);
	std::string fileName("cgiter");
	_model->print(fileName);
      }
      ///////////////////////////////////////////////////

    }
	
    cout << "CG failed to converge after "<<_maxIter<<" iterations:"<<endl
	 << "\t initialNorm = "<<initialNorm<<"; norm = "<<norm
	 << "; energy = "<<_f<<endl;
    return 1;
  
  }

}; // namespace voom

	
