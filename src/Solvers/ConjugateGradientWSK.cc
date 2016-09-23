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
// Revision 1.15  2005/10/09 22:49:20  klug
// Cleaned up print statements.
//
// Revision 1.14  2005/08/22 22:28:28  klug
// Model::setField renamed putField
//
// Revision 1.13  2005/08/20 01:28:32  klug
// acinclude.m4
//
// Revision 1.12  2005/05/23 18:03:17  klug
// Enforced consistency with new Solver interface.  Minor tweaks.
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
  \file ConjugateGradientWSK.cc

  \brief Implementation of a concrete class for nonlinear conjugate
  gradient solver for static equilibrium of a Finite Element model.

*/

#include<iostream>
#include<fstream>
#include<cstdio>
#include<string>

#include "ConjugateGradientWSK.h"

using namespace blitz;

namespace voom
{

  //! update f
  double ConjugateGradientWSK::_computeFunction( Vector_t & x ) {
    for( Vector_t::const_iterator i=x.begin(); i!=x.end(); ++i ) {
      if( std::abs( *i ) > 1.0e5 ) {
	std::cerr << *i << std::endl
		  << x << std::endl;
	_model->print("IncipientDeath");
	blitz::cycleArrays( x, _x );
	_model->putField( *this );
	_model->computeAndAssemble( *this, true, true, false );
	_model->print("Death");
	exit(0);
      }
    }
    blitz::cycleArrays( x, _x );
    double temp1=_f;
    _model->putField( *this );
    _model->computeAndAssemble( *this, true, false, false );
    blitz::cycleArrays( x, _x );    
    double temp2 =_f;
    _f = temp1;
    return temp2;
  }

  //! update gradf
  void ConjugateGradientWSK::_computeGradient( Vector_t & x, Vector_t & grad ) {
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

  double ConjugateGradientWSK::_computeAll( Vector_t & x, Vector_t & grad ) {
    if(_debug) {
      for( Vector_t::const_iterator i=x.begin(); i!=x.end(); ++i ) {
	if( std::abs( *i ) > 1.0e5 ) {
	  std::cerr << *i << std::endl
		    << x << std::endl;
	  _model->print("End");
	  _model->computeAndAssemble( *this, true, true, false );
	  exit(0);
	}
      }
    }
    
    blitz::cycleArrays( x, _x );
    blitz::cycleArrays( grad, _gradf );
    double fsaved=_f;
    _model->putField( *this );
    _model->computeAndAssemble( *this, true, true, false );
    blitz::cycleArrays( x, _x );
    blitz::cycleArrays( grad, _gradf );
    double fnew =_f;
    _f = fsaved;
    return fnew;
  }

  int ConjugateGradientWSK::solve(Model * m)
  {
//     _model->compute(true,true,false);
//     _f = _model->getEnergy();
//     _model->getPositions( _x );
//     _model->getResidual( _gradf );

    _model = m;
    if( _size != _model->dof() || _x.size() != _model->dof() ) resize( _model->dof() );
    const blitz::TinyVector< int, 1 > shp = _x.shape();
    Vector_t gradNew(shp); 
    Vector_t gradOld(shp); 
    Vector_t searchDir(shp); 

    assert( _x.size() == _model->dof() );
    _model->getField( *this );

#if defined(ALICANTE)
    ofstream cgout("/scratch/klug/cg.dat");
#else
    ofstream cgout("cg.dat");
#endif

    //
    // Initial Gradients
    //
    _computeGradient(_x, gradOld);

    //
    // compute initial residual norm 
    //
    double initialNorm=sqrt(sum(sqr(gradOld)));
    double norm=initialNorm;
    double tolerance = std::max(_absTol, _tol*initialNorm);
    //tolerance = std::min( tolerance, initialNorm*1.0e-3 );
//     cout << "CG: initial residual norm = "<< norm //<< endl
// 	 << "    tolerance = " << tolerance << endl;

    //
    // set search direction to steepest descent
    //
    searchDir = -1.0*gradOld;
    
    //
    // CG iteration
    //
      for(unsigned iter=0,totalIter=0; totalIter<_maxIter; iter++,totalIter++ ) {

	//
	// line search to find alpha such that _x += alpha*searchDir
	// is a minimum
	//
	double df0 = sum(gradOld*searchDir);

	//
	// check for descent (PR may not yield this)
	//
	if(df0>0) {
	  if(_debug)cout << "CG: gradOld = "<<gradOld<<endl
			 << "    searchDir = "<<searchDir<<endl
			 << "    f'(0) = "<<df0<<endl;
	  cout << endl << "    f'(0) > 0. Restarting CG."<<endl<<endl;
	  searchDir = -1.0*gradOld;
	}

	double alpha = _lineSearch(searchDir);
// 	if(_debug || alpha < 0.0){
// 	  //_debug = false;
// 	  cout << "CG: alpha = " << alpha << endl;
// 	  Vector_t x(_x.shape());
// 	  Vector_t grad(_x.shape());
// 	  ofstream ofs1("lineSearch_f.dat");
// 	  ofstream ofs2("lineSearch_df.dat");
// 	  x = _x;
// 	  for(int i=0; i<100; i++) {
// 	    double a = alpha*(3.0*i/100.0 - 1.0);
// 	    x = _x + a*searchDir;
// 	    double f = _computeFunction(x);
// 	    //_computeGradient(x,grad);
// 	    //double df = sum(sum(grad*searchDir));
// 	    double df = _computeDirectionalD( x, searchDir, grad );
// 	    ofs1 << setprecision(16) << a << '\t' << f << endl;
// 	    ofs2 << setprecision(16) << a << '\t' << df << endl;
// 	  }
// 	  exit(0);
// 	}
	_x += alpha*searchDir;

	//
	// compute new gradient
	//
// 	_computeGradient(_x, gradNew);
	_f = _computeAll(_x, gradNew);

	//
	// compute norm and check for convergence
	//
// 	_f = _computeFunction(_x);
	double normOld = norm;
	norm = sqrt(sum(sqr(gradNew)));
	//////////////////////////////////////////////////////////
// 	cout << "enter..." << endl;
// 	cout << norm << endl;
// 	cout << _tol << endl;
// 	cout << initialNorm << endl;
// 	cout << totalIter << endl;
// 	cout << _printStride << endl;
	
	if ( norm < tolerance /* || norm < 1.0e-6 HACK!!! */ ) {
#ifdef WITH_MPI
	  if(_processorRank==0) 
#endif
	  {
	    cout << "CG converged with residual norm = "<<norm
		 << " and energy = "<<_f
		 <<" after "<<totalIter<<" iterations."<<endl;
	  }	  
	  _model->print("cg-converged");
	  return 0;
	} else if (  totalIter > 0 && totalIter % _printStride == 0) {
#ifdef WITH_MPI
	  if(_processorRank==0)
#endif
	  {
	    cout << "CG: iteration "<<totalIter 
// 		 << setprecision( 16 ) 
		 << " | alpha = "<<alpha
		 << " | residual norm = "<<norm
		 << " | energy = "<<_f<<endl;
	  }
	  char s[20];
	  sprintf(s,"%04d", totalIter);
	  std::string fileName(s);
	  _model->print(fileName);
	}

	//////////////////////////////////////////////////////////
//	cout << "exit..." << endl;

	//
	// update search direction
	//
	double beta;
	switch (_cgMethod) {
	case SD:
	  // steepest descent
	  beta = 0.0;
	  break;
	case FR:
	  // Fletcher-Reeves
	  beta = sum(sqr(gradNew))/sum(sqr(gradOld));
	  break;
	case PR:
	  // Polak-Ribiere
	  beta = std::max(0.0, 
			  sum(gradNew*(gradNew-gradOld))/sum(sqr(gradOld)));
	  break;
	}
	//searchDir = beta*searchDir - gradNew;
	if(_debug) 
	  cout << "CG: beta = "<<beta<<endl
	       << "    searchDirOld="<<searchDir<<endl
	       << "    gradNew     ="<<gradNew<<endl;
//	double nuMax = 10.0;
#if 0
	double nuMax = 0.5;
	double normSqrdOld = sum(sqr(searchDir));
	double normSqrdNew = sum(sqr(beta*searchDir-gradNew));
	double newDotOld = sum(searchDir*(beta*searchDir-gradNew));
// 	double nu = fabs( newDotOld )/normSqrdOld;
	double nu = std::abs( newDotOld )/sqrt(normSqrdNew*normSqrdOld);
#else
	double nuMax = 0.2;
	double newDotOld = sum(gradOld*gradNew);
	double normSqrdNew = sum(sqr(gradNew));
 	double nu = fabs( newDotOld )/normSqrdNew;
// 	double normSqrdOld = sum(sqr(gradOld));
//  	double nu = fabs( newDotOld )/normSqrdOld;
// 	double nu = std::abs( newDotOld )/sqrt(normSqrdNew*normSqrdOld);
#endif
	if ( nu > nuMax || iter == _restartStride ) {
	  iter = 0;
	  beta = 0.0;
	  searchDir = 0.0;
#ifdef WITH_MPI
	  if(_processorRank==0)
#endif
	    {
	      cout << "CG: iteration "<< totalIter 
// 		   << setprecision( 16 ) 
		   << " | alpha = "<<alpha
		   << " | residual norm = "<<norm
		   << " | initial norm = "<<initialNorm
		   << " | energy = "<<_f
		   << " | nu = " << nu << " | restarting."
		   << endl;
	    }
	}
	searchDir *= beta;
	searchDir -= gradNew;
	if(_debug) cout << "    searchDirNew="<<searchDir<<endl;

	//
	// old <- new
	//
	cycleArrays(gradOld, gradNew);      
      
      }
    
#ifdef WITH_MPI
      if(_processorRank==0) 
#endif
	{
	  cout << "CG failed to converge after "<<_maxIter<<" iterations:"<<endl
	       << "\t initialNorm = "<<initialNorm<<"; norm = "<<norm
	       << "; energy = "<<_f<<endl;
	}
      return 1;      
  }
  
  double ConjugateGradientWSK::_lineSearch(const Vector_t & dir)
  {
    if( _lsMethod == Secant )  return  _secantLineSearch( dir );
    if( _lsMethod == Wolfe  )  return  _wolfeLineSearch( dir );
    exit(0);
    return 0;
  }

  double ConjugateGradientWSK::_secantLineSearch(const Vector_t & dir)
  {
    bool debug=false;

    Vector_t x(_x.shape());

    Vector_t grad(dir.shape());

    x = _x;
//     double f0 = _computeFunction(x);
//     double df0 = _computeDirectionalD( x, dir, grad );
    double f0 = _computeAll(x,grad);
    double df0 = _directionalD( dir, grad );

    double fOld = f0;
    double dfOld = df0;

    double sigma = 1.0e-8*blitz::max(blitz::abs(_x));
//     double normDir=sqrt(sum(sqr(dir)));
//     double tolerance = _tolLS*blitz::max(blitz::abs(_x));
//     double sigma=10.0*tolerance/normDir;

    x = _x + sigma*dir;
//     double fNew = _computeFunction(x);
//     double dfNew = _computeDirectionalD( x, dir, grad );
    double fNew = _computeAll(x,grad);
    double dfNew = _directionalD( dir, grad );

    int maxIter = _maxIterLS;
    double dalpha = sigma;
    double alpha = dalpha;
    
    int iter;
    double tolerance = std::max(_tolLS*std::abs(df0), _absTol*_absTol);
    for( iter=0; iter<maxIter && std::abs(dfNew) > tolerance; iter++ ) {
//     for( iter=0; iter<maxIter && std::abs(dalpha)*normDir > tolerance; iter++ ) {
      dalpha *= dfNew/(dfOld-dfNew);
      x += dalpha*dir;
      alpha += dalpha;
      fOld = fNew;
//       fNew = _computeFunction(x);
      fNew = _computeAll(x,grad);
      dfOld = dfNew;
//       dfNew = _computeDirectionalD( x, dir, grad );
      dfNew = _directionalD( dir, grad );
//       if( dfNew == dfOld) {
// 	cout << "Something's wrong.  Directional derivative is not changing in secant line search." << std::endl
// 	     << "    fOld = " << setw(16) << fOld 
// 	     << "    fNew = " << setw(16) << fNew << endl
// 	     << "   dfOld = " << setw(16) << dfOld 
// 	     << "   dfNew = " << setw(16) << dfNew << endl
// 	     << "   alpha = " << setw(16) << alpha 
// 	     << "  dalpha = " << setw(16) << alpha << endl;
// 	exit(0);
//       }
      if( debug || _debug ) {
	cout << dalpha << "  " << alpha << "  " 
	     << fOld << "  " << dfOld << "  " << fNew << "  " << dfNew << endl;
      }
    }
    if( debug || _debug) cout << "_secantLineSearch: after " << iter 
			      << " iterations, dfNew = " << dfNew 
			      << ", df0 = " << df0 
			      << ", tolerance = "<<tolerance
			      << ", dalpha = "<<dalpha
			      << ", alpha = "<<alpha<<endl;
    return alpha;
  
  }

  //
  // Line search algorithm satisfying strong Wolfe conditions
  //
  double ConjugateGradientWSK::_wolfeLineSearch(const Vector_t & dir)
  {
    bool debug=false;
    //   const double c1 = 1.0e-4, c2 = 0.9;
    const double c1 = _c1, c2 = _c2;

    Vector_t x(_x.shape());
    x = _x;

    double eps = 1.0e-6;
    double 
      alphaOld    = 0.0, 
      alphaNew    = -1.0e-1, 
      alphaMax    = 1.0e2,
      alphaFactor = 100.0; 

    double f0 = _f;

    Vector_t grad(dir.shape());
    //   _computeGradient(x, grad);
    //   if(_debug) {
    //     cout << "lineSearch: x = "<<x<<endl;
    //     cout << "lineSearch: grad = "<<grad<<endl;
    //     cout << "lineSearch: dir = "<<dir<<endl;
    //   }
    //   double df0 = sum(sum(grad*dir));
//     double df0 = _computeDirectionalD( x, dir, grad );
    _computeGradient(x, grad);
    double df0 = _directionalD( dir, grad );
    if(df0 > 0 ) 
      cerr << endl<<"_lineSearch: df0 = "<<df0<<"!!!!!!"<<endl<<endl;
    double 
      fOld  = f0, 
      fNew  = f0,
      dfOld = df0,
      dfNew = df0;

    // quadratic interpolation using df0, f0, & f
    x = _x + ( alphaOld + eps )*dir;
    fNew =_computeFunction(x);
  
    double alo  = alphaOld;
    double ahi  = alphaOld + eps;
    double flo  = fOld;
    double dflo = dfOld;
    double fhi  = fNew;
    alphaNew  = alo - dflo*(ahi-alo)*(ahi-alo)/(2.0*(fhi-flo-dflo*(ahi-alo)));
    for ( int i=0; i<10 && alphaNew < 0.0; i++ ) {
      ahi *= 2.0;
      x = _x + ( ahi )*dir;
      fhi =_computeFunction(x);
      alphaNew  = alo - dflo*(ahi-alo)*(ahi-alo)/(2.0*(fhi-flo-dflo*(ahi-alo)));
    }
    if ( alphaNew < 0.0 ) alphaNew = eps;

    if(_debug||debug) {
      cout << setprecision(16)
	   << "_linesearch: choosing starting value of alpha with quadratic interpolation."
	   << " alo = " << alo << endl
	   << " ahi = " << ahi << endl
	   << " flo = " << flo << endl
	   << " dflo = " << dflo << endl
	   << " fhi = " << fhi << endl
	   << " alphaNew = "<<alphaNew<<endl;
    }
    int maxIter = _maxIterLS;
    for ( unsigned i=1; i<=maxIter; i++ ) {
      //
      // evaluate function at new alpha
      //
      x = _x + alphaNew*dir;
      fNew =_computeFunction(x);

      //
      // check conditions
      //
      if ( ( fNew > f0 + c1*alphaNew*df0 ) /*||
	   ( fNew >= fOld && i > 1 )*/         ) {
	if(_debug || debug) 
	  cout<<"_lineSearch: case 1 iter = "<< i 
	      <<" zooming with fNew="<<fNew<<" f0="<<f0
	      <<" df0="<<df0<<endl
	      <<" fOld="<<fOld<<" dfOld="<<dfOld<<endl
	      <<" alphaOld="<<alphaOld<<" alphaNew="<<alphaNew
	      <<endl;
	double alpha = _zoom(c1, c2, f0, df0, dir, x, grad, alphaOld, alphaNew);
	if(_debug || debug) 
	  cout <<" alpha = "<<alpha<<endl;
	cout << "case 1" << endl;
	return(alpha);
      }
      //
      // evaluate derivative at new alpha
      //
      _computeGradient(x,grad);
      dfNew = _directionalD( dir, grad );
      //
      // check conditions
      //
      if ( std::abs(dfNew) <= -c2*df0 ) {
	if(_debug || debug) 
	  cout<<"_lineSearch: case 2 iter = "<< i 
	      <<" fNew = "<< fNew <<" alpha = "<<alphaNew<<endl;
	cout << "case 2" << endl;
	return(alphaNew);      
      }
      if ( dfNew >= 0 ) {
	double alpha = _zoom(c1, c2, f0, df0, dir, x, grad, alphaNew, alphaOld);
	if(_debug || debug) 
	  cout<<"_lineSearch: case 3 iter = "<< i 
	      <<" fNew = "<< fNew <<" alpha = "<<alpha<<endl;
	cout << "case 3" << endl;
	return(alpha);
      }
      //
      // Choose new alpha using quadratic interpolation using 
      // fOld, fNew, & dfOld
      alo  = alphaOld;
      ahi  = alphaNew;
      flo  = fOld;
      dflo = dfOld;
      fhi  = fNew;
      if(_debug||debug) 
	cout << "_linesearch: choosing new alpha with quadratic interpolation."
	     << endl
	     << " alo="<<alo<<" flo="<<flo<<" dflo="<<dflo<<endl
	     << " ahi="<<ahi<<" fhi="<<fhi<<endl;
      alphaOld = alphaNew;
      alphaNew = alo - dflo*(ahi-alo)*(ahi-alo)/(2.0*(fhi-flo-dflo*(ahi-alo)));
      if ( (alphaNew - alphaOld) < 1.0e-3*alphaOld || i > maxIter/2 ) 
	alphaNew = 10.0*alphaOld;

      if(_debug||debug)
	cout << " alphaNew = "<<alphaNew<<endl;
    
      fOld  = fNew;
      dfOld = dfNew;
    }
    double alpha = alphaOld;
    _debug = true;
   if(debug||_debug){
      cout << "line search reached "<<maxIter<<" iterations." << endl;
      cout <<"alpha = "<<alpha<<"alphaOld = "<<alphaOld<<"alphaNew = "<<alphaNew
	   <<endl;
    }
      cout << "case 4" << endl;
    return(alpha);
  }


  //
  // Specialty function used by Wolfe line search algorithm
  //
  double ConjugateGradientWSK::_zoom(const double c1, const double c2, 
					  const double f0, const double df0, 
					  const Vector_t & dir, 
					  Vector_t & x, 
					  Vector_t & grad,
					  double alo, double ahi)
  {
    bool debug=false;
    bool bisection=false;

    int maxIter = _maxIterLS;

    double a = 0.0;
    double flo;
    double dflo;
    double fhi;
    double f;
    double df;

    x = _x + alo*dir;
//     flo =_computeFunction(x);
//     dflo = _computeDirectionalD( x, dir, grad );
    flo =_computeAll(x,grad);
    dflo = _directionalD( dir, grad );

    x = _x + ahi*dir;
    fhi =_computeFunction(x);

    for( unsigned i=0; i<maxIter; i++) {
      if( std::abs(ahi-alo)/(ahi+alo) < 1.0e-8 || 
	  (ahi+alo)<1.0e-16 ) return((ahi+alo)/2.0);
      //
      // interpolate to find a trial step length alpha between lo and
      // hi.
      //
      if(bisection || i > maxIter/2 ) {
	a = (alo + ahi)/2.0;
      } else { 
	// quadratic interpolation using fhi, flo, & dflo
	a = alo - dflo*(ahi-alo)*(ahi-alo)/(2.0*(fhi-flo-dflo*(ahi-alo)));
      }
      //
      // evaluate function
      //
      x = _x + alo*dir;
      flo =_computeFunction(x);
      x = _x + a*dir;
//       f =_computeFunction(x);
      f =_computeAll(x,grad);
      if(_debug || debug)
	cout << "_zoom: iter = "<<i<<" a="<<a<<" in ["<<alo<<","<<ahi<<"]"<<endl
	     << "\t f="<<f<<" f0="<<f0<<" flo="<<flo<<endl;
      //
      // Check stuff
      //
      if ( ( f > f0 + c1*a*df0 ) || ( f >= flo ) ) {
	ahi = a;
	fhi = f;
      } else {
	//
	// evaluate derivative
	//
	
// 	_computeGradient(x,grad);
	df = _directionalD( dir, grad );

	if(_debug || debug)
	  cout<< "\t df="<<df<<" df0="<<df0<<endl;
	//
	// Check stuff
	//
	if ( std::abs(df) <= -c2*df0 ) {
	  return(a);
	}
	if ( df*(ahi-alo) >= 0 ) {
	  ahi = alo;
	  fhi = flo;
	}
	alo = a;
	flo = f;
	dflo = df;
      }
    }
    //_debug = true;
    if(_debug || debug)
      cout << "_zoom reached "<<maxIter<<" iterations:"
	   << "a="<<a<<" in ["<<alo<<","<<ahi<<"]"
	   <<endl;
    return(a);
  }
}; // namespace voom

	
