// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          William S. Klug
//                University of California Los Angeles
//                   (C) 2004 All Rights Reserved
//
//----------------------------------------------------------------------
//
// $Log$
// Revision 1.2  2005/06/27 04:04:36  klug
// *** empty log message ***
//
// Revision 1.1  2005/05/23 18:05:35  klug
// Initial checkin.
//
//----------------------------------------------------------------------

#include "SimulatedAnnealing.h"

namespace voom {

  // Simulated Annealing Algorithm
  int SimulatedAnnealing::solve(Model * model) 
  {
    _T1 = _T01;
    _T2 = _T02;

    // sched == EXPONENTIAL
    double alpha = pow( _finalTratio, 1.0/_nSteps );

    // sched == LINEAR
    double alpha1 = _T01/_nSteps;
    double alpha2 = _T02/_nSteps;

    if(_size != model->dof() ) resize(model->dof());
    model->getField(*this);
    _compute(model);

    double fBest = _f;
    double fWorst = _f;
    _fSaved = _f;
    _Vector xBest(_x.shape());
    _Vector xWorst(_x.shape());
    xBest = _x;
    xWorst = _x;
    _xSaved = _x;
//     double fAvg = 0;
    std::cout << "\t _fSaved = " << _fSaved << std::endl
	      << "\t _f = " << _f << std::endl
	      << "\t df = " << _f-_fSaved << std::endl
	      << "\t fBest = " << fBest << std::endl
	      << "\t fWorst = " << fWorst << std::endl
// 	      << "\t fAvg = " << fAvg << std::endl
// 	      << "\t 1-fBest/fAvg = " << 1.0-fBest/fAvg << std::endl
	      << "\t _T01 = " << _T1 << std::endl
	      << "\t _T02 = " << _T2 << std::endl;
  
    // seed random number generator
    srand(static_cast<unsigned>(time(0)));

#if defined(ALICANTE)
    _printState( model, "/scratch/klug/SAinit.dat", _x );
#else
    _printState( model, "SAinit.dat", _x );
#endif
    // Begin simulated annealing
    int nPts = _x.size();
    double fSum = _fSaved;
    unsigned accepted = 0;
    for(unsigned step=1; /*1.0-fBest/fAvg > tolerance &&*/ step <= _nSteps; step++ ) {
      accepted = 0;
      for(int pt=0; pt<nPts; pt++) {
	// perform _nSteps metropolis steps at current temperature T
	bool changed = _changeState(model);
	if( changed ) {
	  accepted++;
	  if( _fSaved < fBest ) {
	    fBest = _fSaved;
	    xBest = _xSaved;
	  }
	  if( _fSaved > fWorst ) {
	    fWorst = _fSaved;
	    xWorst = _xSaved;
	  }
	}
	fSum += _fSaved;
      }
//       fAvg = fSum/nPts;
    
      // Reduce temperature
    
      switch( _schedule ) {
      case LINEAR:
	_T1 -= alpha1;
	_T2 -= alpha2;
	break;
      case FAST:
	_T1 = _T01/(1+step);
	_T2 = _T02/(1+step);
	break;
      case EXPONENTIAL:
	_T1 *= alpha;
	_T2 *= alpha;
	break;
      }

      if(step % _printStride == 0) {      
#if defined(ALICANTE)
	_printState( model, "/scratch/klug/SAtemp.dat", _x   );
	_printState( model, "/scratch/klug/SAbest.dat" , xBest );
	_printState( model, "/scratch/klug/SAworst.dat" , xWorst );
#else
	_printState( model, "SAtemp.dat", _x   );
	_printState( model, "SAbest.dat" , xBest );
	_printState( model, "SAworst.dat" , xWorst );
#endif
	// Print acceptance ratio
	std::cout << "Acceptance Ratio: " << accepted << "/" << nPts 
		  << " = " <<(double)accepted/(double)nPts << std::endl
		  << "Temperatures: " << _T1 << ' ' << _T2 << std::endl
		  << "_fSaved: " << _fSaved << std::endl
		  << "_f: " << _f << std::endl
		  << "df = " << _f-_fSaved << std::endl
// 		  << "fAvg: " << fAvg << std::endl
		  << "fBest: " << fBest << std::endl
		  << "fWorst: " << fWorst << std::endl;
      }
    }
#if defined(ALICANTE)
    _printState( model, "/scratch/klug/SAfinal.dat", _xSaved );
    _printState( model, "/scratch/klug/SAbest.dat" , xBest );
    _printState( model, "/scratch/klug/SAworst.dat" , xWorst );
#else
    _printState( model, "SAfinal.dat", _xSaved );
    _printState( model, "SAbest.dat" , xBest );
    _printState( model, "SAworst.dat" , xWorst );
#endif
    //   _x  = _xSaved;
    //   _f = _fSaved;
    _x  = _xSaved  = xBest;
    _f = _fSaved = fBest;

    std::cout << "}" << std::endl;
  }


// --------------------------------------------------------
// Compute a new random trial state for Simulated Annealing
// --------------------------------------------------------
  bool SimulatedAnnealing::_changeState(Model * model)
  {
    double scale=1.0e-1*max(abs(_x));

    // iterate and change field randomly
    for(_Vector::iterator v=_x.begin(); v!=_x.end(); ++v) {
      if(_schedule == FAST ) {
	//
	// generate dv using Cauchy distribution
	// p(dv) = T/(dv^2 + T^2)
	//
	double u = (double)(rand())/RAND_MAX;
	double dv = _T2*tan(M_PI*(u-0.5));
// 	while( false && abs(dv) > u0*(_h[0]+_h[1]+_h[2])/3.0 ) {
// 	  u = (double)(rand())/RAND_MAX;
// 	  dv = _T2*tan(M_PI*(u-0.5));
// 	}
	  //cout << "u = "<<u<<" dv = "<<dv<<endl;
	if( finite(dv) ) *v += dv;
	  
      } else if(_schedule == EXPONENTIAL) {
	//
	// generate dv using gausian distribution
	//
	double u1 = (double)(rand())/RAND_MAX;
	double u2 = (double)(rand())/RAND_MAX;
	double dv = _T2*sqrt(-2.0*log(u1))*cos(2.0*M_PI*u2);
	if( finite(dv) ) *v += dv;
	
      } else {
	//
	// generate dv using uniform distribution
	//
	double dv = 2.0*(double)(rand())/RAND_MAX - 1.0;
	*v += scale*dv;
      }
      
    }
    
    _compute(model);
    
    // decide whether or not to keep the new state
    double df = _f - _fSaved; 
    //       if(df==0.0) std::cout << "df = " << df << std::endl;
      if(_debug)std::cout << "df = " << df << std::endl;
      // metropolis
      double p = ((double)rand())/RAND_MAX;
      if(_debug) {
	std::cout << "p = " << p << std::endl;
	std::cout << "exp(-df/T) = " << exp(-df/_T1) << std::endl;
      }
      
      if( df < 0.0 || p < exp( -df/_T1 ) ) { 
	_xSaved = _x;
	_fSaved = _f;
	return true;
      }
      
      _x = _xSaved;
      _f = _fSaved;

      return false;  
  }

};

