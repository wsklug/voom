// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          William S. Klug
//                           Luigi Perotti
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

#include "MontecarloProtein.h"

namespace voom {

  // MontecarloProtein Algorithm
  void MontecarloProtein::solve(uint ComputeNeighInterval, double Rsearch) 
  {
    stringstream OutputFileStream;
    OutputFileStream << "MC_uAVG_" << 0 << ".dat";
    string OutputFileName = OutputFileStream.str();
    ofstream ofsE;
    if ( _Tsched == STEPWISE )
    {	
      ofsE.open(OutputFileName.c_str());
      if (!ofsE) { std::cout << "Cannot open output file " << OutputFileStream << std::endl;
	exit(0); }
    }

    _body->compute(true, false, false);
    _f = _body->energy();
    std::cout << "Initial energy = " << _f << std::endl;
 
    _fSaved = _f;
    double fBest = _f;
    double fWorst = _f;

    // Temperature
    _T1 = _T01;
    _T2 = _T02;

    // sched == EXPONENTIAL
    double alpha = pow( _FinalTratio, 1.0/_nSteps );

    // sched == LINEAR
    double alpha1 = _T01/_nSteps;
    double alpha2 = _T02/_nSteps;

    

    // seed random number generator
    int RandSeed = int(time(NULL));
    srand(RandSeed);
    cout << "Seed for random number generator is " << RandSeed << endl;
    


    // Print initial configuration
    _printProtein->printMaster(0, 1);



    // Begin simulated annealing
    unsigned accepted = 0,  accepReq = 0, step = 1, pt = 0, j = 0;
    unsigned int StepPerInterval = 0, Interval = 0;
    vector<DeformationNode<3>::Point > OriginalLocations;
    unsigned int adjust = 0;
    ProteinPotential * Mat = _body->getPotential();
    for(pt = 0; pt < _proteins.size(); pt++)
    {
      OriginalLocations.push_back((_proteins[pt]->getHost())->point());
    }

    StepPerInterval = _nSteps/_NT;
    for(step = 1; step <= _nSteps; step++)
    {
      if (_method == 0) { // change one protein at the time
	for(pt = 0; pt < _proteins.size(); pt++)
	{
	  // perform _size metropolis steps at current temperature T
	  bool changed = changeState();
	  if( changed ) {
	    accepted++;
	    if( _fSaved < fBest ) {
	      fBest = _fSaved;
	    }
	    if( _fSaved > fWorst ) {
	      fWorst = _fSaved;
	    }
	  }
	} // End of sub-loop at constant T
      } // method == 0
      else if (_method == 1) { // change all proteins together
	// perform _size metropolis steps at current temperature T
	bool changed = changeAll();
	if( changed ) {
	  accepted++;
	  if( _fSaved < fBest ) {
	    fBest = _fSaved;
	  }
	  if( _fSaved > fWorst ) {
	    fWorst = _fSaved;
	  }
	}
      } // method == 1



      // Additional steps for method == 2 - imposed pressure
      if (_method == 2) 
      { // try to adjust scaling distance and repeat MC moves at constant pressure
	bool changed = changeEqR();
	if( changed ) {
	  accepReq++;
	  if( _fSaved < fBest ) {
	    fBest = _fSaved;
	  }
	  if( _fSaved > fWorst ) {
	    fWorst = _fSaved;
	  }
	} // if changed loop
	
	for(pt = 0; pt < _proteins.size(); pt++)
	{
	  // perform _size metropolis steps at current temperature T
	  bool changed = changeState();
	  if( changed ) {
	    accepted++;
	    if( _fSaved < fBest ) {
	      fBest = _fSaved;
	    }
	    if( _fSaved > fWorst ) {
	      fWorst = _fSaved;
	    }
	  } // if changed loop
	} // for loop over all proteins
      } // End of addtion for method == 2

      

      if ( _Tsched == STEPWISE) {
	// Print energy values on file for every sub-step
	ofsE << step-(StepPerInterval*Interval) << " " <<  this->ComputeUavgSquare(OriginalLocations) << std::endl;
      }
      	
      // Print values of interest (energy, acceptance, Ravg)
      std::cout << "MTS iteration = " << step     << std::endl 
		<< "accepted      = " << accepted << std::endl 
		<< "accepReq      = " << accepReq << std::endl
		<< "R equilibrium = " << Mat->getEquilibriumR() << std::endl
		<< "_f            = " << _f       << std::endl 
		<< "_fSaved       = " << _fSaved  << std::endl 
		<< "fBest         = " << fBest    << std::endl
		<< "fWorst        = " << fWorst   << std::endl
		<< "Temperature   = " << _T1      << std::endl;

      // Print current configuration
      if (step%_printEvery == 0) {
	_printProtein->printMaster(int(step/_printEvery), 0);
      }

      // Recompute Neighbors
      if (step%ComputeNeighInterval == 0) {
	_body->recomputeNeighbors(Rsearch);
      }

      // Lower temperature for next MC iteration
      switch( _Tsched )
      {
      case CONSTANT:
	_T1 = _T01; // it is useless to set it again to T01
	_T2 = _T02; // it is useless to set it again to T02
	break;
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
      case STEPWISE:
	if (step >= StepPerInterval*(Interval+1))
	{
	  _T1 += (_T02-_T01)/double(_NT-1);
	  Interval++;
	  ofsE.close();
	  OutputFileStream.str(string());
	  OutputFileStream << "MC_uAVG_" << Interval << ".dat";
	  OutputFileName = OutputFileStream.str();
	  ofsE.open(OutputFileName.c_str());
	  if (!ofsE) { 
	    std::cout << "Cannot open output file " << OutputFileStream << std::endl;
	    exit(0); 
	  }

	  // Prepare for computing average u square at next temperature
	  for(pt = 0; pt < _proteins.size(); pt++)
	  {
	    OriginalLocations[pt] = (_proteins[pt]->getHost())->point();
	  }

	}
	break;
      }

      if (_T1 < _resetT /* && adjust == 0 */) {
	_body->resetEquilibrium();
	adjust++;
      }

    } // End of simulating annealing loop 

 
    if (_Tsched == STEPWISE) {
      ofsE.close();
    }

    std::cout << "MontecarloProtein annealing done :) " << std::endl;
  }



// --------------------------------------------------------------
// Compute a new random trial state based on chosen distribution
// --------------------------------------------------------------
  bool MontecarloProtein::changeState()
  { 
      uint i = rand()%_proteins.size();
      
      ProteinNode * A = _proteins[i];
      DeformationNode<3> * hostA = A->getHost();
      vector<DeformationNode<3> *> NewHosts = _possibleHosts[hostA];
 
      uint j = rand()%NewHosts.size();
      A->setHost(NewHosts[j]);

      
      // Compute body energy
      _body->compute(true, false, false);
      _f = _body->energy();
    
      // Decide whether or not to keep the new state
      double df = _f - _fSaved; 
       
      // metropolis
      double p = double(rand())/double(RAND_MAX);
      // cout << "df = " << df << endl;
      if( df < 0.0 || p < exp( -df/_T1 ) )
      { 
	_fSaved = _f;
	return true;
      }
      else {
	A->setHost(hostA);
	return false;  
      }
  } // changeState



  bool MontecarloProtein::changeAll()
  { 
    uint Psize = _proteins.size();
    std::map<uint, DeformationNode<3> *> ProteinsChanged;

    for(uint pt = 0; pt < Psize; pt++) {
      uint i = rand()%Psize;

      ProteinNode * A = _proteins[i];
      DeformationNode<3> * hostA = A->getHost();
      ProteinsChanged.insert(make_pair(i, hostA));
      vector<DeformationNode<3> *> NewHosts = _possibleHosts[hostA];
 
      uint j = rand()%NewHosts.size();
      A->setHost(NewHosts[j]);
    } // change position of all proteins at once
      
    // Compute body energy
    _body->compute(true, false, false);
    _f = _body->energy();
    
    // Decide whether or not to keep the new state
    double df = _f - _fSaved; 
    
    // metropolis
    double p = double(rand())/double(RAND_MAX);
    // cout << "df = " << df << endl;
    if( df < 0.0 || p < exp( -df/_T1 ) ) { 
      _fSaved = _f;
      return true;
    }
    else {
      for(map<uint, DeformationNode<3> *>::iterator mP = ProteinsChanged.begin(); 
	  mP != ProteinsChanged.end(); mP++) {
	_proteins[ mP->first ]->setHost( mP->second ); 
      }
      return false;  
    }
  } // changeAll



  bool MontecarloProtein::changeEqR()
  {
    ProteinPotential * Mat = _body->getPotential();
    double CurrentReq = Mat->getEquilibriumR();
    double NewReq = ( (double(rand())/double(RAND_MAX) - 0.5) * 0.1 + 1.0)*CurrentReq;
    
    Mat->setEquilibriumR(NewReq);
    
    // Compute body energy
    _body->compute(true, false, false);
    _f = _body->energy();
    
    // Decide whether or not to keep the new state
    double df = _f - _fSaved; 
    
    // metropolis
    double p = double(rand())/double(RAND_MAX);
    if( df < 0.0 || p < exp( -df/_T1 ) )
    { 
      _fSaved = _f;
      return true;  
    }
    else {
      Mat->setEquilibriumR(CurrentReq);
      return false;  
    }
  }



  double MontecarloProtein::ComputeUavgSquare(vector<DeformationNode<3>::Point > & OriginalLocations) 
  {
    double uSQavg = 0.0;

    uint Psize = _proteins.size();
    for(uint pt = 0; pt < Psize; pt++) {
      DeformationNode<3>::Point b = (_proteins[pt]->getHost())->point();
      uSQavg += pow(tvmet::norm2(OriginalLocations[pt]-b), 2.0);
    }
    
    uSQavg = uSQavg/double(Psize);

    return uSQavg;
  }



}  // namespace voom

