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
    ofstream ofsE, ofsInteractions;
    if (_print)
    {	
      ofsE.open("MontecarloEnergy.dat");
      if (!ofsE) { std::cout << "Cannot open output file MontecarloEnergy.dat" << std::endl;
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
    unsigned accepted = 0, step = 1, pt = 0, j = 0;
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


      if (_print) {
	// Print energy values on file for every sub-step
	ofsE << _f << " " << _fSaved << " " <<  fBest << " " << fWorst << std::endl;
      }
      	
      // Print values of interest (energy, acceptance, Ravg)
      std::cout << "MTS iteration = " << step     << std::endl 
		<< "accepted      = " << accepted << std::endl 
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
      }

    } // End of simulating annealing loop 

 
    if (_print) {
      ofsE.close();
    }

    std::cout << "All done :) " << std::endl;
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



}  // namespace voom

