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

#include "KMCprotein.h"

namespace voom {

  // KMC protein algorithm
  void KMCprotein::solve(uint ComputeNeighInterval, double Rsearch) 
  {
    stringstream OutputFileStream;
    OutputFileStream << "MC_uAVG_" << 0 << ".dat";
    string OutputFileName = OutputFileStream.str();
    ofstream ofsE;
    if ( _Tsched == STEPWISE )
    {	
      ofsE.open(OutputFileName.c_str());
      if (!ofsE) { std::cout << "Cannot open output file " << OutputFileStream.str() << std::endl;
	exit(0); }
    }

    _body->compute(true, false, false);
    _fSaved = _body->energy();
    std::cout << "Initial energy = " << _fSaved << std::endl;
 
    double fBest = _fSaved;
    double fWorst = _fSaved;

    // Temperature
    _T1 = _T01;
    _T2 = _T02;

    // sched == EXPONENTIAL
    double alpha = pow( _FinalTratio, 1.0/_nSteps );

    // sched == LINEAR
    double alpha1 = _T01/_nSteps;
    double alpha2 = _T02/_nSteps;

    // Print initial configuration
    _printProtein->printMaster(0, 1);



    // Begin simulated annealing
    unsigned accepted = 0, step = 1, pt = 0, j = 0;
    unsigned int StepPerInterval = 0, Interval = 0;
    vector<DeformationNode<3>::Point > OriginalLocations;
    vector<DeformationNode<3>* > OriginalHost;
    unsigned int adjust = 0;
    ProteinPotential * Mat = _body->getPotential();
    // std::uniform_real_distribution<double> RealDistribution(0.0, 1.0); // does not include 1.0
    for(pt = 0; pt < _proteinsSize; pt++)
    {
      OriginalLocations.push_back((_proteins[pt]->getHost())->point());
      OriginalHost.push_back( _proteins[pt]->getHost() );
    }

    StepPerInterval = _nSteps/_NT;
    for(step = 1; step <= _nSteps; step++)
    {
      // Compute particles rates in isolation
      this->computeRates();
      // Compute Rmax
      this->computeRmax();

      if (_method == 0) { // Klug MC
	for (int i = 0; i < _proteinsSize; i++) {
	  // p = _randGenerator(RealDistribution);
	  double p = double(rand())/double(RAND_MAX);
	  if (p < _rates[i]) { // accept move
	    _proteins[i]->setHost(_tempHosts[i]);
	    accepted++;
	  }  
	} // loop over all proteins
      } // method 0
      else if (_method == 1) { // KMC adapting Jaime Mariam idea for parallel KMC
	for (int i = 0; i < _proteinsSize; i++) {
	  double r_i = _rates[i]/_rmax;
	  // p = _randGenerator(RealDistribution);
	  double p = double(rand())/double(RAND_MAX);
	  if (p < r_i || r_i > 1.0-1.0e-16) { // accept move
	    _proteins[i]->setHost(_tempHosts[i]);
	    accepted++;
	  }  
	} // loop over all proteins
      } // method 1
      // Update time - correct only for method 1
      double xi = double(rand())/double(RAND_MAX);
      double DeltaTime = -log(xi)/_rmax; // Needed to print following warning
      _time += DeltaTime;
      if ( DeltaTime > 100) { // 100 is an arbitrary number to trigger a warning in the output file - just for debugging
	cout << "DeltaTime = " << DeltaTime << " ; xi = " << xi <<  " ; rmax = " << _rmax << endl;
      }
      _approxTime += 1.0/_rmax;
      
      // Compute body energy
      _body->compute(true, false, false);
      _fSaved = _body->energy();

      if( _fSaved < fBest ) { fBest  = _fSaved; }
      if( _fSaved > fWorst) { fWorst = _fSaved; }

      

      if ( _Tsched == STEPWISE) {
	// Print uAvgSq at each step
	vector<double > uSQavg = this->ComputeUavgSquare(OriginalLocations);
	ofsE << _time << " " << uSQavg[0] << " " << uSQavg[1] << " " << uSQavg[2] << endl;
      }
      	
      // Print values of interest (energy, acceptance, Ravg)
      std::cout << "MTS iteration = " << step        << std::endl 
		<< "accepted      = " << accepted    << std::endl 
		<< "_fSaved       = " << _fSaved     << std::endl 
		<< "fBest         = " << fBest       << std::endl
		<< "fWorst        = " << fWorst      << std::endl
		<< "Temperature   = " << _T1         << std::endl
		<< "time          = " << _time       << std::endl
		<< "approxTime    = " << _approxTime << std::endl;

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
	    std::cout << "Cannot open output file " << OutputFileStream.str() << std::endl;
	    exit(0); 
	  }
	  // Prepare for computing average u square at next temperature or reset proteins to original location
	  for(pt = 0; pt < _proteinsSize; pt++)
	  {
	    // If we do not want to change the original locations then
	    _proteins[pt]->setHost(OriginalHost[pt]);
	    // Otherwise ...
	    // OriginalLocations[pt] = (_proteins[pt]->getHost())->point();
	  }
	  _body->recomputeNeighbors(Rsearch);
	  _printProtein->printMaster(-step, 0);
	  // Reset time to zero
	  _time = 0.0;
	  _approxTime = 0.0;

	  _body->compute(true, false, false);
	  _fSaved = _body->energy();
	  std::cout << "Initial energy = " << _fSaved << std::endl;
	  
	  fBest = _fSaved;
	  fWorst = _fSaved;
	}
	break;
      } // switch

    } // End of simulating annealing loop 

    if (_Tsched == STEPWISE) {
      ofsE.close();
    }

    std::cout << "KMCprotein annealing done :) " << std::endl;
  }



// --------------------------------------------------------------
// Compute a new random trial state based on chosen distribution
// --------------------------------------------------------------
  void KMCprotein::computeRates()
  { 
    for (int i = 0; i < _proteinsSize; i++)
    { 
      ProteinNode * A = _proteins[i];
      DeformationNode<3> * hostA = A->getHost();

      vector<DeformationNode<3> *> NewHosts = _possibleHosts[hostA];
      // std::uniform_int_distribution<int> IntDistribution(0, NewHosts.size()-1);
      // uint j = IntDistribution(generator);
      uint j = rand()%NewHosts.size();
      A->setHost(NewHosts[j]);
      _tempHosts[i] = NewHosts[j]; // If this move is later kept, need to remember where protein went
      
      // Compute body energy
      _body->compute(true, false, false);
      _rates[i] = exp( ( _fSaved - _body->energy())/_T1 );
      A->setHost(hostA); // return protein in position at the beginning of time step
    }
      
  } // computeRates



  vector<double > KMCprotein::ComputeUavgSquare(vector<DeformationNode<3>::Point > & OriginalLocations) 
  {
    double uSQavgTot = 0.0, uSQavgTails = 0.0, uSQavgCenter = 0.0, Increment = 0.0;
    int indTail = 0, indCenter = 0;

    for(uint pt = 0; pt < _proteinsSize; pt++) {
      DeformationNode<3>::Point a = OriginalLocations[pt], b = (_proteins[pt]->getHost())->point();

      // LP: This should be coded in one place only and not everywhere - Body, ProteinPotential, KMCsolver ....
      double DeltaZ = fabs(a(2) - b(2));
      double DeltaZperiodic = fabs(_length - DeltaZ); // 1) if _length < 0 then no periodic BC; 2) assume peridic BC are in Z
      if (DeltaZperiodic < DeltaZ) 
	{ a(2) = 0.0; b(2) = DeltaZperiodic; };
      Increment = pow(tvmet::norm2(a - b), 2.0);

      uSQavgTot += Increment;
      if (b(2) <= _Zmin || b(2) >= _Zmax) {
	uSQavgTails += Increment;
	indTail++;
      } else {
	uSQavgCenter += Increment;
	indCenter++;
      }
	
    }
    
    vector<double > uSQavg(3, 0.0);
    uSQavg[0] = uSQavgTot   /double(_proteinsSize);
    uSQavg[1] = uSQavgTails /double(indTail);
    uSQavg[2] = uSQavgCenter/double(indCenter);
   
    return uSQavg;
  }



}  // namespace voom

