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

#include "MontecarloTwoStages.h"

namespace voom {

  // MontecarloTwoStages Algorithm
  int MontecarloTwoStages::solve(Model * model) 
  {
    ofstream ofsE, ofsV, ofsInteractions;
    if (_print)
    {	
      ofsE.open("MontecarloEnergy.dat");
      ofsV.open("MontecarloVariables.dat");
      if (!ofsE) { std::cout << "Cannot open output file MontecarloEnergy.dat" << std::endl;
	exit(0); }
      if (!ofsV) { std::cout << "Cannot open output file MontecarloVariables.dat" << std::endl;
	exit(0); }
    }

    vector<int > Interactions(2, 1);
    if (_ferroMagnetic)
    {
      ofsInteractions.open("MontecarloInteractions.dat");
      if (!ofsInteractions) { std::cout << "Cannot open output file MontecarloInteractions.dat" << std::endl;
	exit(0); 
      }
    }


    // Solve the problem at the beginning to release pre-existing stresses
    // Then start Montecarlo on the relaxed structure
    _solver->solve(model);

    _f = 0.0;
    Model::BodyContainer bdc = model->bodies();
    for(Model::ConstBodyIterator bdcIt = bdc.begin(); bdcIt != bdc.end(); bdcIt++) {
      _f += (*bdcIt)->energy();
    }
    std::cout << "Initial minimum energy = " << _f << std::endl;
 
    _fSaved = _f;
    double fBest = _f;
    double fWorst = _f;
    _xSaved = _x;
    vector<double > xBest = _x;
    vector<double > xWorst = _x;



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
    vector<double > SymmIndex = _printingStretches->printMaster(0);
    cout << "SymmIndex [nodes, eta, theta] = " << SymmIndex[0] << " " << SymmIndex[1] << " " << SymmIndex[2] << endl;



    // Begin simulated annealing
    unsigned accepted = 0, step = 1, pt = 0, j = 0;
    for(step = 1; step <= _nSteps; step++)
    {
      for(pt = 0; pt < _size; pt++)
      {
	  // perform _size metropolis steps at current temperature T
	  bool changed = changeState(model);
	  if( changed )
	  {
	    accepted++;
	    if( _fSaved < fBest )
	    {
	      fBest = _fSaved;
	      xBest = _xSaved;
	    }
	    if( _fSaved > fWorst )
	    {
	      fWorst = _fSaved;
	      xWorst = _xSaved;
	    }
	  }

	  if (_ferroMagnetic)
	  {
	    Interactions = _printingStretches->FerroMagneticInteractions();
	    ofsInteractions << Interactions[0] << " " << Interactions[1] << endl;
	  }

	  if (_print)
	  {
	    // Solve model before printing to reset displacements to correct values
	    // If changed == false the current displacements are not correct before model is solved again
	    _solver->solve(model); // In theory we should reset also the nodal positions to their original ones 
	                           // but in practice it should not matter since the equilibrium configuration should be fairly close to the current one
	    // Compute average radius
	    double Ravg = 0.0;
	    unsigned int NodesCount = 0;
	    Body::NodeContainer Nodes = bdc[0]->nodes();
	    for(j  = 0; j < Nodes.size(); j++)
	    {
	      DeformationNode<3> * defNodes = dynamic_cast<DeformationNode<3> *>(Nodes[j]);
	      if (defNodes != NULL)
	      {
		Ravg += tvmet::norm2(defNodes->point());
		NodesCount++;
	      }
	    }
	    Ravg /= NodesCount;

	    // Print energy values on file for every sub-step
	    ofsE << _f << " " << _fSaved << " " <<  fBest << " " << Ravg << std::endl;
	    // Print variable values on file for every sub-step
	    for (j = 0; j < _x.size(); j++)
	      ofsV << _x[j] << " ";
	    ofsV << std::endl;
	  }
      } // End of sub-loop at constant T


      // Solve model before printing to reset displacements to correct values
      // Reset mechanics dof in initial position
      Body::NodeContainer Nodes = bdc[0]->nodes();
      for(j = 0; j < Nodes.size(); j++)
      {
	DeformationNode<3> * defNodes = dynamic_cast<DeformationNode<3> *>(Nodes[j]);
	if (defNodes != NULL) {
	  defNodes->setPoint(defNodes->position() );
	}
      }
      _solver->solve(model);
      SymmIndex = _printingStretches->printMaster(step);
      cout << "SymmIndex [nodes, eta, theta] = " << SymmIndex[0] << " " << SymmIndex[1] << " " << SymmIndex[2] << endl;
     


      // Compute average radius
      double Ravg = 0.0;
      unsigned int NodesCount = 0;
      for(j = 0; j < Nodes.size(); j++)
      {
	DeformationNode<3> * defNodes = dynamic_cast<DeformationNode<3> *>(Nodes[j]);
	if (defNodes != NULL) 
	{
	  Ravg += tvmet::norm2(defNodes->point());
	  NodesCount++;
	}
      }
      Ravg /= NodesCount;
	
      // Print values of interest (energy, acceptance, Ravg)
      std::cout << "MTS iteration = " << step     << std::endl 
		<< "accepted      = " << accepted << std::endl 
		<< "_f            = " << _f       << std::endl 
		<< "_fSaved       = " << _fSaved  << std::endl 
		<< "fBest         = " << fBest    << std::endl
		<< "fWorst        = " << fWorst    << std::endl
		<< "Ravg          = " << Ravg     << std::endl;
     

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

 
    if (_print)
    {
      ofsE.close();
      ofsV.close();
    }

    if (_ferroMagnetic == 1)
    {
      ofsInteractions.close();
    }

    std::cout << "All done :) " << std::endl;
  }



// --------------------------------------------------------------
// Compute a new random trial state based on chosen distribution
// --------------------------------------------------------------
  bool MontecarloTwoStages::changeState(Model *model)
  {
    unsigned int i = 0, j = 0, ind = 0;

    double angle = M_PI/3.0;
    double Theta[2] = {-angle, angle};
    double Eta[3] = {-0.01, 0.01}; 
    double tol = 1.0e-8, temp = 0.0;
      
      i = rand()%_montecarloDoF.size();
      for (j = 0; j < i; j++) {
	ind +=  _montecarloDoF[i].size();
      }

      switch (_varType[i])
	{
	case -1: // Reset to zero

	  for (j = 0; j < _montecarloDoF[i].size(); j++)
	  {
	    _montecarloDoF[i][j]->setPoint(0.0);
	  }

	  break;

	case 0: // Change all sequentially
	  _x[ind+_flip] += Theta[rand()%2];
	  _montecarloDoF[i][_flip]->setPoint(_x[ind+_flip]);
	  _flip++;

	  if (_flip ==  _montecarloDoF[i].size() ) 
	  {
	    _flip = 0;
	  }
	  ind += _montecarloDoF[i].size();

	  break; 

	case 1: // Change one randomly - the option we like best for now (12/12/2012, LEP)
	  _flip = rand()%(_montecarloDoF[i].size());
	  _x[ind+_flip] += Theta[rand()%2];
	
	  _montecarloDoF[i][_flip]->setPoint(_x[ind+_flip]);

	  // ind += _montecarloDoF[i].size();
		
	  break;

	case 10:
	  _x[ind+_flip] += Eta[rand()%2];
	  _montecarloDoF[i][_flip]->setPoint(_x[ind+_flip]);
	  _flip++;

	  if (_flip ==  _montecarloDoF[i].size() ) 
	  {
	    _flip = 0;
	  }
	  ind += _montecarloDoF[i].size();

	  break;

	case 11: // Change one randomly - the option we like best for now (12/12/2012, LEP)
	  _flip = rand()%(_montecarloDoF[i].size());
	  temp =  Eta[rand()%2];
	  if ( (_x[ind+_flip] + temp) > (1.0-tol) ) {
	        _x[ind+_flip] += temp;
	  }
	
	  _montecarloDoF[i][_flip]->setPoint(_x[ind+_flip]);

	  // ind += _montecarloDoF[i].size();
		
	  break;
      
	default:
	  cout << "Option not found in MontecarloTwoStages::changeState - exiting " << endl;
	  exit(0);
	} // end of switch
    


    Model::BodyContainer bdc = model->bodies();
    // Reset mechanics dof in initial position
    Body::NodeContainer Nodes = bdc[0]->nodes();
    for(j = 0; j < Nodes.size(); j++)
    {
      DeformationNode<3> * defNodes = dynamic_cast<DeformationNode<3> *>(Nodes[j]);
      if (defNodes != NULL) {
	 defNodes->setPoint(defNodes->position() );
      }
    }
 


    // Compute Model energy after solving for displacement dof
    _solver->solve(model);
    
    // Compute model energy
    _f = 0.0;
    for(Model::ConstBodyIterator bdcIt = bdc.begin(); bdcIt != bdc.end(); bdcIt++)
    {
      (*bdcIt)->compute( true, false, false );
      _f += (*bdcIt)->energy();
    }
    
    // Decide whether or not to keep the new state
    double df = _f - _fSaved; 
       
    // metropolis
    double p = ((double)rand())/RAND_MAX;
    cout << "df = " << df << endl;
    if( df < 0.0 || p < exp( -df/_T1 ) )
    { 
      _xSaved = _x;
      _fSaved = _f;
      return true;
    }
      
    _x = _xSaved;
    ind = 0;
    for (i = 0; i< _montecarloDoF.size(); i++)
    {
      for (j = 0; j < _montecarloDoF[i].size(); j++)
      {
	_montecarloDoF[i][j]->setPoint(_x[ind]);
	ind++;
      }
    }

    return false;  
  }

}  // namespace voom

