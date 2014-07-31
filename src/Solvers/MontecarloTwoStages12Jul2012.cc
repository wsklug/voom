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
    ofstream ofsE, ofsV;
    if (_print)
    {	
      ofsE.open("MontecarloEnergy.dat");
      ofsV.open("MontecarloVariables.dat");
      if (!ofsE) { std::cout << "Cannot open output file MontecarloEnergy.dat" << std::endl;
	exit(0); }
      if (!ofsV) { std::cout << "Cannot open output file MontecarloVariables.dat" << std::endl;
	exit(0); }
    }

    // Solve the problem at the beginning to release pre-existing stresses
    // Then start Montecarlo on the relaxed structure
    _solver->solve(model);

    _f = 0.0;
    Model::BodyContainer bdc = model->bodies();
    for(Model::ConstBodyIterator bdcIt = bdc.begin(); bdcIt != bdc.end(); bdcIt++) {
      _f += (*bdcIt)->energy();
    }
 
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
    // srand(static_cast<unsigned>(time(0)));
    int RandSeed = int(time(NULL));
    // int RandSeed = 36826794;
    srand(RandSeed);
    cout << "Seed for random number generator is " << RandSeed << endl;
    


    // Print initial configuration
    _printingStretches->printMaster(0);


    // Begin simulated annealing
    unsigned accepted = 0, step = 1, pt = 0, printSeqence = 1;
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

	  if (_print)
	  {
	    // Compute average radius
	    double Ravg = 0.0;
	    unsigned int j = 0;
	    Body::NodeContainer Nodes = bdc[0]->nodes();
	    for(j = 0; j < Nodes.size(); j++)
	    {
	      DeformationNode<3> * defNodes = dynamic_cast<DeformationNode<3> *>(Nodes[j]);
	      if (defNodes != NULL)
		Ravg += tvmet::norm2(defNodes->point());
	    }
	    Ravg /= Nodes.size();
	    
	    ofsE << _f << " " << _fSaved << " " <<  fBest << " " << Ravg << std::endl;
	    for (j = 0; j < _x.size(); j++)
	      ofsV << _x[j] << " ";
	    ofsV << std::endl;
	  }

	  if (pt == int(ceil(double(_size)/2.0)))
	  {
	    _solver->solve(model);
	    _printingStretches->printMaster(printSeqence);
	    printSeqence++;
	  }
      }


      // Solve model before printing to reset displacements to correct values
      _solver->solve(model);
      _printingStretches->printMaster(printSeqence);
      printSeqence++;



      std::cout << "MTS iteration = " << step         << std::endl 
		<< "accepted     = " << accepted << std::endl 
		<< "_f            = " << _f        << std::endl 
		<< "_fSaved       = " << _fSaved   << std::endl 
		<< "fBest        = " << fBest    << std::endl;
      /*
      if (_print)
      {
	// Compute average radius
	double Ravg = 0.0;
	unsigned int j = 0;
	Body::NodeContainer Nodes = bdc[0]->nodes();
	for(j = 0; j < Nodes.size(); j++)
	{
	  DeformationNode<3> * defNodes = dynamic_cast<DeformationNode<3> *>(Nodes[j]);
	  if (defNodes != NULL)
	    Ravg += tvmet::norm2(defNodes->point());
	}
	Ravg /= Nodes.size();
	
	ofsE << _f << " " << _fSaved << " " <<  fBest << " " << Ravg << std::endl;
	for (j = 0; j < _x.size(); j++)
	  ofsV << _x[j] << " ";
	ofsV << std::endl;
      }
      */


      switch( _Tsched )
      {
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

    }

 
    if (_print)
    {
      ofsE.close();
      ofsV.close();
    }

    std::cout << "_accepted = " << accepted << std::endl
	      << "_f        = " << _f        << std::endl
	      << "_fSaved   = " << _fSaved   << std::endl
	      << "_fBest    = " << fBest    << std::endl;
  }



// --------------------------------------------------------------
// Compute a new random trial state based on chosen distribution
// --------------------------------------------------------------
  bool MontecarloTwoStages::changeState(Model *model)
  {
    unsigned int i = 0, j = 0, ind = 0;
    double angle = M_PI/3.0;
    // double Theta[3] = {0.0, angle, 2.0*angle}; 
    double Theta[2] = {-angle, angle};
    double Eta[3] = {-0.1, 0.0, 0.1}; 

    for (i = 0; i< _montecarloDoF.size(); i++)
    {
      // cout << _varType[i] << endl;
      switch (_varType[i])
      {
      case -1:
	   for (j = 0; j < _montecarloDoF[i].size(); j++)
	   {
	     _montecarloDoF[i][j]->setPoint(0.0);
	   }
	   break;
      case 0:
	   for (j = 0; j < _montecarloDoF[i].size(); j++)
	   {
	   //_x[ind]  = Theta[rand()%3];
	     _x[ind] += Theta[rand()%3];
	   // cout <<  _x[ind] << endl;
	     _montecarloDoF[i][j]->setPoint(_x[ind]);
	     /*
	     // generate dv using gausian distribution
	     double u1 = (double)(rand())/RAND_MAX;
	     double u2 = (double)(rand())/RAND_MAX;
	     double dv = _T2*AngleScale*sqrt(-2.0*log(u1))*cos(2.0*M_PI*u2);
	     if( finite(dv) ) _x[i+ _shearDoF.size()] += dv;

	     _directionDoF[i]->setPoint(_x[i+ _shearDoF.size()]);
	     */
	     ind++;
	   }
	   break;
      case 1:
	//_x[ind+_flip] = Theta[rand()%3];
	_x[ind+_flip] += Theta[rand()%2];
	
	_montecarloDoF[i][_flip]->setPoint(_x[ind+_flip]);
	_flip++;
	if (_flip ==  _montecarloDoF[i].size() ) {_flip = 0;}
	ind += _montecarloDoF[i].size();
		
	break;
      case 10:
	for (j = 0; j < _montecarloDoF[i].size(); j++)
	   {
	     _x[ind] = Eta[rand()%3];
	     _montecarloDoF[i][j]->setPoint(_x[ind]);
	     ind++;
	   }
	   break;
      case 20:
	for (j = 0; j < _montecarloDoF[i].size(); j++)
	   {
	     // generate dv using gausian distribution
	     double u1 = (double)(rand())/RAND_MAX;
	     double u2 = (double)(rand())/RAND_MAX;
	     double dv = _T2*sqrt(-2.0*log(u1))*cos(2.0*M_PI*u2);
	     if( finite(dv) ) _x[ind] += dv;
	     // cout <<  _x[ind] << endl;
	     _montecarloDoF[i][j]->setPoint(_x[ind]);
	    
	     ind++;
	   }
	   break;
      default:
	cout << "Option not found in MontecarloTwoStages::changeState - exiting " << endl;
	exit(0);
      }
    }  
    

  /*
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
  */


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
    // Minimize energy with respect to the rest of dof
    _solver->solve(model);
    
    // Compute model energy: keep new parameters choice if energy lower than _fSaved
    _f = 0.0;
    for(Model::ConstBodyIterator bdcIt = bdc.begin(); bdcIt != bdc.end(); bdcIt++) {
      (*bdcIt)->compute( true, false, false );
      _f += (*bdcIt)->energy();
    }


    double _ftest = 0.0;
    Body::ElementContainer Elements = bdc[0]->elements();
    for(int el = 0; el < Elements.size(); el++) {
      _ftest += Elements[el]->energy();
    }
    cout << "_f = " << _f << "  _ftest = " << _ftest << endl;











    
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
      
      // _f = _fSaved;

      return false;  
  }

}  // namespace voom

