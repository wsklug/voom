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

    _f = 0.0;
    Model::BodyContainer bdc = model->bodies();
    for(Model::ConstBodyIterator bdcIt = bdc.begin(); bdcIt != bdc.end(); bdcIt++) {
      _f += (*bdcIt)->energy();
    }
    _fSaved = _f;
    unsigned int i = 0;
    double ratio = 1.0;
    while (i < _maxIter && ratio > _finalRatio)
    {
      // Choose and assign new trial values
      changeState();
    
      // Minimize energy with respect to the rest of dof
      _solver->solve(model);
    
      // Compute model energy: keep new parameters choice if energy lower than _fSaved
      _f = 0.0;
      Model::BodyContainer bdc = model->bodies();
      for(Model::ConstBodyIterator bdcIt = bdc.begin(); bdcIt != bdc.end(); bdcIt++) {
	_f += (*bdcIt)->energy();
      }
      
      i++;
      ratio = fabs((_f-_fSaved)/_fSaved);

      std::cout << "MTS iteration = " << i     << std::endl 
		<< "ratio         = " << ratio     << std::endl 
		<< "_f            = " << _f        << std::endl 
		<< "_fSaved       = " << _fSaved   << std::endl 
		<< " _x.size()    = " << _x.size() << std::endl;

      if (_f < _fSaved)
      {
	_fSaved = _f;
	_xSaved = _x;
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
	
	ofsE << _f << " " << Ravg << std::endl;
	for (j = 0; j < _x.size(); j++)
	  ofsV << _x[j] << " ";
	ofsV << std::endl;
      }

    }

    std::cout << "_f      = " << _f      << std::endl
	      << "_fSaved = " << _fSaved << std::endl;

    if (_print)
    {
      ofsE.close();
      ofsV.close();
    }
   
  }



// --------------------------------------------------------------
// Compute a new random trial state based on chosen distribution
// --------------------------------------------------------------
  void MontecarloTwoStages::changeState()
  {
    // Distribution type is now hard wired - need to be modified
    unsigned int i = 0, j = 0;
    double temp;
    srand( time(NULL) );
    
    for (i = 0; i < _shearDoF.size(); i++)
    {
      _x[i] = 0.3*(double)(rand())/RAND_MAX - 0.15;
      _shearDoF[i]->setPoint(_x[i]);
    }

    double Theta[3] = {0, 60, 120}; 
    for (i = 0; i < _directionDoF.size(); i++)
    {
       _x[i+ _shearDoF.size()] = Theta[rand()%3];
      _directionDoF[i]->setPoint(_x[i+ _shearDoF.size()]);
    }

  }

}

