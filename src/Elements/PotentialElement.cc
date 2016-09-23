// -*- C++ -*-
//----------------------------------------------------------------------
//
//                 William S. Klug, Luigi Perotti
//                University of California Los Angeles
//                 (C) 2004-2007 All Rights Reserved
//
//----------------------------------------------------------------------
//

#include <iostream>
#include "PotentialElement.h"

namespace voom {
        
  void PotentialElement::compute(bool fl0, bool fl1, bool fl2)
  {
    if (fl0) {
    _energy = 0.0;
    }
	    
    for (set<DeformationNode<3> *>::iterator pNode = _domain.begin();
	pNode != _domain.end(); pNode++)
    {
	_mat->updateState(_center, *pNode, fl0, fl1, fl2);
	if (fl0) {
	  _energy += _mat->energy();
	}
    }

  }

  void PotentialElement::getTensions(vector<Vector3D > & OneElementsMidPoints, vector<double > & OneElementTension) {

    for (set<DeformationNode<3> *>::iterator pNode = _domain.begin(); pNode != _domain.end(); pNode++)
    {
      OneElementTension.push_back(_mat->computeTension(_center, *pNode));
      Vector3D MidPoint;
      MidPoint = 0.5*( _center->point() + (*pNode)->point() );
      OneElementsMidPoints.push_back(MidPoint);
    }

  }



} // namespace voom
