// -*- C++ -*-
//----------------------------------------------------------------------
//
//                   William S. Klug, Luigi Perotti
//                University of California Los Angeles
//                 (C) 2004-2007 All Rights Reserved
//
//----------------------------------------------------------------------
//
/*!
  \file PotentialElement.h
  \brief Interface for a generic potential element
*/

#ifndef _POTENTIALELEMENT_H_
#define _POTENTIALELEMENT_H_

#include "Element.h"
#include "Potential.h"
#include "VoomMath.h"
#include "Node.h"

#include<set>

using namespace std;

namespace voom {

	class PotentialElement : public Element
	{

	public:
		//! Construct potential from necessary data
		PotentialElement(Potential *mat, DeformationNode<3> * center,
			set<DeformationNode<3> * > domain) : _mat(mat), _center(center), _domain(domain)
		{
			for (set<DeformationNode<3> *>::iterator pNode = _domain.begin();
				pNode != _domain.end(); pNode++)
			{
				_baseNodes.push_back(*pNode);
			}
			_baseNodes.push_back(_center);
		};

		//! Destructor
		~PotentialElement() {};

		// Operators
		//! Based on nodal postitions, calculates energy and nodal forces
		void compute(bool f0, bool f1, bool f2);

		void resetDomain(set<DeformationNode<3> * > NewDomain) { _domain = NewDomain; };

		void getTensions(vector<Vector3D > & OneElementsMidPoints, vector<double > & OneElementTension);

	private:
		//! No access to default constructor
		PotentialElement() {};

		// Members:
		Potential * _mat;
		DeformationNode<3> * _center;
		set<DeformationNode<3> * > _domain;

	}; // Potential element class

}
#endif // _POTENTIALELEMENT_H_
