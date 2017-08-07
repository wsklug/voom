/*
 * OPSBody.h
 *
 *  Created on: Aug 5, 2017
 *      Author: amit
 */

#ifndef _OPSBODY_H_
#define _OPSBODY_H_

#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include "Body.h"
#include "voom.h"
#include "Node.h"
#include "VoomMath.h"

#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkKdTree.h>

#if VTK_MAJOR_VERSION < 6
#define SetInputData SetInput
#endif

using namespace std;

namespace voom {
/*!
 OPSBody is a concrete class derived from Body, implementing
 a oriented-particle system.
 Reference:
 	 Szeliski, Richard and Tonnesen, D. (1992). Surface Modeling with
 	 Oriented Particle Systems.	Siggraph ’92, 26(2), 160.
 	 https://doi.org/10.1017/CBO9781107415324.004
 */

class OPSBody: public Body {
public:

	typedef vector<OPSNode*>::const_iterator opsNodeIterator;

	struct properties{
		//Default values
		properties(){
			alphaM = 1.0; alphaP = 1.0; alphaN = 1.0; alphaC = 1.0;
			epsilon = 1.0; r_e = 1.0; s = 6.9314718056;//Fracture strain = 10%
			K = 1.0; a = 1.0; b = 1.0;
		}
		// Initialize from a vector of doubles
		properties( vector<double> val){
			assert( val.size() > 10 );
			alphaM = val[0]; alphaP = val[1]; alphaN = val[2]; alphaC = val[3];
			epsilon = val[4]; r_e = val[5]; s = val[6];
			K = val[7]; a = val[8]; b = val[9];
		}
		// Weights for the potentials involved in total energy
		double alphaM, alphaP, alphaN, alphaC;
		/* Morse potential parameters:
		 * epsilon = equilibrium energy, r_e = equilibrium separation
		 * s = Morse potential width-controlling parameter
		 */
		double epsilon, r_e, s;
		// Kernel parameters
		double K, a, b;
	};
	//! Construct body from
	OPSBody(const vector<OPSNode*> & nodes, properties p, double r);

	//! Destructor
	~OPSBody() {}

	//!Construct vtkPolyData from nodes
	void updatePolyDataAndKdTree();

	//! Do mechanics on Body
	void compute(bool f0, bool f1, bool f2);

	void updateNeighbors();

	//! Return the energy of the body
	double energy() const { return _energy;}

	//!Return the current OPS properties
	properties getProperties(){return _prop;}

	//!Update the current OPS properties
	void updateProperties(properties p){_prop = p;}

	//!Return the current search radius
	double getSearchRadius(){return _searchR;}

	//!Update the search radius
	void updateSearchRadius( double r ){ _searchR = r;}

	//!Get the polydata associated with the OPS
	vtkSmartPointer<vtkPolyData> getPolyData(){ return _polyData;}

	//!Get the K-d tree associated with the OPS
	vtkSmartPointer<vtkKdTree> getKdTree(){ return _kdTree;}

	//! General printing of a Paraview file
	void printParaview(const string name) const;

	//! Printing of a Paraview file after applying Kabsch algorithm
	void printParaview(const string name, Eigen::Matrix3Xd) const;

	void pushBack(Element* e) { _elements.push_back(e);}

	//! Return the mean squared displacement
	double rmsd();

	//!Get nearest neighbor for rmsd calculation
	std::vector<int> getInitialNearestNeighbor(){
		return _initialNearestNeighbor;
	}

	//!OPS kernel and potential functions
	double morse(Vector3D xi, Vector3D xj);
	double psi(Vector3D vi, Vector3D xi, Vector3D xj);
	double phi_p(Vector3D vi, Vector3D xi, Vector3D xj);
	double phi_n(Vector3D vi, Vector3D vj);
	double phi_c(Vector3D vi, Vector3D vj, Vector3D xi, Vector3D xj);

	//!OPS kernel and potential derivatives
	Vector3D DmorseDxi(Vector3D xi, Vector3D xj);
	Vector3D DmorseDxj(Vector3D xi, Vector3D xj);

	Vector3D DpsiDvi(Vector3D vi, Vector3D xi, Vector3D xj);
	Vector3D DpsiDxi(Vector3D vi, Vector3D xi, Vector3D xj);
	Vector3D DpsiDxj(Vector3D vi, Vector3D xi, Vector3D xj);

	Vector3D Dphi_pDvi(Vector3D vi, Vector3D xi, Vector3D xj);
	Vector3D Dphi_pDxi(Vector3D vi, Vector3D xi, Vector3D xj);
	Vector3D Dphi_pDxj(Vector3D vi, Vector3D xi, Vector3D xj);

	Vector3D Dphi_nDvi(Vector3D vi, Vector3D vj);
	Vector3D Dphi_nDvj(Vector3D vi, Vector3D vj);

	Vector3D Dphi_cDvi(Vector3D vi, Vector3D vj, Vector3D xi, Vector3D xj);
	Vector3D Dphi_cDvj(Vector3D vi, Vector3D vj, Vector3D xi, Vector3D xj);
	Vector3D Dphi_cDxi(Vector3D vi, Vector3D vj, Vector3D xi, Vector3D xj);
	Vector3D Dphi_cDxj(Vector3D vi, Vector3D vj, Vector3D xi, Vector3D xj);

protected:
	const vector<OPSNode*> & _opsNodes; // Nodes
	double _searchR; // Search radius
	properties _prop; // Parameters for the OPS
	vtkSmartPointer<vtkKdTree> _kdTree;
	vtkSmartPointer<vtkPolyData> _polyData;
	vector<vtkSmartPointer<vtkIdList> > _neighbors;
	std::vector<int> _initialNearestNeighbor; // Needed to find rmsd

};// OPS body

}
#endif /* _OPSBODY_H_ */
