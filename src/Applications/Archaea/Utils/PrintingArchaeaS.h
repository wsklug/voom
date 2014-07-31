
//----------------------------------------------------------------------
//                   Luigi Perotti, William S. Klug
//                University of California Los Angeles
//                   (C) 2004 All Rights Reserved
//----------------------------------------------------------------------

#if !defined(__PrintingArchaeaS_h__)
#define __PrintingArchaeaS_h__

#include <string>
#include <iostream>
#include <vector>
#include <fstream>

#include "EvansElastic.h"
#include "SCElastic.h"
#include "TriangleQuadrature.h"
#include "ShapeTri3.h"
#include "LoopShellBody.h"
#include "Potential.h"

using namespace std;

namespace voom
{
  class PrintingArchaeaS
  {
  private:	

    const string _meshFileName;
    const string _defMeshFileName;
    const Model::BodyContainer & _bdc;
    const double _RconnSearch;
    const vector<DeformationNode<3>* > & _proteinPos;
    vector<uint > _dislocations;
    vector<Vector3D > _previousNodalPositions;
    Potential * _mat;


  public:
    
    // Constructor
  PrintingArchaeaS(const string MeshFileName,  
		   const string DefMeshFileName,
		   const Model::BodyContainer & BDC,
		   const double RconnSearch,
		   const vector<DeformationNode<3>* > & ProteinPos,
		   Potential * mat): 
    _meshFileName(MeshFileName), _defMeshFileName(DefMeshFileName),
    _bdc(BDC),
    _RconnSearch(RconnSearch), _proteinPos(ProteinPos),
    _mat(mat)  
    {
      Body::NodeContainer Nodes = _bdc[0]->nodes();
      _previousNodalPositions.resize(Nodes.size(), Vector3D(0.0));
      
      for (uint i = 0; i < Nodes.size(); i++)
      {
	DeformationNode<3> * node = dynamic_cast<DeformationNode<3>* >(Nodes[i]);
	if (node != NULL) 
	{
	  _previousNodalPositions[i] = node->point();
	}
      } // loop over all nodes
    };
    


    // Destructor
    virtual ~ PrintingArchaeaS () {}

    // Printing functions
    vector<uint> printMaster(const int iter)
    {
	stringstream DefMeshFileStream;
	DefMeshFileStream << _defMeshFileName << iter << ".vtk";
	printLattice(DefMeshFileStream.str());
	stringstream InteractionsFileStream;
	InteractionsFileStream << _defMeshFileName << "_Interactions_" << iter << ".vtk";
	printInteractions(InteractionsFileStream.str());
	return _dislocations;
    }



    void printLattice(const string DefMeshFileName)
    {
      ofstream ofs1(DefMeshFileName.c_str());
      if (!ofs1) {
	std::cout << "Error: cannot open output ("
		  << DefMeshFileName
		  << ") file." << std::endl;
	exit(0);
      }

      unsigned int ind = 0;
      // Shell body
      Body::ElementContainer Elements = _bdc[0]->elements();
      Body::NodeContainer Nodes = _bdc[0]->nodes();
      // vector<DeformationNode<3> *> Nodes; // = _shellBody->shellsNodes();
      // vector<LoopShell<EvansElastic> *> Elements = _shellBody->shells();
      
      
      // Node Section
      ofs1 << "# vtk DataFile Version 3.0" << endl
	   << "Test example" << endl
	   << "ASCII" << endl
	   << "DATASET POLYDATA" << endl
	   << "POINTS  " << Nodes.size() << "  double" << endl;

      // Output nodal postions
      for (Body::NodeIterator pn = Nodes.begin(); pn!= Nodes.end(); pn ++)
      {
	DeformationNode<3> * node = dynamic_cast<DeformationNode<3>* >(*pn);
	if (node != NULL) 
	{
	  const Vector3D & nodalPos =  node->point();
	  ofs1 << std::setprecision(16) 
	       << nodalPos(0) << "  "
	       << nodalPos(1) << "  "
	       << nodalPos(2) << std::endl;
	}
      } // loop over all nodes



      // Create vector to hold node connectivity
      vector<uint > NodesConn(Nodes.size(), 0);



      // Element Section
      int nel = Elements.size();
      ofs1 << "POLYGONS  " << nel << "  "
	   << 4*nel << endl;
      for(int e = 0; e < nel; e++)
      {
	if(!(_bdc[0]->active(e)) ) exit(0); // All elements should be active in this application
	uint a = Elements[e]->baseNodes()[0] -> id();
	uint b = Elements[e]->baseNodes()[1] -> id();
	uint c = Elements[e]->baseNodes()[2] -> id();
	ofs1 << 3 << "  "
	     << setw(10) << a
	     << setw(10) << b
	     << setw(10) << c
	     << endl;
	NodesConn[a] += 1;
	NodesConn[b] += 1;
	NodesConn[c] += 1;
      } // All elements are printed in the deformed mesh file



      // Cell data section
      ofs1 << "CELL_DATA    " << nel << endl;
  
      // Element elastic energy
      ofs1 << "SCALARS    ElEnDensity    double    1" << endl;
      ofs1 << "LOOKUP_TABLE default" << endl;
      for(int e = 0; e < nel; e++) {
	ofs1 << (Elements[e]->energy()) << endl;
      }
      ofs1 << std::endl;


      
      // Node data section
      vector<uint > Valence20(_proteinPos.size(), 0),  Valence50(_proteinPos.size(), 0),  Valence80(_proteinPos.size(), 0);

      for (uint i = 0; i < _proteinPos.size(); i++)
      {
	vector<double> NeighborDist;
	double MinDist = _RconnSearch;
	for (uint j = 0; j < _proteinPos.size(); j++)
	{
	  double R = tvmet::norm2(_proteinPos[i]->point() - _proteinPos[j]->point());
	  if ( i!=j && R < _RconnSearch )
	  {
	    NeighborDist.push_back(R);
	    if (MinDist > R) {
	      MinDist = R;
	    }
	  } 
	}
	
	for (uint j = 0; j < NeighborDist.size(); j++)
	{
	  if (NeighborDist[j] <= MinDist*1.2) {
	    Valence20[i] += 1;
	  }
	  if (NeighborDist[j] <= MinDist*1.5) {
	    Valence50[i] += 1;
	  }
	  if (NeighborDist[j] <= MinDist*1.8) {
	    Valence80[i] += 1;
	  }
	}
      } // i loop
	  
      ofs1 << "POINT_DATA    " << Nodes.size() << endl;
      ofs1 << "SCALARS    NodesConn    double    1" << endl;
      ofs1 << "LOOKUP_TABLE default" << endl;
      _dislocations.assign(3, 0);
      for (uint i = 0; i < Nodes.size(); i++)
      {
	ofs1 << NodesConn[i] << endl;
	if (NodesConn[i] == 5) {
	  _dislocations[0] += 1; }
	else if (NodesConn[i] == 6) {
	  _dislocations[1] += 1; }
	else if (NodesConn[i] == 7) {
	  _dislocations[2] += 1; }
      }

      ofs1 << "SCALARS    Valence20    double    1" << endl;
      ofs1 << "LOOKUP_TABLE default" << endl;
      for (uint i = 0; i < _proteinPos.size(); i++)
      {
	ofs1 << Valence20[i] << endl;
      }

      ofs1 << "SCALARS    Valence50    double    1" << endl;
      ofs1 << "LOOKUP_TABLE default" << endl;
      for (uint i = 0; i < _proteinPos.size(); i++)
      {
	ofs1 << Valence50[i] << endl;
      }
      
      ofs1 << "SCALARS    Valence80    double    1" << endl;
      ofs1 << "LOOKUP_TABLE default" << endl;
      for (uint i = 0; i < _proteinPos.size(); i++)
      {
	ofs1 << Valence80[i] << endl;
      }
      
      // ofs1 << "POINT_DATA    " << Nodes.size() << endl;
      ofs1 << "VECTORS    IncU     double" << endl;
      for (uint i = 0; i < Nodes.size(); i++)
      {
	DeformationNode<3> * node = dynamic_cast<DeformationNode<3>* >(Nodes[i]);
	if (node != NULL) 
	{
	  const Vector3D & nodalPos =  node->point();
	  ofs1 << std::setprecision(16) 
	       << nodalPos(0) - _previousNodalPositions[i](0) << "  "
	       << nodalPos(1) - _previousNodalPositions[i](1) << "  "
	       << nodalPos(2) - _previousNodalPositions[i](2) << std::endl;
	  _previousNodalPositions[i] = nodalPos; // Update previous position for next time step
	}
      }
 
      // This is all for the def mesh file
      ofs1.close();
    };



    void printInteractions(const string InteractionsFileName)
    {
      ofstream ofs1(InteractionsFileName.c_str());
      if (!ofs1) {
	std::cout << "Error: cannot open output ("
		  << InteractionsFileName
		  << ") file." << std::endl;
	exit(0);
      }

      Body::ElementContainer Elements = _bdc[0]->elements();
      Body::NodeContainer Nodes = _bdc[0]->nodes();

      vector<Vector3D > MidNodesPositions;
      vector<double > tension;
      for(int e = 0; e < Elements.size(); e++)
      {
	uint a = Elements[e]->baseNodes()[0] -> id();
	uint b = Elements[e]->baseNodes()[1] -> id();
	uint c = Elements[e]->baseNodes()[2] -> id();
	
	DeformationNode<3> * nodeA = dynamic_cast<DeformationNode<3>* >(Nodes[a]);
	DeformationNode<3> * nodeB = dynamic_cast<DeformationNode<3>* >(Nodes[b]);
	DeformationNode<3> * nodeC = dynamic_cast<DeformationNode<3>* >(Nodes[c]);

	Vector3D MidPoint;

	MidPoint = 0.5*(nodeA->point() + nodeB->point());
	MidNodesPositions.push_back(MidPoint);	
	tension.push_back(_mat->computeTension(nodeA, nodeB));

	MidPoint = 0.5*(nodeA->point() + nodeC->point());
	MidNodesPositions.push_back(MidPoint);	
	tension.push_back(_mat->computeTension(nodeA, nodeC));

	MidPoint = 0.5*(nodeB->point() + nodeC->point());
	MidNodesPositions.push_back(MidPoint);	
	tension.push_back(_mat->computeTension(nodeB, nodeC));
	
      } 

      /*
	vector<PotentialElement* > PotentialElements = PotBody->getPotentialElements();
 
	vector<Vector3D > MidNodesPositions;
	vector<double > tension;
	for (uint i = 0; i < PotentialElements.size(); i++)
	{
	vector<Vector3D > OneElementMidPoints;
	vector<double > OneElementTension;
	PotentialElements[i]->getTensions(OneElementMidPoints, OneElementTension);
	MidNodesPositions.insert(MidNodesPositions.end(), OneElementMidPoints.begin(), OneElementMidPoints.end());
	tension.insert(tension.end(), OneElementTension.begin(), OneElementTension.end());
	}
      */

      // Node Section
      ofs1 << "# vtk DataFile Version 3.0" << endl
	   << "Test example" << endl
	   << "ASCII" << endl
	   << "DATASET POLYDATA" << endl
	   << "POINTS  " << MidNodesPositions.size() << "  double" << endl;

      // Output nodal postions
      for (uint i = 0; i < MidNodesPositions.size(); i++)
      {
	ofs1 << std::setprecision(16) 
	     << MidNodesPositions[i](0) << "  "
	     << MidNodesPositions[i](1) << "  "
	     << MidNodesPositions[i](2) << std::endl;
      }
     
      ofs1 << "POINT_DATA    " << tension.size() << endl;
      ofs1 << "SCALARS    tension    double    1" << endl;
      ofs1 << "LOOKUP_TABLE default" << endl;
      for (uint i = 0; i < tension.size(); i++)
      {
	ofs1 << tension[i] << endl;
      }

    } // print Interactions

  }; // PrintingArchaeaS
  
}; // namespace voom

#endif // __PrintingArchaeaS_h__


