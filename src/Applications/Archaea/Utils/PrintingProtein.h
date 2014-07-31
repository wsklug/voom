
//----------------------------------------------------------------------
//                   Luigi Perotti, William S. Klug
//                University of California Los Angeles
//                   (C) 2004 All Rights Reserved
//----------------------------------------------------------------------

#if !defined(__PrintingProtein_h__)
#define __PrintingProtein_h__

#include <string>
#include <iostream>
#include <vector>
#include <fstream>

#include "ProteinBody.h"

using namespace std;

namespace voom
{
  class PrintingProtein
  {
  private:	
    
    const string _meshFileName;
    const string _defMeshFileName;
    const vector<NodeBase* > & _nodes;
    const vector< tvmet::Vector<int,3> > & _connectivities;

    const vector<ProteinNode* > & _proteins;
    double _RconnSearch;
    
    
  public:
    
    // Constructor
    PrintingProtein(const string MeshFileName,  
		    const string DefMeshFileName,
		    const vector<NodeBase* > & Nodes,
		    const vector< tvmet::Vector<int,3> > & Connectivities,
		    const vector<ProteinNode* > & Proteins,
		    double RconnSearch): 
    _meshFileName(MeshFileName), _defMeshFileName(DefMeshFileName),
      _nodes(Nodes), _connectivities(Connectivities),
      _proteins(Proteins), _RconnSearch(RconnSearch) {};
    
    // Destructor
    virtual ~PrintingProtein() {};
    
    // Printing functions
    void printMaster(const int iter, int printAll = 1)
    {
      if (printAll == 1) {
	stringstream DefMeshFileStream;
	DefMeshFileStream << _defMeshFileName << iter << ".vtk";
	printLattice(DefMeshFileStream.str());
      }
      stringstream ParticleFileStream;
      ParticleFileStream << _defMeshFileName << "_Pr_" << iter << ".vtk";
      printParticle(ParticleFileStream.str());
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
      
      // Node Section
      ofs1 << "# vtk DataFile Version 3.0" << endl
	   << "Test example" << endl
	   << "ASCII" << endl
	   << "DATASET POLYDATA" << endl
	   << "POINTS  " << _nodes.size() << "  double" << endl;

      // Output nodal postions
      for (uint i = 0; i < _nodes.size(); i++)
      {
	DeformationNode<3> * Node = dynamic_cast<DeformationNode<3>* >(_nodes[i]);
	if (Node != NULL) 
	{
	  const Vector3D & NodalPos =  Node->point();
	  ofs1 << std::setprecision(16) 
	       << NodalPos(0) << "  "
	       << NodalPos(1) << "  "
	       << NodalPos(2) << std::endl;
	}
      } // loop over all nodes

      int nel = _connectivities.size();
      ofs1 << "POLYGONS  " << nel << "  "
	   << 4*nel << endl;
      for(int e = 0; e < nel; e++)
      {
	ofs1 << 3 << "  "
	     << setw(10) << _connectivities[e](0)
	     << setw(10) << _connectivities[e](1)
	     << setw(10) << _connectivities[e](2)
	     << endl;
      } // All elements are printed in the deformed mesh file


      // This is all for the def mesh file
      ofs1.close();
    };



    void printParticle(const string ParticleFileName)
    {
      ofstream ofs1(ParticleFileName.c_str());
      if (!ofs1) {
	std::cout << "Error: cannot open output ("
		  << ParticleFileName
		  << ") file." << std::endl;
	exit(0);
      }

      uint NumProteins = _proteins.size();
      vector<DeformationNode<3>::Point > ProteinPositions(NumProteins, Vector3D(0.0) );
      for (uint i = 0; i < NumProteins; i++)
      {
	ProteinPositions[i] = (_proteins[i]->getHost())->point();
      }

      // Node data section
      vector<uint > Valence20(NumProteins, 0),  Valence50(NumProteins, 0),  Valence80(NumProteins, 0);

      for (uint i = 0; i < NumProteins; i++)
      {
	vector<double> NeighborDist;
	double MinDist = _RconnSearch;
	for (uint j = 0; j < NumProteins; j++)
	{
	  double R = tvmet::norm2(ProteinPositions[i] - ProteinPositions[j]);
	  if ( i != j && R < _RconnSearch )
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
	  

      // Node Section
      ofs1 << "# vtk DataFile Version 3.0" << endl
	   << "Test example" << endl
	   << "ASCII" << endl
	   << "DATASET POLYDATA" << endl
	   << "POINTS  " << NumProteins << "  double" << endl;

      // Output nodal postions
      for (uint i = 0; i < NumProteins; i++)
      {
	DeformationNode<3>::Point A = ProteinPositions[i];
	ofs1 << std::setprecision(16) 
	     << A(0) << "  " << A(1) << " " << A(2) << std::endl;
      }
     
      ofs1 << "POINT_DATA    " << NumProteins << endl;
      ofs1 << "SCALARS    Valence20    double    1" << endl;
      ofs1 << "LOOKUP_TABLE default" << endl;
      for (uint i = 0; i < NumProteins; i++)
      {
	ofs1 << Valence20[i] << endl;
      }

      ofs1 << "SCALARS    Valence50    double    1" << endl;
      ofs1 << "LOOKUP_TABLE default" << endl;
      for (uint i = 0; i < NumProteins; i++)
      {
	ofs1 << Valence50[i] << endl;
      }

      ofs1 << "SCALARS    Valence80    double    1" << endl;
      ofs1 << "LOOKUP_TABLE default" << endl;
      for (uint i = 0; i < NumProteins; i++)
      {
	ofs1 << Valence80[i] << endl;
      }

    } // print Interactions

  }; // PrintingArchaeaS
  
}; // namespace voom

#endif // __PrintingArchaeaS_h__


