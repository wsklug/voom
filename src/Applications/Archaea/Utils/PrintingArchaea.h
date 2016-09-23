
//----------------------------------------------------------------------
//                   Luigi Perotti, William S. Klug
//                University of California Los Angeles
//                   (C) 2004 All Rights Reserved
//----------------------------------------------------------------------

#if !defined(__PrintingArchaea_h__)
#define __PrintingArchaea_h__

#include <string>
#include <iostream>
#include <vector>
#include <fstream>

#include "EvansElastic.h"
#include "SCElastic.h"
#include "TriangleQuadrature.h"
#include "ShapeTri3.h"

using namespace std;

namespace voom
{
  class PrintingArchaea
  {
  private:	

    const string _meshFileName;
    const string _outlineFileName;  
    const unsigned int _nPout; 
    
    const string _defMeshFileName;
    const string _defOutlineFileName;
    
    const Model::BodyContainer & _bdc;
    const int _refinement;


  public:
    
    // Constructor
  PrintingArchaea(const string MeshFileName, const string OutlineFileName,  const unsigned int NPout, 
		  const string DefMeshFileName, const string DefOutlineFileName, 
		  const Model::BodyContainer & BDC,
		  const int refinement): 
    _meshFileName(MeshFileName), _outlineFileName(OutlineFileName), _nPout(NPout),
      _defMeshFileName(DefMeshFileName), _defOutlineFileName(DefOutlineFileName),
      _bdc(BDC), _refinement(refinement) {};
    
    // Destructor
    virtual ~ PrintingArchaea () {}

    // Printing functions
    void printMaster(const int iter)
    {
	stringstream DefMeshFileStream, DefOutlineFileStream;
	DefMeshFileStream << _defMeshFileName << iter << ".vtk";
	DefOutlineFileStream << _defOutlineFileName << iter << ".vtk";
	printLattice(DefMeshFileStream.str(), DefOutlineFileStream.str(), iter);
    }

    void printLattice(const string DefMeshFileName,
		      const string DefOutlineFileName,
		      const unsigned int iter)
    {
      ofstream ofs1(DefMeshFileName.c_str());
      if (!ofs1) {
	std::cout << "Error: cannot open output ("
		  << DefMeshFileName
		  << ") file." << std::endl;
	exit(0);
      }

      ofstream ofs2(DefOutlineFileName.c_str());
      if (!ofs2) {
	std::cout << "Error: cannot open output ("
		  << DefOutlineFileName
		  << ") file." << std::endl;
	exit(0);
      }

      unsigned int ind = 0;
      // Shell body
      Body::ElementContainer Elements = _bdc[0]->elements();
      Body::NodeContainer Nodes = _bdc[0]->nodes();
      
      // Node Section
      ofs1 << "# vtk DataFile Version 3.0" << endl
	   << "Test example" << endl
	   << "ASCII" << endl
	   << "DATASET POLYDATA" << endl
	   << "POINTS  " << Nodes.size() << "  double" << endl;

      ofs2 << "# vtk DataFile Version 3.0" << endl
	   << "vtk output" << endl
	   << "ASCII" << endl
	   << "DATASET POLYDATA" << endl
	   << "POINTS " << _nPout << " float" << endl;

      // Output nodal postions
      Body::NodeIterator pn; 
      for (pn = Nodes.begin(); pn!= Nodes.end(); pn ++)
      {
	DeformationNode<3> * node = dynamic_cast<DeformationNode<3>* >(*pn);

	if (node != NULL) 
	{
	  const Vector3D & nodalPos =  node->point();
	  ofs1 << std::setprecision(16) 
	       << nodalPos(0) << "  "
	       << nodalPos(1) << "  "
	       << nodalPos(2) << std::endl;
	  
	  // Update nodal position for the outline 
	  // and considering that the nodes are always in the same order (all the mesh files need to be consistent)
	  if (ind < _nPout)
	  {
	    ofs2 << std::setprecision(16) 
		 << nodalPos(0) << "  "
		 << nodalPos(1) << "  "
		 << nodalPos(2) << std::endl;
	    ind++;
	  } // outline nodes 
	} 
      } // loop over all nodes



      // Element Section
      int nel = Elements.size();
      ofs1 << "POLYGONS  " << nel << "  "
	   << 4*nel << endl;
      for(int e = 0; e < nel; e++)
      {
	if(!(_bdc[0]->active(e)) ) exit(0); // All elements should be active in this application
	ofs1 << 3 << "  "
	     << setw(10) << Elements[e]->baseNodes()[0] -> id()
	     << setw(10) << Elements[e]->baseNodes()[1] -> id()
	     << setw(10) << Elements[e]->baseNodes()[2] -> id()
	     << endl;
      } // All elements are printed in the deformed mesh file

      // Polygons in the outline file
      std::string OutlineFile = _outlineFileName+".vtk";
      // Outline for pentamers and examers
      ifstream ifs(OutlineFile.c_str());
      if (!ifs) {
	cout << "Error: Cannot open input file: " << OutlineFile << endl;
	exit(0);
      }

      string token;
      ifs >> token; 
      while( token != "LINES") ifs >> token; // Set cursor to reading position
      // Print header
      ofs2 << token << " ";
      ifs >> token;
      ofs2 << token << " ";
      ifs >> token;
      ofs2 << token << endl;
      ifs >> token;
      // Print polygons
      while(!ifs.eof())
	{
	  ofs2 << token << " ";
	  ifs >> token;
	}
      ifs.close();

      // This is all for the outline file 
      ofs2.close();

      

      



      // Cell section
      ofs1 << "CELL_DATA    " << nel << endl;
  
      // Element elastic energy
      ofs1 << "SCALARS    ElasticEnergy    double    1" << endl;
      ofs1 << "LOOKUP_TABLE default" << endl;
      for(int e = 0; e < nel; e++) {
	ofs1 << Elements[e]->energy() << endl;
      }
      ofs1 << std::endl;

      // This is all for the def mesh file
      ofs1.close();
    };
    
  }; // PrintingArchaea
  
}; // namespace voom

#endif // __PrintingArchaea_h__


