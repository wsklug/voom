
//----------------------------------------------------------------------
//                   Luigi Perotti, William S. Klug
//                University of California Los Angeles
//                   (C) 2004 All Rights Reserved
//----------------------------------------------------------------------

#if !defined(__PrintingStretches_h__)
#define __PrintingStretches_h__

#include "voom.h"
#include <string>
#include <iostream>
#include <vector>
#include <fstream>


#include "C0MembraneShear.h"
#include "C0MembraneStretch.h"
#include "C0Membrane.h"
#include "EvansElastic.h"
#include "SCElastic.h"
// #include "EvansElastic_SkewedMin.h"
#include "TriangleQuadrature.h"
#include "ShapeTri3.h"

using namespace std;

namespace voom
{
 class PrintingStretches
  {
  public:
   
    // Constructor
    PrintingStretches(const string MeshFileName, const string OutlineFileName,  const unsigned int NPout, 
		      const string FileNameA, const string FileNameB,  const string FileNameC, const string FileNameD, const unsigned int Type, 
		      const Model::BodyContainer & BDC,
		      const vector<ScalarFieldNode<3>* > & StretchNodes, const vector<ScalarFieldNode<3>* > & DirectionNodes,  
		      const vector<unsigned int> & CapsomersNum,
		      const vector<double > & VectorAngleShift,
		      const vector<double > VectorReferenceAngle,
		      const double HexonShift,
		      const double epsilon,
		      const double AvgEdgeLength,
		      const int refinement,
		      const bool FinalConfiguration,
		      vector<vector<unsigned int> > T9ThreeFold,
		      const unsigned int MT=1): 
    _meshFileName(MeshFileName), _outlineFileName(OutlineFileName), _nPout(NPout),
    _fileNameA(FileNameA), _fileNameB(FileNameB), _fileNameC(FileNameC), _fileNameD(FileNameD), _type(Type),
    _bdc(BDC), _stretchNodes(StretchNodes), _directionNodes(DirectionNodes),
    _capsomersNum(CapsomersNum), _vectorAngleShift(VectorAngleShift), _vectorReferenceAngle(VectorReferenceAngle),
      _hexonShift(HexonShift),  _epsilon(epsilon), _avgEdgeLength(AvgEdgeLength), _refinement(refinement),
      _finalConfiguration(FinalConfiguration), _T9ThreeFold(T9ThreeFold), _MTW(MT) {}
    
    // Destructor
    virtual ~ PrintingStretches () {}

    // Printing functions
    vector<double > printMaster(const int iter)
    {
      // type = 0: Hexagonal lattice
      // type = 1: T virus
	stringstream OutA, OutB;
	string OutFileNameA, OutFileNameB;
	OutA << _fileNameA << iter << ".vtk";
	OutFileNameA = OutA.str();
	OutB << _fileNameB << iter << ".vtk";
	OutFileNameB = OutB.str();
	printLattice(OutFileNameA, OutFileNameB, _type, iter);

	vector<double > SymmIndA, SymmIndB, SymmInd;
	if(_type == 1) { // T virus
	  stringstream OutC, OutD;
	  string OutFileNameC, OutFileNameD;
	  OutC << _fileNameC << iter << ".vtk";
	  OutD << _fileNameD << iter << ".vtk";
	  OutFileNameC = OutC.str();
	  OutFileNameD = OutD.str();
	  
	  SymmIndA = printStretchDirections(OutFileNameC, _epsilon);
	  // SymmIndB = printStretchElemDir(OutFileNameD, _epsilon);

	  // SymmInd.reserve(6);
	  // SymmInd.insert( SymmInd.end(), SymmIndA.begin(), SymmIndA.end() );
	  // SymmInd.insert( SymmInd.end(), SymmIndB.begin(), SymmIndB.end() );
	  //  printStretchEigenvector(OutFileNameD);
	  
	}
	// return SymmInd;
	return SymmIndA;
    }

    void printLattice(const std::string OutFileNameA,
		      const std::string OutFileNameB, 
		      const unsigned int type,
		      const unsigned int iter)
    {
  // OutPut bending and stretching energy in current configuration
  ofstream ofs1(OutFileNameA.c_str());
  if (!ofs1) {
    std::cout << "Error: cannot open output ("
	      << OutFileNameA
	      << ") file." << std::endl;
    exit(0);
  }

  std::string OutlineFile = _outlineFileName+".vtk";
  
  ofstream ofs2(OutFileNameB.c_str());
  if (!ofs2) {
    std::cout << "Error: cannot open output ("
	      << OutFileNameB
	      << ") file." << std::endl;
    exit(0);
  }
  unsigned int ind = 0;

  // Start from bending body which has all elements
  Body::ElementContainer Elements = _bdc[0]->elements();
  Body::NodeContainer Nodes = _bdc[0]->nodes();
  int NumHexon = 0;
  if(_type == 0) {
    NumHexon = _capsomersNum.back()+1;
  }
  
    

  // Node Section
  ofs1 << "# vtk DataFile Version 3.0" << endl
	<< "Test example" << endl
	<< "ASCII" << endl
	<< "DATASET POLYDATA" << endl
        << "POINTS  " << Nodes.size()-2*NumHexon << "  double" << endl;

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

    if (node != NULL) {
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
      } 
    }
  }


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
  }


  // Outline for pentamers and examers
  ifstream ifs(OutlineFile.c_str());
  if (!ifs) {
    cout << "Error: Cannot open input file: " << OutlineFile << endl;
    exit(0);
  }

  string token;
  ifs >> token; 
  while( token != "LINES") ifs >> token; 
  ofs2 << token << " ";
  ifs >> token;
  ofs2 << token << " ";
  ifs >> token;
  ofs2 << token << endl;
  ifs >> token;
  while(!ifs.eof())
  {
    ofs2 << token << " ";
    ifs >> token;
  }
  ifs.close();


  // Cell section
  ofs1 << "CELL_DATA    " << nel << endl;
  
  // Energy of body 1
  if (type == 0) {
    ofs1 << "SCALARS    StretchingEnergy    double    1" << endl;
    ofs1 << "LOOKUP_TABLE default" << endl;
  }
  else if (type == 1) {
    ofs1 << "SCALARS    BendingEnergy    double    1" << endl;
    ofs1 << "LOOKUP_TABLE default" << endl;
  }

  for(int e = 0; e < nel; e++) {
    ofs1 << Elements[e]->energy() << endl;
  }
  ofs1 << std::endl;
  
  // Energy of body > 1 if any
  int nPel = 0;
  if (type == 1) {
    Body::ElementContainer PentElements = _bdc[1]->elements();
    Body::ElementContainer HexElements = _bdc[2]->elements();

    nel = HexElements.size();
    nPel = PentElements.size();

    ofs1 << "SCALARS    StretchingEnergy    double    1" << endl;
    ofs1 << "LOOKUP_TABLE default" << endl;
  
    for(int e = 0; e < nel; e++) {
      ofs1 << HexElements[e]->energy() << endl;
    }

    for(int e = 0; e < nPel; e++) {
      ofs1 << PentElements[e]->energy() << endl;
    }
    ofs1 << endl;
    
    // Conformational energy
    ofs1 << "SCALARS    ConformationalEnergy    double    1" << endl;
    ofs1 << "LOOKUP_TABLE default" << endl;
  
    for(int e = 0; e < nel; e++) {
      C0MembraneStretch * TempE = dynamic_cast<C0MembraneStretch *>(HexElements[e]);
      if (TempE != NULL) {
	 ofs1 << TempE->conformationalEnergy() << endl;
      }
      else {
	ofs1 << 0.0 << endl;
	// cout << "Error in PrintingStretches.h: HexElement does not have conformationalEnergy() " << endl;
	// exit(0);
      }
    }

    for(int e = 0; e < nPel; e++) {
      ofs1 << 0.0 << endl;
    }
    ofs1 << endl;
    

  }


  // Stretch magnitude
  ofs1 << "SCALARS    MaxPrStr    double    1" << endl;
  ofs1 << "LOOKUP_TABLE default" << endl;
  double etaPrint = 0.0; 
  
  if (_stretchNodes.size() == 1) {
    etaPrint = _stretchNodes[0]->point(); 
    for(int e = 0; e < nel; e++) {
      if (etaPrint > 1.0) {
	ofs1 << etaPrint << endl;
      }
      else {
	ofs1 << 1.0/etaPrint << endl;
      }
    }
  }
  else {
    for(int e = 0; e < nel; e++) {
      etaPrint = _stretchNodes[_capsomersNum[e]]->point(); 
      if (etaPrint > 1.0) {
	ofs1 << etaPrint << endl;
      }
      else {
	ofs1 << 1.0/etaPrint << endl;
      }
    }
  }

  for(int e = 0; e < nPel; e++) {
      ofs1 << 1.0 << endl;
  }
  ofs1 << endl;



  ofs1 << "SCALARS    MinPrStr    double    1" << endl;
  ofs1 << "LOOKUP_TABLE default" << endl;
  
  if (_stretchNodes.size() == 1) {
    etaPrint = _stretchNodes[0]->point(); 
    for(int e = 0; e < nel; e++) {
      if (etaPrint < 1.0) {
	ofs1 << etaPrint << endl;
      }
      else {
	ofs1 << 1.0/etaPrint << endl;
      }
    }
  }
  else {
    for(int e = 0; e < nel; e++) {
      etaPrint = _stretchNodes[_capsomersNum[e]]->point(); 
      if (etaPrint < 1.0) {
	ofs1 << etaPrint << endl;
      }
      else {
	ofs1 << 1.0/etaPrint << endl;
      }
    }
  }

  for(int e = 0; e < nPel; e++) {
      ofs1 << 1.0 << endl;
  }
  ofs1 << endl;



  // Eta
  ofs1 << "SCALARS    eta    double    1" << endl;
  ofs1 << "LOOKUP_TABLE default" << endl;
  
  if (_stretchNodes.size() == 1) {
    etaPrint = _stretchNodes[0]->point(); 
    for(int e = 0; e < nel; e++) {
      ofs1 << etaPrint << endl;
    }
  }
  else {
    for(int e = 0; e < nel; e++) {
      etaPrint = _stretchNodes[_capsomersNum[e]]->point();
      ofs1 << etaPrint << endl;
    }
  }

  for(int e = 0; e < nPel; e++) {
      ofs1 << 1.0 << endl;
  }
  ofs1 << endl;



  // Theta
  ofs1 << "SCALARS    theta    double    1" << endl;
  ofs1 << "LOOKUP_TABLE default" << endl;
  
  double thetaPrint = 0.0;
  if (_directionNodes.size() == 1) {
    thetaPrint = _directionNodes[0]->point(); 
    for(int e = 0; e < nel; e++) {
      ofs1 << thetaPrint << endl;
    }
  }
  else {
    for(int e = 0; e < nel; e++) {
      thetaPrint = _directionNodes[_capsomersNum[e]]->point();
      ofs1 << thetaPrint << endl;
    }
  }

  for(int e = 0; e < nPel; e++) {
      ofs1 << 1.0 << endl;
  }
  ofs1 << endl;



  if (_MTW == 1) {
    // Stretch direction
    ofs1 << "SCALARS    StretchDirection    double    1" << endl;
    ofs1 << "LOOKUP_TABLE default" << endl;

    double angle = 0.0, tol = 1.0e-3, PIthird = M_PI/3.0, TwoPIthird = 2.0*M_PI/3.0;
    if (_directionNodes.size() == 1) {
      for(int e = 0; e < nel; e++) {
	// Change reference if needed with _vectorReferenceAngle
	angle = fmod(_directionNodes[0]->point() + _vectorAngleShift[_capsomersNum[e]] - _vectorReferenceAngle[_capsomersNum[e]], M_PI);
	
	if ( (M_PI - fabs(angle)) < tol) {
	  angle = 0.0;
	}
	else if ( fabs(-TwoPIthird - angle) < tol) {
	  angle = PIthird;
	}
	else if ( fabs(TwoPIthird - angle) < tol) {
	  angle = -PIthird;
	}
	
	ofs1 << angle << endl;
      }
    }
    else {
      for(int e = 0; e < nel; e++) {
	angle = fmod(_directionNodes[_capsomersNum[e]]->point() + _vectorAngleShift[_capsomersNum[e]] - _vectorReferenceAngle[_capsomersNum[e]], M_PI);
	
	if ( (M_PI - fabs(angle)) < tol) {
	  angle = 0.0;
	}
	else if ( fabs(-TwoPIthird - angle) < tol) {
	  angle = PIthird;
	}
	else if ( fabs(TwoPIthird - angle) < tol) {
	  angle = -PIthird;
	}
	
	ofs1 << angle << endl;
      }
    }
    
    for(int e = 0; e < nPel; e++) {
      ofs1 << 0.0 << endl;
    }
    ofs1 << endl;
  }



  if (type == 1) {
    Body::ElementContainer PentElements = _bdc[1]->elements();
    Body::ElementContainer HexElements = _bdc[2]->elements();

    nel = HexElements.size();
    nPel = PentElements.size();
    unsigned int indP = 0;
    // Calculate the eigenvalues of the cauchy stress tensor and the 2D invariants of C_hat(stretching part only)
    vector<double > eval1, eval2, eval3;
    vector<double > trChat, Jhat;
    evec1.clear();
    evec2.clear();
    evec3.clear();

    for(int e = 0; e < nel+nPel; e++) {
      if (iter > 0)
      {
	Tensor3D sigma;
	vector<double > invariants(2, 0.0);
	if (e < nel) {
	  C0MembraneStretch * TempEstretch = dynamic_cast<C0MembraneStretch *>(HexElements[e]);
	  if (TempEstretch != NULL) {
	    sigma = TempEstretch->cauchyStress();
	    invariants = TempEstretch->matInvariants();
	    trChat.push_back(invariants[0]);
	    Jhat.push_back(invariants[1]);
	  }
	  else {
	    C0MembraneShear * TempEshear = dynamic_cast<C0MembraneShear *>(HexElements[e]);
	    if (TempEshear != NULL) {
	      sigma = TempEshear->cauchyStress();
	      invariants = TempEshear->matInvariants();
	      trChat.push_back(invariants[0]);
	      Jhat.push_back(invariants[1]);
	    }
	    else {
	      cout << "Something wrong printing cauchy stress principal directions for C0MembraneStretch - C0MembraneShear element." << endl;
	      exit(1);
	    }
	  }
	}
	else {
	  C0Membrane<TriangleQuadrature, EvansElastic, ShapeTri3> * TempEpent = dynamic_cast<C0Membrane<TriangleQuadrature, EvansElastic, ShapeTri3> *>(PentElements[indP]);
	  indP++;
	  if (TempEpent != NULL) {
	    sigma = TempEpent->cauchyStress();
	    invariants = TempEpent->matInvariants();
	    trChat.push_back(invariants[0]);
	    Jhat.push_back(invariants[1]);
	  }
	  else {
	    cout << "Something wrong printing cauchy stress principal directions for C0Membrane element." << endl;
	    exit(1);
	  }
	}
	
	
	// compute Eigenvalues and Eigenvectors by calling LAPACK library
	char jobz = 'V';
	char uplo = 'L';
	int  n    = 3;
	int  lda  = n;
	int  lwork = 3*n-1;
	int  info = 0;
	double A[9];
	double evalues[3];
	double work[lwork];

	// calling lapack function here
	dsyev_(&jobz, &uplo, &n, sigma.data(),&lda, evalues, work, &lwork, &info);
	if (info != 0) {
	  cout << sigma << endl;
	  std::cout << "Error: Something is wrong in DSYEV_ - element = " << e << std::endl;
	  // exit(0);
	}
	double TOL=1e-12;

	Vector3D vec;
	// Stretch direction
	//if(std::fabs(evalues[0])<TOL) {
	vec = row(sigma,0);   evec1.push_back(vec);
	eval1.push_back(evalues[0]);
	vec = row(sigma,1);   evec2.push_back(vec);
	eval2.push_back(evalues[1]);
	vec = row(sigma,2);   evec3.push_back(vec);
	eval3.push_back(evalues[2]);
	
	  
	  /*}
	else if(std::fabs(evalues[1])<TOL) {
	  vec=row(sigma,1);      evec1.push_back(vec);
	  eval1.push_back(evalues[1]);
	  vec=row(sigma,0);      evec2.push_back(vec);
	  eval2.push_back(evalues[0]);
	  vec=row(sigma,2);      evec3.push_back(vec);
	  eval3.push_back(evalues[2]);
	}
	else {
	  vec=row(sigma,2);      evec1.push_back(vec);
	  eval1.push_back(evalues[2]);
	  vec=row(sigma,0);      evec2.push_back(vec);
	  eval2.push_back(evalues[0]);
	  vec=row(sigma,1);      evec3.push_back(vec);
	  eval3.push_back(evalues[1]);
	  }*/
      }
      else {
	Vector3D veczero(0.0);
	evec1.push_back(veczero);
	eval1.push_back(0.0);
	evec2.push_back(veczero);
	eval2.push_back(0.0);
	evec3.push_back(veczero);
	eval3.push_back(0.0);
	trChat.push_back(2.0);
	Jhat.push_back(1.0);
      }
    }
      

    
    ofs1 << "SCALARS    CauchyEigenvalue_1    double    1" << endl;
    ofs1 << "LOOKUP_TABLE default" << endl;
    for(int e = 0; e < nel+nPel; e++) {
      ofs1 << eval1[e] << endl;
    }
    ofs1 << std::endl;
    
    ofs1 << "SCALARS    CauchyEigenvalue_2    double    1" << endl;
    ofs1 << "LOOKUP_TABLE default" << endl;
    for(int e = 0; e < nel+nPel; e++) {
      ofs1 << eval2[e] << endl;
    }
    ofs1 << std::endl;
    
    ofs1 << "SCALARS    CauchyEigenvalue_3    double    1" << endl;
    ofs1 << "LOOKUP_TABLE default" << endl;
    for(int e = 0; e < nel+nPel; e++) {
      ofs1 << eval3[e] << endl;
    }
    ofs1 << std::endl;

    ofs1 << "SCALARS    trChat    double    1" << endl;
    ofs1 << "LOOKUP_TABLE default" << endl;
    for(int e = 0; e < nel+nPel; e++) {
      ofs1 << trChat[e] << endl;
    }
    ofs1 << std::endl;

    ofs1 << "SCALARS    Jhat    double    1" << endl;
    ofs1 << "LOOKUP_TABLE default" << endl;
    for(int e = 0; e < nel+nPel; e++) {
      ofs1 << Jhat[e] << endl;
    }
    ofs1 << std::endl;
  }
    
  
  ofs1.close();
  ofs2.close();

 }










 vector<double > printStretchDirections(const std::string OutFileNameC, double _epsilon)
 {
  ofstream ofs(OutFileNameC.c_str());
  if (!ofs) {
    std::cout << "Error: cannot open output ("
	      << OutFileNameC
	      << ") file." << std::endl;
    exit(0);
  }



  // Use Hexon body which contain all the hexamers
  // Body::ElementContainer ShellElements = _bdc[0]->elements();
  Body::ElementContainer HexElements = _bdc[2]->elements();
  Body::NodeContainer Nodes = _bdc[2]->nodes();

  // As many points as are the hexons 
  int NumHexon = _capsomersNum.back()-11;
  ofs << "# vtk DataFile Version 2.0\n"
      << "Test example" << endl
      << "ASCII" << endl
      << "DATASET POLYDATA" << endl
      << "POINTS  " << 2*NumHexon << "  double" << endl;



  // Compute geometric center of each hexon and vector in the final configuration
  int CurrentCapsomer =  _capsomersNum.front();
  int e = 0, i = 0, j = 0;
  double ind = 0.0, alpha = 0.0;
  Vector3D Hexon(0.0), Nbar(0.0), nbar(0.0), nbarTemp(0.0), NbarPerp(0.0), nbarPerp(0.0), nbarTempPerp(0.0);
  vector<Vector3D > HexonCenters, PreStretchDirHex, PreStretchDirPerpHex;
  for (e = 0; e < HexElements.size(); e++ )
  {
    DeformationNode<3> * nodeA = dynamic_cast<DeformationNode<3>* >(HexElements[e]->baseNodes()[0]);
    DeformationNode<3> * nodeB = dynamic_cast<DeformationNode<3>* >(HexElements[e]->baseNodes()[1]);
    DeformationNode<3> * nodeC = dynamic_cast<DeformationNode<3>* >(HexElements[e]->baseNodes()[2]);

    const Vector3D & nodalPosA =  nodeA->point();
    const Vector3D & nodalPosB =  nodeB->point();
    const Vector3D & nodalPosC =  nodeC->point();
    
   
      Hexon(0) += nodalPosA(0) + nodalPosB(0) + nodalPosC(0);
      Hexon(1) += nodalPosA(1) + nodalPosB(1) + nodalPosC(1);
      Hexon(2) += nodalPosA(2) + nodalPosB(2) + nodalPosC(2);
      ind += 3.0;
      
      // Pre-stretch direction in the reference configuration
      alpha = _vectorAngleShift[CurrentCapsomer] + _directionNodes[CurrentCapsomer]->point();
      Nbar(0) = cos(alpha);
      Nbar(1) = sin(alpha);
      Nbar(2) = 0.0;
      // Perpendicular to pre-stretch direction in the reference configuration
      NbarPerp(0) =-sin(alpha);
      NbarPerp(1) = cos(alpha);
      NbarPerp(2) = 0.0;

      C0MembraneStretch * TempEstretch = dynamic_cast<C0MembraneStretch *>(HexElements[e]);
      // LoopShell<SCElastic > * TempEstretch = dynamic_cast<LoopShell<SCElastic> *>(ShellElements[e]);
      nbarTemp = TempEstretch->PushForwardOperator(Nbar);
      nbarTempPerp = TempEstretch->PushForwardOperator(NbarPerp);
     
      if(tvmet::dot(nbarTemp,nbar)>0) { // parallel
	nbar += nbarTemp;
      } 
      else { // anti-parallel
	nbar -= nbarTemp;
      }

      if(tvmet::dot(nbarTempPerp,nbarPerp)>0) { // parallel
	nbarPerp += nbarTempPerp;
      } 
      else { // anti-parallel
	nbarPerp -= nbarTempPerp;
      }
      
      if (_capsomersNum[e+1] != CurrentCapsomer)
      {
	Hexon(0) /= ind;
	Hexon(1) /= ind;
	Hexon(2) /= ind;
	HexonCenters.push_back(Hexon);
	
	nbar /= norm2(nbar);
	PreStretchDirHex.push_back(nbar);
	nbarPerp /= norm2(nbarPerp);
	PreStretchDirPerpHex.push_back(nbarPerp);
      
      
	Hexon(0) = 0.0;
	Hexon(1) = 0.0;
	Hexon(2) = 0.0;
	ind = 0.0;

	nbar(0) = 0.0;
	nbar(1) = 0.0;
	nbar(2) = 0.0;

	nbarPerp(0) = 0.0;
	nbarPerp(1) = 0.0;
	nbarPerp(2) = 0.0;

	CurrentCapsomer = _capsomersNum[e+1];
      }
      

  }




  // Print hexon position and stretch directions
  for(i = 0; i < NumHexon; i++) {
    ofs << HexonCenters[i](0) << " " << HexonCenters[i](1) << " " << HexonCenters[i](2) << endl
	<< HexonCenters[i](0) << " " << HexonCenters[i](1) << " " << HexonCenters[i](2) << endl;
  }
  ofs << std::endl;
 


  vector<double > Etas;
  vector<Vector3D > Directions;
  // Cell section
  double eta = 0.0;
  ofs << "POINT_DATA    " << 2*NumHexon << endl;
  // Initial stretch directions
  ofs << "VECTORS    MaxStretchDirections    double" << endl;
  for(i = 0; i < NumHexon; i++) {
    
    if (_stretchNodes.size() == 1) {
      eta = _stretchNodes[0]->point(); }
    else {
      eta = _stretchNodes[i]->point(); };

    if ( fabs(1.0-eta) < _epsilon*1.2) {
      ofs << 0.0 << " " <<0.0  << " " << 0.0 << endl
	  << 0.0 << " " <<0.0  << " " << 0.0 << endl;
      if ( eta > 1.0) {	
	Etas.push_back(eta);
	Directions.push_back(PreStretchDirHex[i]); 
      }
      else { 
	Etas.push_back(1.0/eta);
	Directions.push_back(PreStretchDirPerpHex[i]); 
      }
      
    }
    else if ( eta > 1.0) {
      ofs << PreStretchDirHex[i](0)  << " " << PreStretchDirHex[i](1)  << " " << PreStretchDirHex[i](2) << endl
	  << -PreStretchDirHex[i](0) << " " << -PreStretchDirHex[i](1) << " " << -PreStretchDirHex[i](2) << endl;
      Etas.push_back(eta);
      Directions.push_back(PreStretchDirHex[i]);
    }
    else {
      ofs << PreStretchDirPerpHex[i](0)  << " " << PreStretchDirPerpHex[i](1)  << " " << PreStretchDirPerpHex[i](2) << endl
	  << -PreStretchDirPerpHex[i](0) << " " << -PreStretchDirPerpHex[i](1) << " " << -PreStretchDirPerpHex[i](2) << endl;
      Etas.push_back(1.0/eta);
      Directions.push_back(PreStretchDirPerpHex[i]);
    }
	   
  }
  ofs << endl;



  ofs << "VECTORS    EtaDirections    double" << endl;
  for(i = 0; i < NumHexon; i++) {
    
    if (_stretchNodes.size() == 1) {
      eta = _stretchNodes[0]->point(); }
    else {
      eta = _stretchNodes[i]->point(); };

    if ( fabs(1.0-eta) < _epsilon*1.2) {
      ofs << 0.0 << " " <<0.0  << " " << 0.0 << endl
	  << 0.0 << " " <<0.0  << " " << 0.0 << endl;
    }
    else {
      ofs << PreStretchDirHex[i](0)  << " " << PreStretchDirHex[i](1)  << " " << PreStretchDirHex[i](2) << endl
	  << -PreStretchDirHex[i](0) << " " << -PreStretchDirHex[i](1) << " " << -PreStretchDirHex[i](2) << endl;
    }
    	   
  }
  ofs << endl;

  ofs << "VECTORS    warp    double" << endl;
  for(i = 0; i < NumHexon; i++) {
    ofs << HexonCenters[i](0) << " " << HexonCenters[i](1) << " " << HexonCenters[i](2) << endl
	<< HexonCenters[i](0) << " " << HexonCenters[i](1) << " " << HexonCenters[i](2) << endl;
  }
  ofs << endl;

  
  ofs.close();

  return ComputeSymmIndex(HexonCenters, Etas, Directions, _avgEdgeLength*_refinement*sqrt(3.0)*0.5);
 }










vector<double > printStretchElemDir(const std::string OutFileNameD, double _epsilon)
{
  ofstream ofs(OutFileNameD.c_str());
  if (!ofs) {
    std::cout << "Error: cannot open output ("
	      << OutFileNameD
	      << ") file." << std::endl;
    exit(0);
  }



  // Use Hexon body which contain all the hexamers
  Body::ElementContainer HexElements = _bdc[2]->elements();

  // As many points as are the hexons 
  int NumElem = HexElements.size();
  ofs << "# vtk DataFile Version 2.0\n"
      << "Test example" << endl
      << "ASCII" << endl
      << "DATASET POLYDATA" << endl
      << "POINTS  " << 2*NumElem << "  double" << endl;



  // Compute geometric center of each hexon and vector in the final configuration
  int e = 0, i = 0, j = 0;
  double alpha = 0.0;
  Vector3D ElemCenter(0.0), Nbar(0.0), nbar(0.0), NbarPerp(0.0), nbarPerp(0.0);
  vector<Vector3D > ElemCenters, PreStretchDirEl, PreStretchDirPerpEl;
  for (e = 0; e < NumElem; e++)
  {
    DeformationNode<3> * nodeA = dynamic_cast<DeformationNode<3>* >(HexElements[e]->baseNodes()[0]);
    DeformationNode<3> * nodeB = dynamic_cast<DeformationNode<3>* >(HexElements[e]->baseNodes()[1]);
    DeformationNode<3> * nodeC = dynamic_cast<DeformationNode<3>* >(HexElements[e]->baseNodes()[2]);

    const Vector3D & nodalPosA =  nodeA->point();
    const Vector3D & nodalPosB =  nodeB->point();
    const Vector3D & nodalPosC =  nodeC->point();
    
   
    ElemCenter(0) = (nodalPosA(0) + nodalPosB(0) + nodalPosC(0))/3.0;
    ElemCenter(1) = (nodalPosA(1) + nodalPosB(1) + nodalPosC(1))/3.0;
    ElemCenter(2) = (nodalPosA(2) + nodalPosB(2) + nodalPosC(2))/3.0;

    ElemCenters.push_back(ElemCenter);
      
      // Pre-stretch direction in the reference configuration
      alpha = _vectorAngleShift[_capsomersNum[e]] + _directionNodes[_capsomersNum[e]]->point();
      Nbar(0) = cos(alpha);
      Nbar(1) = sin(alpha);
      Nbar(2) = 0.0;
      // Perpendicular to pre-stretch direction in the reference configuration
      NbarPerp(0) = -sin(alpha);
      NbarPerp(1) = cos(alpha);
      NbarPerp(2) = 0.0;

      C0MembraneStretch * TempEstretch = dynamic_cast<C0MembraneStretch *>(HexElements[e]);
      nbar = TempEstretch->PushForwardOperator(Nbar);
      nbarPerp = TempEstretch->PushForwardOperator(NbarPerp);
           
      nbar /= norm2(nbar);
      PreStretchDirEl.push_back(nbar);
      nbarPerp /= norm2(nbarPerp);
      PreStretchDirPerpEl.push_back(nbarPerp);    

  }




  // Print hexon position and stretch directions
  for(e = 0; e < NumElem; e++) {
    ofs << ElemCenters[e](0) << " " << ElemCenters[e](1) << " " << ElemCenters[e](2) << endl
	<< ElemCenters[e](0) << " " << ElemCenters[e](1) << " " << ElemCenters[e](2) << endl;
  }
  ofs << std::endl;
 

  
  vector<double > Etas;
  vector<Vector3D > Directions;
  // Cell section
  double eta = 0.0;
  ofs << "POINT_DATA    " << 2*NumElem << endl;
  // Initial stretch directions
  ofs << "VECTORS    MaxStretchDirections    double" << endl;
  for(e = 0; e < NumElem; e++) {
    
    if (_stretchNodes.size() == 1) {
      eta = _stretchNodes[0]->point(); }
    else {
      eta = _stretchNodes[_capsomersNum[e]]->point(); };

    if ( fabs(1.0-eta) < _epsilon*1.2) {
      ofs << 0.0 << " " << 0.0  << " " << 0.0 << endl
	  << 0.0 << " " << 0.0  << " " << 0.0 << endl;
      if ( eta > 1.0) {	
	Etas.push_back(eta); 
	Directions.push_back(PreStretchDirEl[e]);
	  }
      else { 
	Etas.push_back(1.0/eta); 
	Directions.push_back(PreStretchDirPerpEl[e]);
      }
    }
    else if ( eta > 1.0) {
      ofs <<  PreStretchDirEl[e](0) << " " <<  PreStretchDirEl[e](1) << " " <<  PreStretchDirEl[e](2) << endl
	  << -PreStretchDirEl[e](0) << " " << -PreStretchDirEl[e](1) << " " << -PreStretchDirEl[e](2) << endl;
      Etas.push_back(eta);
      Directions.push_back(PreStretchDirEl[e]);
    }
    else {
      ofs <<  PreStretchDirPerpEl[e](0) << " " <<  PreStretchDirPerpEl[e](1) << " " <<  PreStretchDirPerpEl[e](2) << endl
	  << -PreStretchDirPerpEl[e](0) << " " << -PreStretchDirPerpEl[e](1) << " " << -PreStretchDirPerpEl[e](2) << endl;
      Etas.push_back(1.0/eta);
      Directions.push_back(PreStretchDirPerpEl[e]);
    }
	   
  }
  ofs << endl;



  ofs << "VECTORS    EtaDirections    double" << endl;
  for(e = 0; e < NumElem; e++) {
    
    if (_stretchNodes.size() == 1) {
      eta = _stretchNodes[0]->point(); }
    else {
      eta = _stretchNodes[_capsomersNum[e]]->point(); };

    if ( fabs(1.0-eta) < _epsilon*1.2) {
      ofs << 0.0 << " " <<0.0  << " " << 0.0 << endl
	  << 0.0 << " " <<0.0  << " " << 0.0 << endl;
    }
    else {
      ofs <<  PreStretchDirEl[e](0)  << " "<<  PreStretchDirEl[e](1)  << " "<<  PreStretchDirEl[e](2) << endl
	  << -PreStretchDirEl[e](0) << " " << -PreStretchDirEl[e](1) << " " << -PreStretchDirEl[e](2) << endl;
    }
    	   
  }
  ofs << endl;
  
  ofs << "VECTORS    warp    double" << endl;
  for(e = 0; e < NumElem; e++) {
    ofs << ElemCenters[e](0) << " " << ElemCenters[e](1) << " " << ElemCenters[e](2) << endl
	<< ElemCenters[e](0) << " " << ElemCenters[e](1) << " " << ElemCenters[e](2) << endl;
  }
  ofs << endl;
  

 
  ofs.close();



  return ComputeSymmIndex(ElemCenters, Etas, Directions, _avgEdgeLength*0.2);

 }










void printStretchEigenvector(const std::string OutFileNameD)
 {
  ofstream ofs(OutFileNameD.c_str());
  if (!ofs) {
    std::cout << "Error: cannot open output ("
	      << OutFileNameD
	      << ") file." << std::endl;
    exit(0);
  }

  // Start from bending body which has all elements
  Body::ElementContainer Elements = _bdc[0]->elements();
  Body::NodeContainer Nodes = _bdc[0]->nodes();

  // As many points as are the hexons 
  int NumCapsomers = _capsomersNum.back()+1;
  ofs << "# vtk DataFile Version 2.0\n"
      << "Test example" << endl
      << "ASCII" << endl
      << "DATASET POLYDATA" << endl
      << "POINTS  " << NumCapsomers << "  double" << endl;

  // Compute geometric center of each hexon (and pentons)
  int CurrentCapsomer =  _capsomersNum.front();
  int i = -1, j = 0;
  double ind = 0.0;
  vector<Vector3D > HexonCenters;
  Vector3D Hexon(0.0);
  while (CurrentCapsomer < NumCapsomers)
  {
    i++;
    DeformationNode<3> * nodeA = dynamic_cast<DeformationNode<3>* >(Elements[i]->baseNodes()[0]);
    DeformationNode<3> * nodeB = dynamic_cast<DeformationNode<3>* >(Elements[i]->baseNodes()[1]);
    DeformationNode<3> * nodeC = dynamic_cast<DeformationNode<3>* >(Elements[i]->baseNodes()[2]);

    const Vector3D & nodalPosA =  nodeA->point();
    const Vector3D & nodalPosB =  nodeB->point();
    const Vector3D & nodalPosC =  nodeC->point();
    
    if (_capsomersNum[i] == CurrentCapsomer)
    {
      Hexon(0) += nodalPosA(0) + nodalPosB(0) + nodalPosC(0);
      Hexon(1) += nodalPosA(1) + nodalPosB(1) + nodalPosC(1);
      Hexon(2) += nodalPosA(2) + nodalPosB(2) + nodalPosC(2);
      ind += 3.0;
    }
    else
    {
      Hexon(0) /= ind;
      Hexon(1) /= ind;
      Hexon(2) /= ind;
      HexonCenters.push_back(Hexon);
      
      Hexon(0) = nodalPosA(0) + nodalPosB(0) + nodalPosC(0);
      Hexon(1) = nodalPosA(1) + nodalPosB(1) + nodalPosC(1);
      Hexon(2) = nodalPosA(2) + nodalPosB(2) + nodalPosC(2);
      ind = 3.0;
      CurrentCapsomer = _capsomersNum[i];
    }
  }

  // Print hexon position and stretch directions
  for(i = 0; i < NumCapsomers; i++) {
    ofs << HexonCenters[i](0) << " " 
	<< HexonCenters[i](1) << " "
	<< HexonCenters[i](2) << endl;
  }
  ofs << std::endl;
 
  // Cell section
  ofs << "POINT_DATA    " << NumCapsomers << endl;
  ofs << "VECTORS    StretchEigenvector_1    double" << endl;

  for(i = 0; i < NumCapsomers; i++) {
    ofs << evec1[i](0) << " " 
	<< evec1[i](1) << " "
	<< evec1[i](2) << endl;
  }
  ofs << endl;

  ofs << "POINT_DATA    " << NumCapsomers << endl;
  ofs << "VECTORS    StretchEigenvector_2    double" << endl;

  for(i = 0; i < NumCapsomers; i++) {
    ofs << evec2[i](0) << " " 
	<< evec2[i](1) << " "
	<< evec2[i](2) << endl;
  }
  ofs << endl;

  ofs << "POINT_DATA    " << NumCapsomers << endl;
  ofs << "VECTORS    StretchEigenvector_3    double" << endl;

  for(i = 0; i < NumCapsomers; i++) {
    ofs << evec3[i](0) << " " 
	<< evec3[i](1) << " "
	<< evec3[i](2) << endl;
  }
  ofs << endl;
  
  ofs.close();

 }


    
    void setPrintingConfig(const bool PrintingConfig)
    {
      _finalConfiguration = PrintingConfig;
    }



    Vector3D Rotate(Vector3D A, vector<double > R)
    {
      Vector3D B(0.0);
  
      B(0) = R[0]*A(0) + R[1]*A(1) + R[2]*A(2);
      B(1) = R[3]*A(0) + R[4]*A(1) + R[5]*A(2);
      B(2) = R[6]*A(0) + R[7]*A(1) + R[8]*A(2);
  
      return B;
    }



    double crossProductAngle(Vector3D A, Vector3D B)
    {
      double normA = sqrt( pow(A(0), 2.0) + pow(A(1), 2.0) + pow(A(2), 2.0) );
      double normB = sqrt( pow(B(0), 2.0) + pow(B(1), 2.0) + pow(B(2), 2.0) );
      return fabs( asin( sqrt( pow(A(1)*B(2) - A(2)*B(1), 2.0) +
			       pow(A(2)*B(0) - A(0)*B(2), 2.0) +
			       pow(A(0)*B(1) - A(1)*B(0), 2.0) )
			 /(normA*normB)
			 )
		   )*(180/M_PI); 
    }



    vector<double > ComputeSymmIndex(vector<Vector3D > & Centers, vector<double > & Etas,  vector<Vector3D > & Directions, double tolCenter)
    {
      vector<double > SymmIndex(3, 0.0);
      int i = 0, j = 0;
      double tolNode = _avgEdgeLength*0.1;

      string RotationFile = "RotationMatrices.dat";
      // Read in rotations file
      ifstream ifs(RotationFile.c_str(), ios::in);
      if (!ifs) {
	cout << "Cannot open rotation matrices file: " << RotationFile << endl;
	exit(0);
      }

      vector<vector<double > > Rotations(60, vector<double > (9, 0.0));
      for (i=0; i<60; i++) {
	for (j=0; j<9; j++) {
	  ifs >> Rotations[i][j];
	}
      }
      ifs.close();

      // All necessary data have been read in at this point

      // Compute deviation from symmetry
      Body::NodeContainer Nodes = _bdc[0]->nodes(); // Nodes of LoopShell body
      int NumNodes = Nodes.size();
      unsigned int k = 0, Mapped = 0, indk = 0, npts = Centers.size();
      Vector3D CenterRotated(0.0);
      Vector3D DirectionRotated(0.0);
      Vector3D NodeRotated(0.0);
      for (i=0; i<60; i++)
	{
	  // Rotate centers
	  for  (j=0; j<npts; j++) 
	    {
	      Mapped = 0;
	      indk = 0;
	      CenterRotated = Rotate(Centers[j],  Rotations[i]);
	      
	      DirectionRotated = Rotate(Directions[j],  Rotations[i]);
	      
	      // Find corresponding rotated center
	      for (k=0; k<npts; k++)
		{
		  if ( sqrt( pow(CenterRotated(0) - Centers[k](0), 2.0) +
			     pow(CenterRotated(1) - Centers[k](1), 2.0) + 
			     pow(CenterRotated(2) - Centers[k](2), 2.0) ) < tolCenter )
		    {
		      indk = k;
		      Mapped++;
		    }
		}
	      
	      if (Mapped != 1) {
		cout << "Error with center mapping. Mapped = " << Mapped << " at i = " << i << endl;
		// exit(0);
	      }
	      
	      // Check if corresponding hexamer has the same stretch magnitude and direction
	      SymmIndex[1] += fabs(Etas[j] - Etas[indk]);
	      SymmIndex[2] += fabs(1.0-Etas[j])*crossProductAngle(DirectionRotated, Directions[indk]);
	    }

	  /*
	  // Rotate nodes
	  for  (j=0; j<NumNodes; j++) 
	    {
	      Mapped = 0;
	      indk = 0;
	      DeformationNode<3> * nodeA = dynamic_cast<DeformationNode<3>* >(Nodes[j]);
    
	      NodeRotated = Rotate(nodeA->point(),  Rotations[i]);
	      
	      // Find corresponding rotated center
	      for (k=0; k<NumNodes; k++)
	      {
		DeformationNode<3> * nodeK = dynamic_cast<DeformationNode<3>* >(Nodes[k]);
		if ( sqrt( pow(NodeRotated(0) - (nodeK->point())(0), 2.0) +
			   pow(NodeRotated(1) - (nodeK->point())(1), 2.0) + 
			   pow(NodeRotated(2) - (nodeK->point())(2), 2.0) ) < tolNode )
		  {
		    indk = k;
		    Mapped++;
		  }
	      }
	      
	      if (Mapped != 1) {
		cout << "Error with node mapping. Mapped = " << Mapped << " at i = " << i << endl;
		// exit(0);
	      }
	      DeformationNode<3> * nodeK = dynamic_cast<DeformationNode<3>* >(Nodes[indk]);
	      SymmIndex[0] +=  sqrt( pow(NodeRotated(0) - (nodeK->point())(0), 2.0) +
				     pow(NodeRotated(1) - (nodeK->point())(1), 2.0) + 
				     pow(NodeRotated(2) - (nodeK->point())(2), 2.0) );
	    }
	  */
	  
	}
    
      
      SymmIndex[0] /= double(60*NumNodes);
      SymmIndex[1] /= 60.0;
      SymmIndex[2] /= 60.0;
      
      return SymmIndex;
    }

    vector<int > FerroMagneticInteractions()
    {
      vector<int > Interactions(2, 0);
      Body::ElementContainer HexElements = _bdc[2]->elements();
      Body::NodeContainer Nodes = _bdc[2]->nodes();

      // As many points as are the hexons 
      int NumHexon = _capsomersNum.back()-11;

      // Compute geometric center of each hexon and vector in the final configuration
      int CurrentCapsomer =  _capsomersNum.front();
      int e = 0, i = 0, j = 0;
      double ind = 0.0, alpha = 0.0;
      Vector3D Hexon(0.0), Nbar(0.0), nbar(0.0), nbarTemp(0.0), NbarPerp(0.0), nbarPerp(0.0), nbarTempPerp(0.0);
      vector<Vector3D > HexonCenters, PreStretchDirHex, PreStretchDirPerpHex;
      for (e = 0; e < HexElements.size(); e++ )
	{
	  DeformationNode<3> * nodeA = dynamic_cast<DeformationNode<3>* >(HexElements[e]->baseNodes()[0]);
	  DeformationNode<3> * nodeB = dynamic_cast<DeformationNode<3>* >(HexElements[e]->baseNodes()[1]);
	  DeformationNode<3> * nodeC = dynamic_cast<DeformationNode<3>* >(HexElements[e]->baseNodes()[2]);

	  const Vector3D & nodalPosA =  nodeA->point();
	  const Vector3D & nodalPosB =  nodeB->point();
	  const Vector3D & nodalPosC =  nodeC->point();
    
   
	  Hexon(0) += nodalPosA(0) + nodalPosB(0) + nodalPosC(0);
	  Hexon(1) += nodalPosA(1) + nodalPosB(1) + nodalPosC(1);
	  Hexon(2) += nodalPosA(2) + nodalPosB(2) + nodalPosC(2);
	  ind += 3.0;
      
	  // Pre-stretch direction in the reference configuration
	  alpha = _vectorAngleShift[CurrentCapsomer] + _directionNodes[CurrentCapsomer]->point();
	  Nbar(0) = cos(alpha);
	  Nbar(1) = sin(alpha);
	  Nbar(2) = 0.0;
	  // Perpendicular to pre-stretch direction in the reference configuration
	  NbarPerp(0) =-sin(alpha);
	  NbarPerp(1) = cos(alpha);
	  NbarPerp(2) = 0.0;

	  C0MembraneStretch * TempEstretch = dynamic_cast<C0MembraneStretch *>(HexElements[e]);
	  // LoopShell<SCElastic > * TempEstretch = dynamic_cast<LoopShell<SCElastic> *>(ShellElements[e]);
	  nbarTemp = TempEstretch->PushForwardOperator(Nbar);
	  nbarTempPerp = TempEstretch->PushForwardOperator(NbarPerp);
	  
	  if(tvmet::dot(nbarTemp,nbar)>0) { // parallel
	    nbar += nbarTemp;
	  } 
	  else { // anti-parallel
	    nbar -= nbarTemp;
	  }
	  
	  if(tvmet::dot(nbarTempPerp,nbarPerp)>0) { // parallel
	    nbarPerp += nbarTempPerp;
	  } 
	  else { // anti-parallel
	    nbarPerp -= nbarTempPerp;
	  }
      
	  if (_capsomersNum[e+1] != CurrentCapsomer)
	    {
	      Hexon(0) /= ind;
	      Hexon(1) /= ind;
	      Hexon(2) /= ind;
	      HexonCenters.push_back(Hexon);
	
	      nbar /= norm2(nbar);
	      PreStretchDirHex.push_back(nbar);
	      nbarPerp /= norm2(nbarPerp);
	      PreStretchDirPerpHex.push_back(nbarPerp);
      
	      
	      Hexon(0) = 0.0;
	      Hexon(1) = 0.0;
	      Hexon(2) = 0.0;
	      ind = 0.0;
	      
	      nbar(0) = 0.0;
	      nbar(1) = 0.0;
	      nbar(2) = 0.0;
	      
	      nbarPerp(0) = 0.0;
	      nbarPerp(1) = 0.0;
	      nbarPerp(2) = 0.0;

	      CurrentCapsomer = _capsomersNum[e+1];
	    }
     
	}



      vector<double > Etas;
      vector<Vector3D > Directions;
     
      double eta = 0.0;
      
      for(i = 0; i < NumHexon; i++) {
    
	if (_stretchNodes.size() == 1) {
	  eta = _stretchNodes[0]->point(); }
	else {
	  eta = _stretchNodes[i]->point(); };

	if ( eta > 1.0) {	
	  Directions.push_back(PreStretchDirHex[i]); 
	}
	else { 
	  Directions.push_back(PreStretchDirPerpHex[i]); 
	}
	
      }

    


      // for (i = 0; i<_T9ThreeFold.size(); i++)
      for (i = 0; i<20; i++)
      {
	Vector3D Vi, Di;
	Vi = HexonCenters[i]/norm2(HexonCenters[i]); Di = Directions[i];
	//for (j = 0; j<(_T9ThreeFold[i]).size(); j++)
	for (j = 0; j<3; j++)
	{
	  Vector3D Vj, Dj;
	  Vj = HexonCenters[j]/norm2(HexonCenters[j]); Dj = Directions[j];
	  Vector3D k = crossProduct(Vi, Vj);
	  k /= norm2(k);

	  double theta = acos(dotProduct(Vi, Vj));

	  Tensor3D Omega(0.0), Rotation(0.0), Identity(0.0);
	  Omega(0,0) = 0.0;  Omega(0,1) =-k(2); Omega(0,2) = k(2);
	  Omega(1,0) = k(3); Omega(1,1) = 0.0;  Omega(1,2) =-k(1); 
	  Omega(2,0) =-k(2); Omega(2,1) = k(1); Omega(2,2) = 0.0;

	  Identity(0,0) = 1.0; Identity(1,1) = 1.0; Identity(2,2) = 1.0; 

	  Rotation = Identity + Omega*sin(theta) + Omega*Omega*(1-cos(theta));

	  Vector3D DiRotated;
	  DiRotated = Rotation*Di;

	  double check =  acos(dotProduct(DiRotated,Dj));
	  if (check > (M_PI/6.0) && check < (5.0*M_PI/6.0)) {
	    Interactions[1] +=1;
	  }
	  else {
	    Interactions[0] +=1;
	  }
	   
	}

      }

      return Interactions;

    }





    Vector3D crossProduct(Vector3D Vi, Vector3D Vj)
    {
      Vector3D Vk(0.0);
      Vk(0) = Vi(1)*Vj(2) - Vi(2)*Vj(1);
      Vk(1) = Vi(2)*Vj(0) - Vi(0)*Vj(2);
      Vk(2) = Vi(0)*Vj(1) - Vi(1)*Vj(0);
     
      return Vk;
    }

    double dotProduct(Vector3D Vi, Vector3D Vj)
    {    
      return Vi(0)*Vj(0) + Vi(1)*Vj(1) + Vi(2)*Vj(2);
    }

    

    



  private:	

    const string _meshFileName;
    const string _outlineFileName;  
    const unsigned int _nPout; 

    const string _fileNameA;
    const string _fileNameB;
    const string _fileNameC;
    const string _fileNameD;
    const unsigned int _type;
    const unsigned int _MTW;

    const Model::BodyContainer & _bdc;
    const vector<ScalarFieldNode<3>* > & _stretchNodes;
    const vector<ScalarFieldNode<3>* > & _directionNodes;  
    const vector<unsigned int> & _capsomersNum;

    const vector<double> & _vectorAngleShift;
    const vector<double> _vectorReferenceAngle;
    const double _hexonShift;
    const double _epsilon;
    const double _avgEdgeLength;
    const int _refinement;
    bool _finalConfiguration;
    vector<vector<unsigned int> > & _T9ThreeFold;

    vector<Vector3D> evec1, evec2, evec3;

  };
  
}; // namespace voom

#endif // __PrintingStretches_h__


