
//----------------------------------------------------------------------
//                   Luigi Perotti, William S. Klug
//                University of California Los Angeles
//                   (C) 2004 All Rights Reserved
//----------------------------------------------------------------------

#if !defined(__PrintingStretches_h__)
#define __PrintingStretches_h__

#include <string>
#include <iostream>
#include <vector>
#include <fstream>


#include "C0MembraneShear.h"
#include "C0MembraneStretch.h"
#include "C0Membrane.h"
#include "EvansElastic.h"
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
		      const vector<double > VectorAngleShift,
		      const vector<double > VectorReferenceAngle,
		      const double HexonShift,
		      const bool FinalConfiguration,
		      const unsigned int MT=1): 
    _meshFileName(MeshFileName), _outlineFileName(OutlineFileName), _nPout(NPout),
    _fileNameA(FileNameA), _fileNameB(FileNameB), _fileNameC(FileNameC), _fileNameD(FileNameD), _type(Type),
    _bdc(BDC), _stretchNodes(StretchNodes), _directionNodes(DirectionNodes),
    _capsomersNum(CapsomersNum), _vectorAngleShift(VectorAngleShift), _vectorReferenceAngle(VectorReferenceAngle), _hexonShift(HexonShift),
    _finalConfiguration(FinalConfiguration), _MT(MT) {}
    
    // Destructor
    virtual ~ PrintingStretches () {}

    // Printing functions
    void printMaster(const int iter)
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

	if(_type == 1) { // T virus
	  stringstream OutC, OutD;
	  string OutFileNameC, OutFileNameD;
	  OutC << _fileNameC << iter << ".vtk";
	  OutD << _fileNameD << iter << ".vtk";
	  OutFileNameC = OutC.str();
	  OutFileNameD = OutD.str();
	  
	  printStretchDirections(OutFileNameC);
	  //  printStretchEigenvector(OutFileNameD);
	  
	}
    }

    void printLattice(const std::string OutFileNameA,
		      const std::string OutFileNameB, 
		      const unsigned int type,
		      const unsigned int iter)
    {
  // OutPut bending and stretching energy in current configuration
  ofstream ofs1(OutFileNameA.c_str());
  if (!ofs1) {
    std::cout << "can not open output ("
	      << OutFileNameA
	      << ") file." << std::endl;
    exit(0);
  }

  std::string OutlineFile = _outlineFileName+".vtk";
  
  ofstream ofs2(OutFileNameB.c_str());
  if (!ofs2) {
    std::cout << "can not open output ("
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
    cout << "Cannot open input file: " << OutlineFile << endl;
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



  if (_MT == 1) {
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
	  std::cout << "Something is wrong in DSYEV_ - element = " << e << std::endl;
	  exit(0);
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



 void printStretchDirections(const std::string OutFileNameC)
 {
  ofstream ofs(OutFileNameC.c_str());
  if (!ofs) {
    std::cout << "can not open output ("
	      << OutFileNameC
	      << ") file." << std::endl;
    exit(0);
  }

  // Start from bending body which has all elements
  Body::ElementContainer Elements = _bdc[0]->elements();
  Body::NodeContainer Nodes = _bdc[0]->nodes();

  // As many points as are the hexons 
  int NumHexon = _capsomersNum.back()-11;
  ofs << "# vtk DataFile Version 2.0\n"
      << "Test example" << endl
      << "ASCII" << endl
      << "DATASET POLYDATA" << endl
      << "POINTS  " << 2*NumHexon << "  double" << endl;
   // << "POINTS  " << 4*NumHexon << "  double" << endl;

  // Compute geometric center of each hexon (and pentons)
  int CurrentCapsomer =  _capsomersNum.front();
  int i = -1, j = 0;
  double ind = 0.0;
  vector<Vector3D > HexonCenters;
  Vector3D Hexon(0.0);
  while (CurrentCapsomer < NumHexon)
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

  



  // Read in the initial stretch directions\
  // Create input stream
  ifstream ifs;
  ifs.open( _meshFileName.c_str(), ios::in);
  if (!ifs) {
    cout << "Cannot open input file: " << _meshFileName << endl;
    exit(0);
  }
  // Create vector of initial stretch directions
  vector<tvmet::Vector<double,3> > InitialStretchDirections;
  vector<tvmet::Vector<double,3> > RotationAxis;
  
  string token;
  ifs >> token; 
  while( token != "shear_direction") ifs >> token; 
  ifs >> token;

  int SavedCapsomer = -1;
  CurrentCapsomer = _capsomersNum.front();
  i = 0;
  while (CurrentCapsomer < NumHexon)
  {
    if (CurrentCapsomer != SavedCapsomer)
    {
      tvmet::Vector<double,3> StretchDir;
      ifs >>  StretchDir(0) >> StretchDir(1) >> StretchDir(2);
      StretchDir = StretchDir/tvmet::norm2(StretchDir);
      InitialStretchDirections.push_back(StretchDir);
      // cout << StretchDir << CurrentCapsomer << endl;
      SavedCapsomer = CurrentCapsomer;

      // Computing axis of rotation
      DeformationNode<3> * nodeA = dynamic_cast<DeformationNode<3>* >(Elements[i]->baseNodes()[0]);
      DeformationNode<3> * nodeB = dynamic_cast<DeformationNode<3>* >(Elements[i]->baseNodes()[1]);
      DeformationNode<3> * nodeC = dynamic_cast<DeformationNode<3>* >(Elements[i]->baseNodes()[2]);

      tvmet::Vector<double,3> e31(nodeA->position()-nodeC->position()), 
 	                      e32(nodeB->position()-nodeC->position());
      
      tvmet::Vector<double,3> n = tvmet::cross(e31,e32);
      n = n/tvmet::norm2(n);

      RotationAxis.push_back(n);
    }
    else
    {
      ifs >> token >> token >> token;
    }
    i++;
    CurrentCapsomer = _capsomersNum[i];
    // cout << i << " " << CurrentCapsomer << " " << SavedCapsomer << endl;
  }





  // Compute stretch directions and positions
  // vector<Vector3D > HexonPoint(4*NumHexon, Vector3D(0.0)),  Stretch(4*NumHexon, Vector3D(0.0));
  vector<Vector3D > HexonPoint(2*NumHexon, Vector3D(0.0)),  Stretch(2*NumHexon, Vector3D(0.0)), 
                    StretchScaled(2*NumHexon, Vector3D(0.0)), EtaDir(2*NumHexon, Vector3D(0.0));
  double beta = 0.0, eta = 0.0, factEta = 0.0;
  for(i = 0; i < NumHexon; i++) {
    // Imposed initial angle shift (IN RADIANTS!!)
    if (_finalConfiguration) {
      if (_directionNodes.size() == 1) {
	beta = _vectorAngleShift[i] + _directionNodes[0]->point();
      }
      else {
	beta = _vectorAngleShift[i] + _directionNodes[i]->point();
      }
    }
    else {
      beta = _vectorAngleShift[i];
    }

    tvmet::Vector<double, 3> StretchLoc, StretchLocP;
    StretchLoc = InitialStretchDirections[i]*cos(beta) + tvmet::cross(RotationAxis[i],InitialStretchDirections[i])*sin(beta);
    StretchLoc /= tvmet::norm2(StretchLoc);


    // Move point of application of InitStretch in the Yloc direction
    StretchLocP = tvmet::cross(RotationAxis[i],StretchLoc);
    if (_stretchNodes.size() == 1) {
      eta = _stretchNodes[0]->point();
    }
    else {
      eta = _stretchNodes[i]->point();
    }
    /*
    StretchLoc *= eta;
    StretchLocP /= eta; */



    // if (eta > 1.0-1.0e-12) {
    if (eta > 1.0) {
      HexonPoint[i] = HexonCenters[i];
      HexonPoint[i+NumHexon] = HexonCenters[i];
      // HexonPoint[i+2*NumHexon] = HexonCenters[i] + _hexonShift*StretchLocP;
      // HexonPoint[i+3*NumHexon] = HexonCenters[i] - _hexonShift*StretchLocP;

      Stretch[i] = StretchLoc;
      Stretch[i+NumHexon] = -StretchLoc;
      // Stretch[i+2*NumHexon] = -StretchLocP;
      // Stretch[i+3*NumHexon] =  StretchLocP;
      factEta = eta-1.0;
      StretchScaled[i] = StretchLoc*factEta;
      StretchScaled[i+NumHexon] = -StretchLoc*factEta;
    }
    else {
      HexonPoint[i] = HexonCenters[i];
      HexonPoint[i+NumHexon] = HexonCenters[i];
      // HexonPoint[i+2*NumHexon] = HexonCenters[i] + _hexonShift*StretchLoc;
      // HexonPoint[i+3*NumHexon] = HexonCenters[i] - _hexonShift*StretchLoc;

      Stretch[i] = StretchLocP;
      Stretch[i+NumHexon] = -StretchLocP;
      // Stretch[i+2*NumHexon] = -StretchLoc;
      // Stretch[i+3*NumHexon] =  StretchLoc;
      factEta = (1.0/eta)-1.0;
      StretchScaled[i] = StretchLocP*factEta;
      StretchScaled[i+NumHexon] = -StretchLocP*factEta;
    }

    EtaDir[i] = StretchLoc;
    EtaDir[i+NumHexon] = -StretchLoc;
    
  }



  // Print hexon position and stretch directions
  //for(i = 0; i < 4*NumHexon; i++) {
  for(i = 0; i < 2*NumHexon; i++) {
    ofs << HexonPoint[i](0) << " " 
	<< HexonPoint[i](1) << " "
	<< HexonPoint[i](2) << endl;
  }
  ofs << std::endl;
 


  // Cell section
  // ofs << "POINT_DATA    " << 4*NumHexon << endl;
  ofs << "POINT_DATA    " << 2*NumHexon << endl;
  // Initial stretch directions
  ofs << "VECTORS    MaxStretchDirections    double" << endl;
  // for(i = 0; i < 4*NumHexon; i++) {
  for(i = 0; i < 2*NumHexon; i++) {
    ofs << Stretch[i](0) << " " 
	<< Stretch[i](1) << " "
	<< Stretch[i](2) << endl;
  }
  ofs << endl;



  ofs << "VECTORS    MaxStrDirScaled    double" << endl;
  // for(i = 0; i < 4*NumHexon; i++) {
  for(i = 0; i < 2*NumHexon; i++) {
    ofs << StretchScaled[i](0) << " " 
	<< StretchScaled[i](1) << " "
	<< StretchScaled[i](2) << endl;
  }
  ofs << endl;



  ofs << "VECTORS    EtaDirections    double" << endl;
  for(i = 0; i < 2*NumHexon; i++) {
    ofs << EtaDir[i](0) << " " 
	<< EtaDir[i](1) << " "
	<< EtaDir[i](2) << endl;
  }
  ofs << endl;
  


  ofs.close();

 }



void printStretchEigenvector(const std::string OutFileNameD)
 {
  ofstream ofs(OutFileNameD.c_str());
  if (!ofs) {
    std::cout << "can not open output ("
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

	
  private:	

    const string _meshFileName;
    const string _outlineFileName;  
    const unsigned int _nPout; 

    const string _fileNameA;
    const string _fileNameB;
    const string _fileNameC;
    const string _fileNameD;
    const unsigned int _type;
    const unsigned int _MT;

    const Model::BodyContainer & _bdc;
    const vector<ScalarFieldNode<3>* > & _stretchNodes;
    const vector<ScalarFieldNode<3>* > & _directionNodes;  
    const vector<unsigned int> & _capsomersNum;

    const vector<double> _vectorAngleShift;
    const vector<double> _vectorReferenceAngle;
    const double _hexonShift;
    bool _finalConfiguration;

    vector<Vector3D> evec1, evec2, evec3;
  };
  
}; // namespace voom

#endif // __PrintingStretches_h__


