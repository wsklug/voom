// -*- C++ -*-
//----------------------------------------------------------------------
//
//                    	    Ankush Aggarwal 
//                          William S.Klug
//                          Luigi Perotti
//                University of California Los Angeles
//                   (C) 2010 All Rights Reserved
//
//----------------------------------------------------------------------

#include <string>
#include <fstream>
#include <blitz/array-impl.h>

// #include "StretchHexonBody.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

#ifdef WITH_MPI
#include <mpe.h>
#endif

// define prototype of LAPACK functions
extern "C" void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda,
		       double *w, double *work, int *lwork, int *info);

namespace voom
{
  StretchHexonBody::StretchHexonBody(const ConnectivityContainer & connectivities,
				     const vector<DeformationNode<3>* > & DefNodes,
				     const StretchNodeContainer & StretchNodes,
				     const StretchNodeContainer & DirectionNodes,
				     const std::vector<double> & stretch_angle,
				     const double mu,
				     const double kS,
				     const int WcType,
				     double* WcConst,
				     const Quadrature<2> & Quad,
				     Shape<2> * shape,
				     const vector<unsigned int > & CapsomersNum,
				     ScalarFieldNode<3>* phiNode)
  {
    
#ifdef WITH_MPI
  MPI_Comm_size( MPI_COMM_WORLD, &_nProcessors );
  MPI_Comm_rank( MPI_COMM_WORLD, &_processorRank );
#endif

    _area = 0.0;
    _volume = 0.0;
    _energy = 0.0;
    
    // initialize _nDOF and find shell (position) nodes
    _dof = 0;
    _membraneNodes = DefNodes;
    _stretchNodes = StretchNodes;
    _directionNodes = DirectionNodes;
    // _phiNode = phiNode;
    _nodes.insert(_nodes.begin(), DefNodes.begin(), DefNodes.end() );
    _nodes.insert(_nodes.end(), StretchNodes.begin(), StretchNodes.end() );
    _nodes.insert(_nodes.end(), DirectionNodes.begin(), DirectionNodes.end() );
    // _nodes.push_back(phiNode);
    for(ConstNodeIterator n = _nodes.begin(); n != _nodes.end(); n++) {
      _dof+=(*n)->dof();
    }

    cout << "StretchHexonBody: initialized dof " << _dof << endl;
    
    

    // create elements
    int Cs = connectivities.size();
    _membraneElements.reserve(Cs);
    _elements.reserve(Cs);
    _active.reserve(Cs);
    if(connectivities.size() != stretch_angle.size())
      cout << "StretchHexonBody::Error, number of elements and stretch-angle are different." << endl;

    int ElCount=0;
    for(ConstConnectivityIterator c = connectivities.begin(); c != connectivities.end(); c++)
    {
      MembraneNodeContainer nds;   // build elemental node container
      C0MembraneStretch * MemEl;
      for(ElementConnectivity::const_iterator n = c->begin(); n != c->end(); n++)
      {
	nds.push_back(_membraneNodes[*n]);
      }
      // create element
      if (StretchNodes.size() == 1)
      {
	if (DirectionNodes.size() == 1)
	{
	  EvansElastic_Stretch * Mat = new  EvansElastic_Stretch(mu, kS, StretchNodes[0]->point(), stretch_angle[ElCount]+DirectionNodes[0]->point(),
							       WcType, WcConst);
	  MemEl = new C0MembraneStretch(nds, StretchNodes[0], DirectionNodes[0], phiNode, stretch_angle[ElCount], Quad, Mat, shape);
	}
	else
	{
	  int cap = CapsomersNum[ElCount];
	  EvansElastic_Stretch * Mat = new  EvansElastic_Stretch(mu, kS, StretchNodes[0]->point(), stretch_angle[ElCount]+DirectionNodes[cap]->point(),
							       WcType, WcConst);
	  MemEl = new C0MembraneStretch(nds, StretchNodes[0], DirectionNodes[cap], phiNode, stretch_angle[ElCount], Quad, Mat, shape);
	}
      }
      else
      {
	if (DirectionNodes.size() == 1)
	{
	  int cap = CapsomersNum[ElCount];
	  EvansElastic_Stretch * Mat = new  EvansElastic_Stretch(mu, kS, StretchNodes[cap]->point(), stretch_angle[ElCount]+DirectionNodes[0]->point(),
							       WcType, WcConst);
	  MemEl = new C0MembraneStretch(nds, StretchNodes[cap], DirectionNodes[0], phiNode, stretch_angle[ElCount], Quad, Mat, shape);
	}
	else
	{
	  int cap = CapsomersNum[ElCount];
	  EvansElastic_Stretch * Mat = new  EvansElastic_Stretch(mu, kS, StretchNodes[cap]->point(), stretch_angle[ElCount]+DirectionNodes[cap]->point(),
							       WcType, WcConst);
	  MemEl = new C0MembraneStretch(nds, StretchNodes[cap], DirectionNodes[cap], phiNode, stretch_angle[ElCount], Quad, Mat, shape);
	}
      }

      _membraneElements.push_back(MemEl);
      _elements.push_back(MemEl);
      _active.push_back(true);
      ElCount++;
    }

    cout << "StretchHexonBody: initialized elements " << ElCount << endl;
    // compute mechanics/geometry
    compute(false, false, false);

    // No constraints
    _constraints.clear();
    
  }



  //! compute
  void StretchHexonBody::compute( bool f0, bool f1, bool f2 )
  {
    int i = 0;
    int ElCount = _membraneElements.size();
    // Predictor/corrector approach for constraint
    for(ConstraintIterator c = _constraints.begin(); c != _constraints.end(); c++) {
      (*c)->predict();
    }

    // Initialize energy and forces
    if(f0) _energy = 0.0;

#ifdef _OPENMP	
#pragma omp parallel for 			\
  schedule(static) default(shared)		
#endif
    for(i = 0; i < ElCount; i++)
    {
       if( !_active[i] ) continue;
       _membraneElements[i]->compute(false,false,false);
    }

    // compute element areas and volumes
    double volume = 0.0;
    double area = 0.0;
	      
    for(i = 0; i < ElCount; i++)
    {
      if( !_active[i] ) continue;
      _volume += _membraneElements[i]->volume();
      _area += _membraneElements[i]->area();
    }

#ifdef _OPENMP	
#pragma omp parallel for 			\
  schedule(static) default(shared)		
#endif
    for(i = 0; i < ElCount; i++)
    {
      if( !_active[i] ) continue;
      _membraneElements[i]->compute( f0, f1, f2 );
    }

    if(f0) { 	
      for(i = 0; i < ElCount; i++)
      {
	if( !_active[i] ) continue;
	_energy += _membraneElements[i]->energy();
      }
    }

    // Predictor/corrector approach for constraint
    for(ConstraintIterator c = _constraints.begin(); c != _constraints.end(); c++)
    {
      (*c)->correct();
    }

    return;
  }
	

	
  //! create input file used by Paraview, a 3D viewer
  void StretchHexonBody::printParaview(const std::string name) const
  {
    int ElCount = _membraneElements.size(), e = 0, npts = 0;
    
    blitz::Array<double,1> energy(blitz::shape(ElCount));
    energy = 0.0;

    ConstMembraneElementIterator pe;
    MembraneElement::ConstQuadPointIterator p;
    for (pe = _membraneElements.begin(), e = 0; pe != _membraneElements.end(); pe++, e++) 
    {
      npts = 0;
      for(p = (*pe)->quadPoints().begin() ; p != (*pe)->quadPoints().end(); p++)
      {
	energy(e) += p->material->energyDensity();
	npts++;
      }
      energy(e) /= (double)( npts );
    }



    std::string fileName = name + ".vtk";
    std::ofstream ofs(fileName.c_str());

    if (!ofs) {
      cout << "Cannot open output (" << fileName << ") file." << endl;
      exit(0);
    }
    
    //  node Section
    ofs << "# vtk DataFile Version 2.0\n"
	<< "Test example" << endl
	<< "ASCII" << endl
	<< "DATASET POLYDATA" << endl
	<< "POINTS  " << _membraneNodes.size() << "  double" << endl;

    // output nodal postions
    ConstMembraneNodeIterator pn;
    for (pn = _membraneNodes.begin(); pn != _membraneNodes.end(); pn++)
    {
      const Vector3D & nodalPos = (*pn)->position();
      ofs << std::setprecision(16) 
	  << nodalPos(0) << "  "
	  << nodalPos(1) << "  "
	  << nodalPos(2) << endl;
    }

    // Element Section
    int nActiveMembranes = 0;
    for(e = 0; e < ElCount; e++) {
      if(_active[e]) nActiveMembranes++;
    }
    ofs << "POLYGONS  " << nActiveMembranes << "  " << 4*nActiveMembranes << endl;

    for(e = 0; e < ElCount; e++) {
      if(!_active[e])  continue;
      const MembraneNodeContainer & pnc =  _membraneElements[e]->nodes();
      ofs << 3 << "  "
	  << setw(10) << pnc[0] -> id()
	  << setw(10) << pnc[1] -> id()
	  << setw(10) << pnc[2] -> id()
	  << endl;
    }


     	
    // Output color for each element (corresponding to mean curvature, energy...)
    ofs << "CELL_DATA    " << nActiveMembranes << endl;
    // Output for strain energy
    ofs << "SCALARS    strainEnergy    double    1" << endl;
    ofs << "LOOKUP_TABLE default" << endl;
    for (e = 0; e < ElCount; e++)
    {
      if(!_active[e])  continue;
      ofs << energy(e) << endl;
    }
    ofs << endl;

    ofs << endl << "POINT_DATA " << _membraneNodes.size() << endl
	<< "VECTORS displacements double" << endl
        << "LOOKUP_TABLE displacements" << endl;

    // output nodal postions
    for (pn = _membraneNodes.begin(); pn != _membraneNodes.end(); pn ++)
    {
      Vector3D nodalDisp;
      nodalDisp = (*pn)->point() - (*pn)->position();
      ofs << setprecision(16) 
	  << nodalDisp(0)
	  << '\t' <<nodalDisp(1)
	  << '\t' <<nodalDisp(2) << endl;
    }
    
    ofs << endl << "VECTORS forces double" << endl
        << "LOOKUP_TABLE displacements" << endl; 

    // output nodal forces
   
    for (pn = _membraneNodes.begin(); pn != _membraneNodes.end(); pn ++)
    {
      const Vector3D & nodalForce = (*pn)->force();
      ofs << setprecision(16) 
	  << nodalForce(0)
	  << '\t' <<nodalForce(1)
	  << '\t' <<nodalForce(2) << endl;
    }
    ofs.close();

    return;
  }



  //! create output file used by Paraview to display CELL vectors
  void StretchHexonBody::printParaviewEigVec(const string name) const
  {
    string fileName = name + "-eigvec.vtk";
    ofstream ofs(fileName.c_str());
    if (!ofs) { 
      cout << "can not open output (" << fileName << ") file." << endl;
      exit(0);
    }
    
    int nActiveMembranes = 0, e = 0, ElCount = _membraneElements.size(), iNode = 0, npts = 0;
    for(e = 0; e < ElCount; e++) {
      if(_active[e]) nActiveMembranes++;
    }

    // Node Section
    // output nodal postions (positions of the element centroids)
    ofs << "# vtk DataFile Version 2.0\n"
	<< "Test example" << endl
	<< "ASCII" << endl
	<< "DATASET POLYDATA" << endl
	<< "POINTS  " << nActiveMembranes << "  double" << endl;
    
    for(e = 0; e < ElCount; e++)
    {
      if(!_active[e])  continue;
      const MembraneElement::NodeContainer & pnc = _membraneElements[e]->nodes();
      Vector3D center(0);
      for(iNode = 0; iNode < pnc.size(); iNode++)
        center += _membraneNodes[pnc[iNode]->id()]->position();

      center /= pnc.size();
      ofs << setprecision(16) 
	  << center(0) << "  "
	  << center(1) << "  "
	  << center(2) << endl;
    }

    // output nodal displacements (displacements of the element centroids)
    ofs << "POINT_DATA    " << nActiveMembranes << endl;
    ofs << "VECTORS displacements double" << endl;
    for(e = 0; e < ElCount; e++) 
    {
      if(!_active[e])  continue;
      const MembraneElement::NodeContainer & pnc = _membraneElements[e]->nodes();
      Vector3D center(0);
      for(iNode = 0; iNode < pnc.size(); iNode++)
        center += _membraneNodes[pnc[iNode]->id()]->point()-_membraneNodes[pnc[iNode]->id()]->position();

      center /= pnc.size();
      ofs << setprecision(16) 
	  << center(0) << "  "
	  << center(1) << "  "
	  << center(2) << endl;
    }
    
    vector<Vector3D> evec1, evec2, evec3;
    vector<double> eval1, eval2, eval3;

    // calculate the eigenvalues of B=F*F^T
    ConstMembraneElementIterator pe;
    MembraneElement::ConstQuadPointIterator p; 
    for (pe = _membraneElements.begin(), e = 0; pe != _membraneElements.end(); pe++, e++)
    {
      if(!_active[e])  continue;
      npts = 0;
      Tensor3D F;
      for(p = (*pe)->quadPoints().begin(); p != (*pe)->quadPoints().end(); p++)
      {
	F = p->material->DefGradient();
	npts++;
      }
      if( npts!=1 ) cerr << "printParaviewEigVec::Assumption of one GQ pt. not true" << endl;

      Tensor3D B(F*trans(F));

      // compute Eigenvalues and Eigenvectors by calling LAPACK library
      char jobz = 'V';
      char uplo = 'L';
      int  n    = 3;
      int  lda  = n;
      int  lwork = 3*n-1;
      int  info;
      double evalues[3];
      double work[lwork];
      
      // calling lapack function here to compute
      dsyev_(&jobz, &uplo, &n, B.data(),&lda, evalues, work, &lwork, &info);
      if (info != 0) {
	cout << "Something is wrong in DSYEV_" << endl;
	exit(0);
      }
      Vector3D vec(row(B,0));      
                         evec1.push_back(vec);
      vec=row(B,1);      evec2.push_back(vec);
      vec=row(B,2);      evec3.push_back(vec);
      eval1.push_back(evalues[0]);
      eval2.push_back(evalues[1]);
      eval3.push_back(evalues[2]);
    }
     	
    //  output vectors for each element (point in this file)
    ofs << "VECTORS    Beigvec1    double    " << endl;
    for (e = 0; e < nActiveMembranes; e++)
    {
      ofs << setprecision(16) 
	  << evec1[e](0)
	  << '\t' <<evec1[e](1)
	  << '\t' <<evec1[e](2) << endl;
    }
    ofs << endl;
    
    ofs << "VECTORS    Beigvec2    double    " << endl;
    for (e = 0; e < nActiveMembranes; e++)
    {
      ofs << setprecision(16) 
	  << evec2[e](0)
	  << '\t' <<evec2[e](1)
	  << '\t' <<evec2[e](2) << endl;
    }
    ofs << endl;
    
    ofs << "VECTORS    Beigvec3    double    " << endl;
    for (e = 0; e < nActiveMembranes; e++) {
      ofs << setprecision(16) 
	  << evec3[e](0)
	  << '\t' <<evec3[e](1)
	  << '\t' <<evec3[e](2) << endl;
    }
    ofs << endl;

    //  output eigen values for each element (point in this file)
    ofs << "SCALARS    B-Eval1    float    1" << endl;
    ofs << "LOOKUP_TABLE default" << endl;
    for (e = 0; e < nActiveMembranes; e++) 
      ofs << eval1[e] << endl; 
    
    ofs << "SCALARS    B-Eval2    float    1" << endl;
    ofs << "LOOKUP_TABLE default" << endl;
    for (e = 0; e < nActiveMembranes; e++) 
      ofs << eval2[e] << endl; 

    ofs << "SCALARS    B-Eval3    float    1" << endl;
    ofs << "LOOKUP_TABLE default" << endl;
    for (e = 0; e < nActiveMembranes; e++) 
      ofs << eval3[e] << endl; 

    // Calculate the eigenevectors and eigenvalues of cauchy stress and output that too
    evec1.clear(); evec2.clear(); evec3.clear();
    eval1.clear(); eval2.clear(); eval3.clear();

    pe = _membraneElements.begin();
    for (e = 0; pe != _membraneElements.end(); pe++, e++)
    {
      if(!_active[e])  continue;
      npts=0;
      Tensor3D sigma;
      for(p = (*pe)->quadPoints().begin(); p != (*pe)->quadPoints().end(); p++)
      {
	sigma = p->material->cauchyStress();
	npts++;
      }
      if(npts!=1) cerr << "printParaviewEigVec::Assumption of one GQ pt. not true" << endl;
      // compute Eigenvalues and Eigenvectors by calling LAPACK library
      char jobz = 'V';
      char uplo = 'L';
      int  n    = 3;
      int  lda  = n;
      int  lwork = 3*n-1;
      int  info;
      double evalues[3];
      double work[lwork];
      
      // calling lapack function here to compute
      dsyev_(&jobz, &uplo, &n, sigma.data(),&lda, evalues, work, &lwork, &info);
      if (info != 0) {
	std::cout << "Something is wrong in DSYEV_" << std::endl;
	exit(0);
      }
      double TOL=1e-12;
      Vector3D vec;

      if(std::fabs(evalues[0])<TOL){
      vec=row(sigma,0);      evec1.push_back(vec);
      eval1.push_back(evalues[0]);
      vec=row(sigma,1);      evec2.push_back(vec);
      eval2.push_back(evalues[1]);
      vec=row(sigma,2);      evec3.push_back(vec);
      eval3.push_back(evalues[2]);}

      else if(std::fabs(evalues[1])<TOL) {
      vec=row(sigma,1);      evec1.push_back(vec);
      eval1.push_back(evalues[1]);
      vec=row(sigma,0);      evec2.push_back(vec);
      eval2.push_back(evalues[0]);
      vec=row(sigma,2);      evec3.push_back(vec);
      eval3.push_back(evalues[2]);}

      else {
      vec=row(sigma,2);      evec1.push_back(vec);
      eval1.push_back(evalues[2]);
      vec=row(sigma,0);      evec2.push_back(vec);

      eval2.push_back(evalues[0]);
      vec=row(sigma,1);      evec3.push_back(vec);
      eval3.push_back(evalues[1]);}
    }
     	
    //  output vectors for each element (point in this file)
    ofs << "VECTORS    cauchy-eigvec1    double    " << std::endl;
    for ( int e = 0; e<nActiveMembranes; e++) {
      ofs << setprecision(16) 
	  << evec1[e](0)
	  << '\t' <<evec1[e](1)
	  << '\t' <<evec1[e](2) << endl;
    }
    ofs << std::endl;
    
    ofs << "VECTORS    cauchy-eigvec2    double    " << endl;
    for ( int e = 0; e<nActiveMembranes; e++) {
      ofs << setprecision(16) 
	  << evec2[e](0)
	  << '\t' <<evec2[e](1)
	  << '\t' <<evec2[e](2) << endl;
    }
    ofs << std::endl;
    
    ofs << "VECTORS    cauchy-eigvec3    double    " << endl;
    for ( int e = 0; e<nActiveMembranes; e++) {
      ofs << setprecision(16) 
	  << evec3[e](0)
	  << '\t' <<evec3[e](1)
	  << '\t' <<evec3[e](2) << endl;
    }
    ofs << std::endl;
    //  output eigen values for each element (point in this file)
    ofs << "SCALARS    cauchy-Eval1    float    1" << endl;
    ofs << "LOOKUP_TABLE default" << endl;
    for ( int e = 0; e<nActiveMembranes; e++) 
      ofs << eval1[e] << endl; 
    
    ofs << "SCALARS    cauchy-Eval2    float    1" << endl;
    ofs << "LOOKUP_TABLE default" << endl;
    for ( int e = 0; e<nActiveMembranes; e++) 
      ofs << eval2[e] << endl; 

    ofs << "SCALARS    cauchy-Eval3    float    1" << endl;
    ofs << "LOOKUP_TABLE default" << std::endl;
    for ( int e = 0; e<nActiveMembranes; e++) 
      ofs << eval3[e] << endl; 
    ofs.close();

    return;
  }



  //! create output file used by Paraview to display CELL vectors
  void StretchHexonBody::printParaview2(const string name) const
  {
    string fileName = name + "-eigvec.vtk";
    ofstream ofs(fileName.c_str());
    if (!ofs) {
      cout << "can not open output (" << fileName << ") file." << endl;
      exit(0);
    }
    
    int nActiveMembranes = 0, e = 0, ElCount = _membraneElements.size(), npts = 0;
    for(e = 0; e < ElCount; e++) 
    {
      if(_active[e]) nActiveMembranes++;
    }

    // Node Section 
    // output nodal postions (positions of the element centroids)
    ofs << "# vtk DataFile Version 2.0\n"
	<< "Test example" << endl
	<< "ASCII" << endl
	<< "DATASET POLYDATA" << endl
	<< "POINTS  " << _membraneNodes.size() << "  double" << endl;
    
    ConstMembraneNodeIterator pn;
    for (pn = _membraneNodes.begin(); pn != _membraneNodes.end(); pn++)
    {
      const Vector3D & nodalPos =  (*pn)->position();
      ofs << setprecision(16) 
	  << nodalPos(0) << "  "
	  << nodalPos(1) << "  "
	  << nodalPos(2) << endl;
    }

    // Element Section
    ofs << "POLYGONS  " << nActiveMembranes << "  " << 4*nActiveMembranes << endl;
    for(e = 0; e < ElCount; e++)
    {
      if(!_active[e]) continue;
      const MembraneElement::NodeContainer & pnc = _membraneElements[e]->nodes();
      ofs << 3 << "  "
	  << setw(10) << pnc[0] -> id()
	  << setw(10) << pnc[1] -> id()
	  << setw(10) << pnc[2] -> id()
	  << endl;
    }

    // Output nodal displacements (displacements of the element centroids)
    ofs << "POINT_DATA    " << _membraneNodes.size() << endl;
    ofs << "VECTORS displacements double" << endl;
    
    for (pn = _membraneNodes.begin(); pn != _membraneNodes.end(); pn++)
    {
      Vector3D nodalDisp;
      nodalDisp = (*pn)->point() - (*pn)->position();
      ofs << setprecision(16) 
	  << nodalDisp(0)
	  << '\t' <<nodalDisp(1)
	  << '\t' <<nodalDisp(2) << endl;
    }
    
    vector<Vector3D> evec1, evec2, evec3;
    vector<double> eval1, eval2, eval3;

    //calculate the eigenvalues of B=F*F^T
    ConstMembraneElementIterator pe;
    MembraneElement::ConstQuadPointIterator p;
    for (pe = _membraneElements.begin(), e = 0; pe != _membraneElements.end(); pe++, e++)
    {
      if(!_active[e])  continue;
      npts = 0;
      Tensor3D F;
      for(p = (*pe)->quadPoints().begin(); p != (*pe)->quadPoints().end(); p++)
      {
	F = p->material->DefGradient();
	npts++;
      }
      if(npts>1) cerr << "printParaview2::Assumption of one GQ pt. not true" << endl;
      Tensor3D B(F*trans(F));
      // compute Eigenvalues and Eigenvectors by calling LAPACK library
      char jobz = 'V';
      char uplo = 'L';
      int  n    = 3;
      int  lda  = n;
      int  lwork = 3*n-1;
      int  info;
      double evalues[3];
      double work[lwork];
      
      // calling lapack function here to compute
      dsyev_(&jobz, &uplo, &n, B.data(),&lda, evalues, work, &lwork, &info);
      if (info != 0) {
	cout << "Something is wrong in DSYEV_" << endl;
	exit(0);
      }
      Vector3D vec(row(B,0));      
                         evec1.push_back(vec);
      vec=row(B,1);      evec2.push_back(vec);
      vec=row(B,2);      evec3.push_back(vec);
      eval1.push_back(evalues[0]);
      eval2.push_back(evalues[1]);
      eval3.push_back(evalues[2]);
    }
     	
    //  output vectors for each element (point in this file)
    ofs << "CELL_DATA    " << nActiveMembranes << endl;
    ofs << "VECTORS    Beigvec1    double    " << endl;
    for (e = 0; e < nActiveMembranes; e++) {
      ofs << setprecision(16) 
	  << evec1[e](0)
	  << '\t' <<evec1[e](1)
	  << '\t' <<evec1[e](2) << endl;
    }
    ofs << endl;
    
    ofs << "VECTORS    Beigvec2    double    " << endl;
    for (e = 0; e < nActiveMembranes; e++) {
      ofs << setprecision(16) 
	  << evec2[e](0)
	  << '\t' <<evec2[e](1)
	  << '\t' <<evec2[e](2) << endl;
    }
    ofs << endl;
    
    ofs << "VECTORS    Beigvec3    double    " << endl;
    for (e = 0; e < nActiveMembranes; e++) {
      ofs << setprecision(16) 
	  << evec3[e](0)
	  << '\t' <<evec3[e](1)
	  << '\t' <<evec3[e](2) << endl;
    }
    ofs << endl;
    //  output eigen values for each element (point in this file)
    ofs << "SCALARS    B-Eval1    float    1" << endl;
    ofs << "LOOKUP_TABLE default" << endl;
    for (e = 0; e < nActiveMembranes; e++) 
      ofs << eval1[e] << endl; 
    
    ofs << "SCALARS    B-Eval2    float    1" << endl;
    ofs << "LOOKUP_TABLE default" << std::endl;
    for (e = 0; e < nActiveMembranes; e++) 
      ofs << eval2[e] << endl; 

    ofs << "SCALARS    B-Eval3    float    1" << endl;
    ofs << "LOOKUP_TABLE default" << endl;
    for (e = 0; e < nActiveMembranes; e++) 
      ofs << eval3[e] << endl; 

    // Calculate the eigenevectors and eigenvalues of cauchy stress and output that
    evec1.clear(); evec2.clear(); evec3.clear();
    eval1.clear(); eval2.clear(); eval3.clear();

    for (pe = _membraneElements.begin(), e = 0; pe != _membraneElements.end(); pe++, e++)
    {
      if(!_active[e])  continue;
      npts = 0;
      Tensor3D sigma;
      for(p = (*pe)->quadPoints().begin(); p != (*pe)->quadPoints().end(); p++) {
	sigma = p->material->cauchyStress();
	npts++;
      }
      if(npts != 1) cerr << "printParaview2::Assumption of one GQ pt. not true" << endl;
      // compute Eigenvalues and Eigenvectors by calling LAPACK library
      char jobz = 'V';
      char uplo = 'L';
      int  n    = 3;
      int  lda  = n;
      int  lwork = 3*n-1;
      int  info;
      double evalues[3];
      double work[lwork];
      
      // calling lapack function here to compute
      dsyev_(&jobz, &uplo, &n, sigma.data(),&lda, evalues, work, &lwork, &info);
      if (info != 0) {
	std::cout << "Something is wrong in DSYEV_" << std::endl;
	exit(0);
      }
      double TOL=1e-12;
      Vector3D vec;

      if( fabs(evalues[0]) < TOL ){
      vec=row(sigma,0);      evec1.push_back(vec);
      eval1.push_back(evalues[0]);
      vec=row(sigma,1);      evec2.push_back(vec);
      eval2.push_back(evalues[1]);
      vec=row(sigma,2);      evec3.push_back(vec);
      eval3.push_back(evalues[2]);}

      else if( fabs(evalues[1]) < TOL ) {
      vec=row(sigma,1);      evec1.push_back(vec);
      eval1.push_back(evalues[1]);
      vec=row(sigma,0);      evec2.push_back(vec);
      eval2.push_back(evalues[0]);
      vec=row(sigma,2);      evec3.push_back(vec);
      eval3.push_back(evalues[2]);}

      else {
      vec=row(sigma,2);      evec1.push_back(vec);
      eval1.push_back(evalues[2]);
      vec=row(sigma,0);      evec2.push_back(vec);
      eval2.push_back(evalues[0]);
      vec=row(sigma,1);      evec3.push_back(vec);
      eval3.push_back(evalues[1]);}
    }
     	
    //  output vectors for each element (point in this file)
    ofs << "VECTORS    cauchy-eigvec1    double    " << endl;
    for (e = 0; e < nActiveMembranes; e++) {
      ofs << setprecision(16) 
	  << evec1[e](0)
	  << '\t' <<evec1[e](1)
	  << '\t' <<evec1[e](2) << endl;
    }
    ofs << endl;
    
    ofs << "VECTORS    cauchy-eigvec2    double    " << endl;
    for (e = 0; e < nActiveMembranes; e++) {
      ofs << setprecision(16) 
	  << evec2[e](0)
	  << '\t' <<evec2[e](1)
	  << '\t' <<evec2[e](2) << endl;
    }
    ofs << endl;
    
    ofs << "VECTORS    cauchy-eigvec3    double    " << endl;
    for (e = 0; e < nActiveMembranes; e++) {
      ofs << setprecision(16) 
	  << evec3[e](0)
	  << '\t' <<evec3[e](1)
	  << '\t' <<evec3[e](2) << endl;
    }
    ofs << endl;
    //  output eigen values for each element (point in this file)
    ofs << "SCALARS    cauchy-Eval1    float    1" << endl;
    ofs << "LOOKUP_TABLE default" << endl;
    for (e = 0; e < nActiveMembranes; e++) 
      ofs << eval1[e] << endl; 
    
    ofs << "SCALARS    cauchy-Eval2    float    1" << endl;
    ofs << "LOOKUP_TABLE default" << endl;
    for (e = 0; e < nActiveMembranes; e++) 
      ofs << eval2[e] << endl; 

    ofs << "SCALARS    cauchy-Eval3    float    1" << endl;
    ofs << "LOOKUP_TABLE default" << endl;
    for (e = 0; e < nActiveMembranes; e++) 
      ofs << eval3[e] << endl; 
    ofs.close();

    return;
  }
} // namespace voom
