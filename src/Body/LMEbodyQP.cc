#include <string>
#include <fstream>
#include <blitz/array-impl.h>
#include "VoomMath.h"


#if defined(_OPENMP)
#include <omp.h>
#endif

#ifdef WITH_MPI
#include <mpi.h>
#endif

namespace voom
{

  template<class Material_t, class Shape_t >
  LMEbodyQP<Material_t,Shape_t>::LMEbodyQP(Material_t material,
					 const NodeContainer & nodes,
					 const NodeContainer & QP,
					 const std::vector<double> & node_volume,
					 const std::vector<double> & supp_size,
					 const double beta, const double searchR, const double tol, const int nItMax): 
    _nodes(nodes), _beta(beta), _searchR(searchR), _tol(tol), _nItMax(nItMax) 
  {

#ifdef WITH_MPI
    MPI_Comm_size( MPI_COMM_WORLD, &_nProcessors );
    MPI_Comm_rank( MPI_COMM_WORLD, &_processorRank );
#endif
        
    // initialize _nDOF
    _quadPoints.reserve(QP.size());

    bool FindNeigh = true;
    int i = 0, j = 0, Nfailed = 0, n = 0, a = 0;
    for(ConstNodeIterator nIt = QP.begin(); nIt != QP.end(); nIt++)
    {
      double weight = node_volume[i];
      Shape_t shapeN = LMEshape(beta, searchR*supp_size[i], tol, nItMax);
      LMEshape::CoordinateArray p((*nIt)->getPosition(0),(*nIt)->getPosition(1),(*nIt)->getPosition(2));
      // LMEshape::CoordinateArray p;
      cout << p << endl;
      // p(0) = 0.0;
      // p(1) = 0.0;
      // p(2) = 33.0;
 
      bool Store = shapeN.compute(p,nodes,FindNeigh, true);
      
      if (Store == true)
      {
	// Store linear system to compute actual nodal values
	const std::vector<int > & neighbours =  shapeN.neighb();
	n = neighbours.size();
	double *A;	
	A = (double*) malloc (n*n*sizeof(double));
	
	for(int a = 0; a < n; a++)
	{
	  const Vector3D & xa = nodes[neighbours[a]]->point();	
	  Shape_t shapeS = LMEshape(beta, searchR*supp_size[i], tol, nItMax);
	  shapeS.neighb() = neighbours;
	  bool Status = shapeS.compute(xa, nodes, false, false);
	  
	  if (Status == true)
	  {
	    const typename LMEshape::FunctionContainer &  N = shapeS.functions();
	    for(int j = 0; j < n; j++)
	    {
	      A[a+(j*n)] = N[j];
	    }
	  }
	  else 
	  {
	    free(A);
	    Nfailed++;
	    goto NewPoint;
	  }
	}
	
	_quadPoints.push_back(QuadPointStruct(weight, material, shapeN, A));
      }
      else
      {
	Nfailed++;
      }
      NewPoint:
      i++;
      if (i%10 == 0) cout << i << " nodes have been checked. " << endl;
    }
    cout << endl << "Failed points = " << Nfailed << " Percentage of failed points = " << double(Nfailed)/double(i) <<  endl;
  }

  //! Constructor for the quadrature points (same as nodes)
  template<class Material_t, class Shape_t >
  LMEbodyQP<Material_t,Shape_t>::QuadPointStruct::
  QuadPointStruct(double w, const Material_t & m, const Shape_t & s, double *_A): weight(w), material(m), A(_A)
  {
    shapeFunctions = s.functions();
    shapexDerivatives = s.xderivative();
    shapeyDerivatives = s.yderivative();
    shapezDerivatives = s.zderivative();
    neighbours = s.neighb();
  }

  
  
  //! Compute E0, E1, E2
  template<class Material_t, class Shape_t >
  void LMEbodyQP<Material_t,Shape_t>::compute( bool f0, bool f1, bool f2 )
  {
    // Initialize energy, forces and stiffness to be all zero
    if(f0) _energy = 0.0;
    if(f1) {
      for(ConstNodeIterator n = _nodes.begin(); n !=_nodes.end(); n++) 
	for(int i = 0; i < (*n)->dof(); i++)
	  (*n)->setForce(i,0.0);
    }
    
#ifdef _OPENMP	
#pragma omp parallel for			\
  schedule(static) default(shared) 
#endif	

    int iinod=0;
    for(ConstQuadPointIterator p = _quadPoints.begin(); p != _quadPoints.end(); p++) 
    {
      // Compute F
      Tensor3D F(0.0);
      const typename LMEshape::FunctionContainer &  xDN = p->shapexDerivatives;
      const typename LMEshape::FunctionContainer &  yDN = p->shapeyDerivatives;
      const typename LMEshape::FunctionContainer &  zDN = p->shapezDerivatives;
      const typename LMEshape::NodeNContainer & neighbours = p->neighbours;
      
      for(int a = 0; a < neighbours.size(); a++) 
      {
	const Vector3D & xa = _nodes[neighbours[a]]->point();
	for(int i = 0; i < 3; i++) 
	{
	  F(i,0) += xa(i)*xDN[a];
	  F(i,1) += xa(i)*yDN[a];
	  F(i,2) += xa(i)*zDN[a];
	}
      }
      
      Material_t material = p->material;
      // Send updated deformation gradient to material
      material.setDeformationGradient(F);
      // Compute strain energy and/or 1st PK stress
      material.updateState(f0, f1, f2); 
      
      // Compute energy
      if ( f0 )	_energy += material.energyDensity()*p->weight;
      
      // Compute forces
      if ( f1 ) {
	const Tensor3D & P = material.piolaStress();
	for (int a = 0; a < neighbours.size(); a++) {
	  for(int i = 0; i < 3; i++) {
	    _nodes[neighbours[a]]->addForce( i, (P(i,0)*xDN[a] + P(i,1)*yDN[a] + P(i,2)*zDN[a])*p->weight ); 
	  }
	}
      }
      
      // Compute stiffness matrix
      if( f2 ) {
	// Not implemented
      }
      
    } // end quad points loop
    
    return;
  }
  


  //! Compute invariants
  template<class Material_t, class Shape_t >
  void LMEbodyQP<Material_t, Shape_t>::cal_invariants( std::vector<double> & I1, std::vector<double> & I2, std::vector<double> & I3 )
  {
    // computes invariants of the deformation gradient at all the quadrature points (which for this class are assumed to be the same as nodes)
    /* if(_nodes.size()!=_quadPoints.size()){
      std::cout << "Number of nodes not equal to the number of quad points. Exiting"<<std::endl;
      exit(0);
    }
    */
#ifdef _OPENMP	
#pragma omp parallel for			\
  schedule(static) default(shared) 
#endif	

    bool FindNeigh = false;
    int inode=0;
    // cout << "Size of _quadPoints = " << _quadPoints.size() << endl;
    for(ConstQuadPointIterator p=_quadPoints.begin(); p!=_quadPoints.end(); p++){
      
      // Compute deformation gradient
      Tensor3D F(0.0), C(0.0), Csquare(0.0);
      const typename LMEshape::FunctionContainer &  xDN = p->shapexDerivatives;
      const typename LMEshape::FunctionContainer &  yDN = p->shapeyDerivatives;
      const typename LMEshape::FunctionContainer &  zDN = p->shapezDerivatives;
      const typename LMEshape::NodeNContainer & neighbours = p->neighbours;
      
      // Solve linear system to compute actual nodal values
      int  n = neighbours.size();
      int  nrhs = 3;
      int  lda  = n;
      int  ldb  = n;
      int  info;
      double *A = p->A, *Alap, *blap;
      Alap = (double*) malloc (n*n*sizeof(double));
      blap = (double*) malloc (n*3*sizeof(double));
      int ipiv[n];
    
      for(int a = 0; a < n; a++) {
	const Vector3D & xa = _nodes[neighbours[a]]->point();
	
	for(int i = 0; i < 3; i++) {
	  blap[a+i*n] = xa(i);
	}
	//	cout << "Point = " << xa << " " << b[a] << " " << b[a+n] << " " << b[a+2*n]  << endl;
      }

      for(int a = 0; a < n*n; a++) {
	Alap[a] = A[a];
      }

      // Solve for actual displacement
      dgesv_(&n, &nrhs, Alap, &lda, ipiv, blap, &ldb, &info);

 

      for(int a = 0; a < n; a++) {
	// const Vector3D & xa = _nodes[neighbours[a]]->point();
	/*	if ( fabs(blap[a]-xa(0)) > 0.1 )
	{
            cout << "Point = " << xa << " " << blap[a] << " " << blap[a+n] << " " << blap[a+2*n]  << endl;
	    }*/
	for(int i = 0; i < 3; i++) {
	  F(i,0) += blap[a+i*n]*xDN[a];
	  F(i,1) += blap[a+i*n]*yDN[a];
	  F(i,2) += blap[a+i*n]*zDN[a];

	  // F(i,0) += xa(i)*xDN[a];
	  // F(i,1) += xa(i)*yDN[a];
	  // F(i,2) += xa(i)*zDN[a];
	}
	// cout << xDN[a] << " " << yDN[a] << " " << zDN[a] << endl;
      }

      free(Alap);
      free(blap);
      // cout << F << endl;
      // assert(0);



      // compute the invariants
      double detF = F(0,0)*( F(1,1)*F(2,2) - F(1,2)*F(2,1) ) + F(0,1)*( F(2,0)*F(1,2) - F(1,0)*F(2,2)) + F(0,2)*( F(1,0)*F(2,1) - F(2,0)*F(1,1) );
      double detF_TwoThird = pow(detF, -2.0/3.0);
      for (int i = 0; i < 3; i++) {
	for (int j = 0; j < 3; j++) {
	  for (int k = 0; k < 3; k++) {
	    C(i,j) += F(k,i)*F(k,j);
	  }
	}
      }

      for (int i = 0; i < 3; i++) {
	for (int j = 0; j < 3; j++) {
	  for (int k = 0; k < 3; k++) {
	    Csquare(i,j) += C(i,k)*C(k,j);
	  }
	}
      }

      double trC = C(0,0) + C(1,1) + C(2,2);

      I1[inode] = trC*detF_TwoThird;
      I2[inode] = 0.5*(pow(trC,2.0) - Csquare(0,0) - Csquare(1,1) - Csquare(2,2))*pow(detF_TwoThird, 2.0);
      // I3[inode] = detF*detF; 
      I3[inode] = detF;

      if(I3[inode] < 0.0) std::cout<<"The jacobian is negative at ."<< inode << endl;

      inode++;
      
    } // end quad points loop

    return;
  }

 

  //! create input file used by Paraview, a 3D viewer
  template<class Material_t, class Shape_t >
  void LMEbodyQP<Material_t,Shape_t>::printParaview(const std::string name) const
  {
    std::string fileName = name + ".vtk";
    std::ofstream ofs(fileName.c_str());
    if (!ofs) {
      std::cout << "Can not open output ("
		<< fileName
		<< ") file." << std::endl;
      exit(0);
    }
    
    // Node Section
    ofs << "# vtk DataFile Version 2.0\n"
	<< "Test example" << std::endl
	<< "ASCII" << std::endl
	<< "DATASET UNSTRUCTURED_GRID" << std::endl
	<< "POINTS  " << _nodes.size() << "  double" << std::endl;
    
    // Output nodal reference postions
    for ( ConstNodeIterator pn = _nodes.begin(); pn!= _nodes.end(); pn ++) {
      const Vector3D & nodalPos =  (*pn)->position();
      ofs << std::setprecision(16) 
	  << nodalPos(0) << "  "
	  << nodalPos(1) << "  "
	  << nodalPos(2) << std::endl;
    }

    // Append the mesh data from file mesh.vtk into this file before appending displacements
    ofs.close();
    char temp[100];
    sprintf(temp,"cat mesh.vtk >> %s",fileName.c_str());
    system(temp);
    ofs.open(fileName.c_str(),std::ios::app); 
    
    // Output nodal final positions
    ofs << endl << "POINT_DATA " << _nodes.size() << endl
	<< "VECTORS displacements float" << endl;
    
    // Output nodal final postions
    for (ConstNodeIterator pn = _nodes.begin(); pn!= _nodes.end(); pn ++) {
      Vector3D nodalDisp;
      nodalDisp = (*pn)->point() - (*pn)->position();
      ofs << std::setprecision(16) 
	  << nodalDisp(0)
	  << '\t' <<nodalDisp(1)
	  << '\t' <<nodalDisp(2) << std::endl;
    }
    
    ofs.close();
    
    return;
  }
  


} // namespace voom
