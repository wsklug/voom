//----------------------------------------------------------------------
//
//                    William S. Klug, Feng Feng
//                University of California Los Angeles
//                 (C) 2004-2007 All Rights Reserved
//
//----------------------------------------------------------------------



#include <string>
#include <fstream>
#include <ctime>
#include <blitz/array-impl.h>
#include "HalfEdgeMesh.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

#ifdef WITH_MPI
#include <mpe.h>
#endif

namespace voom
{
  using std::endl;
  using std::cout;

  template<class Material_t>
  void GLBody<Material_t>::initializeBody
  (Material_t material,
   ConnectivityContainer & connectivities,
   const NodeContainer & nodes,
   const unsigned quadOrder,
   const int numberOfBoundaries,
   const double pressure,
   const double tension,
   const double penaltyVolume,
   const double penaltyArea,
   GlobalConstraint volumeConstraint,
   GlobalConstraint areaConstraint ) {

#ifdef WITH_MPI
  MPI_Comm_size( MPI_COMM_WORLD, &_nProcessors );
  MPI_Comm_rank( MPI_COMM_WORLD, &_processorRank );
#endif

    _penaltyArea = penaltyArea;
    _penaltyVolume = penaltyVolume;
    _areaConstraint = areaConstraint;
    _volumeConstraint = volumeConstraint;
    _volume = 0.0;
    _area = 0.0;

	  
    // initialize _nDOF and find shell (position) nodes

    _dof = 0;
    _nodes = nodes;
    for(ConstNodeIterator n=_nodes.begin(); n!=_nodes.end(); n++) {
      _dof+=(*n)->dof();
      FeNode_t * sn = dynamic_cast<FeNode_t*>(*n);
      if( sn ) _shellNodes.push_back(sn); 
    }
        
    NodeBase::DofIndexMap idx(1);
    idx[0]=-1;
    _pressureNode = new MultiplierNode(_nodes.size(), idx, pressure );
    _fixedPressure = pressure;

    _tensionNode = new MultiplierNode(_nodes.size(), idx, tension );
    _fixedTension = tension;
    
    std::cout << "GLBody::initializeBody(): dof = "<< _dof 
	      << std::endl;

    //------------------------------------------------------------------------
    // build elements using simple data structure
    std::cout << "Building elements..." << std::endl;
    clock_t t1=clock();

    HalfEdgeMesh mesh(connectivities, _shellNodes.size());
       
    for(int f=0; f<mesh.faces.size(); f++) {
      Face * F = mesh.faces[f];
      // get valences of corners and find one-ring of
      // next-nearest-neighbor vertices
      typename GLElement<Material_t>::NodeContainer nds;
      CornerValences v;
      //check if any of the vertices is on the boundary
      if( F->halfEdges[0]->vertex->boundary ) continue;
      if( F->halfEdges[1]->vertex->boundary ) continue;
      if( F->halfEdges[2]->vertex->boundary ) continue;

      for(int a=0; a<3; a++) {	
	// add corner vertices
	HalfEdge * H = F->halfEdges[a];
	nds.push_back( _shellNodes[H->vertex->id] );
	// valence is the number of incident halfedges
	v(a) = H->vertex->halfEdges.size();
      }
      // walk around the one-ring adding vertices
      //
      // Start by walking CCW through the vertices incident to the
      // second corner vertex.  Stop at first vertex incident to the
      // third corner vertex.
      int v1=3;
      for(HalfEdge * H = F->halfEdges[1]->opposite->next; 
	  H->vertex->id != F->halfEdges[2]->opposite->next->vertex->id; 
	  H = H->next->opposite->next  ) {
	nds.push_back( _shellNodes[H->vertex->id] );
	v1++;
      }
      // Now continue CCW through the vertices incident to the third
      // corner vertex. Stop at first vertex incident to the first
      // corner vertex.
      int v2=3;
      for(HalfEdge * H = F->halfEdges[2]->opposite->next; 
	  H->vertex->id != F->halfEdges[0]->opposite->next->vertex->id; 
	  H = H->next->opposite->next  ) {
	nds.push_back( _shellNodes[H->vertex->id] );
	v2++;
      }
      // Now continue CCW through the vertices incident to the first
      // corner vertex. Stop at first vertex incident to the second
      // corner vertex.
      int v0=3;
      for(HalfEdge * H = F->halfEdges[0]->opposite->next; 
	  H->vertex->id != F->halfEdges[1]->opposite->next->vertex->id; 
	  H = H->next->opposite->next  ) {
	nds.push_back( _shellNodes[H->vertex->id] );
	v0++;
      }
      assert(v0==v(0));
      assert(v1==v(1));
      assert(v2==v(2));
      const unsigned npn = v(0) + v(1) + v(2) - 6;		
      FeElement_t * elem 
	= new FeElement_t( TriangleQuadrature(quadOrder), 
			   material,
			   nds,
			   v,
			   _pressureNode,
			   _tensionNode,
			   _volumeConstraint,
			   _areaConstraint
			   );
      _shells.push_back(elem);
    }
    _elements.reserve(_shells.size());
    for(ConstFeElementIterator e=_shells.begin(); e!=_shells.end(); e++) {
      _elements.push_back(*e);
    }

    clock_t t2=clock();
    std::cout << "Done building elements.  Elapsed time: "
	      << (double)(t2-t1)/CLOCKS_PER_SEC
	      << std::endl;

    // compute mechanics
    compute(false, false, false);
    _prescribedVolume = _volume;
    _prescribedArea = _area;

  }

  //! compute
  template< class Material_t >
  void GLBody<Material_t>::compute( bool f0, bool f1, bool f2 )
  {
    
    int eBegin=0, eEnd=_elements.size();
    int sBegin=0, sEnd=_shells.size();
#ifdef WITH_MPI
    MPE_Decomp1d( _elements.size(), _nProcessors, _processorRank, 
		  &eBegin, &eEnd );
    eBegin--;
    MPE_Decomp1d( _shells.size(), _nProcessors, _processorRank, 
		  &sBegin, &sEnd );
    sBegin--;
//     std::cout << _processorRank << ": [" << eBegin << "," << eEnd << ")"
// 	      << " of " << _elements.size() << " elements." << std::endl;    
//     MPI_Finalize();
//     exit(0);
#endif

    // Predictor/corrector approach for contact
    for(ContactIterator c=_contacts.begin(); c!=_contacts.end(); c++) {
      (*c)->predict();
    }


    //
    // compute element areas and volumes

#ifdef _OPENMP	
#pragma omp parallel for 			\
  schedule(static) default(shared)		
#endif	      
     for(int si=sBegin; si<sEnd; si++) {
       FeElement_t* s=_shells[si];
       s->compute( false,false,false );
//       volume = volume + s->volume();     
//       area = area + s->area();
     }
//     _volume=volume;
//     _area=area;


      double volume = 0.0;
      double area = 0.0;
      
// For some reason, icc-9.1 doesn't like the reduction clauses here.
// I think it is a bug in icc. -WSK

// #ifdef _OPENMP	
// #pragma omp parallel for 			\
//   schedule(static) default(shared)		\
//   reduction(+:V) reduction(+:A) 
// #endif	      
      for(int si=sBegin; si<sEnd; si++) {
	FeElement_t* s=_shells[si];
	volume += s->volume();     
	area += s->area();	
      }

      _volume=volume;
      _area=area;


#ifdef WITH_MPI
    double myVolume=_volume;
    MPI_Allreduce(&myVolume, &_volume, 1, MPI_DOUBLE, 
		  MPI_SUM, MPI_COMM_WORLD);
    double myArea=_area;
    MPI_Allreduce(&myArea, &_area, 1, MPI_DOUBLE, 
		  MPI_SUM, MPI_COMM_WORLD);
#endif 
      //std::cout << "Body volume = " << _volume << std::endl;

      // if total volume is less than 0, the element faces may be
      // incorrectly oriented.
//       if (_volume <= 0.0){ 
// 	std::cout << "The current volume is less than Zero and it = "
// 		  << _volume
// 		  << std::endl
// 		  << "Element faces may be incorrectly oriented or mesh may be inverted." 
// 		  << std::endl;
// 	print("NonPositiveVolume");

	//exit(0);
//       }


    double dV = _volume - _prescribedVolume;	
    double dA = _area - _prescribedArea;
    
    if(f0) {
      _energy = 0.0;
      
      if( _volumeConstraint == penalty || _volumeConstraint == augmented ) {
	_energy += 0.5 * _penaltyVolume * sqr(dV/_prescribedVolume);
      }
      if( _volumeConstraint == multiplier || _volumeConstraint == augmented ) {
	_energy += - _fixedPressure * dV;
      }
      

      if( _areaConstraint == penalty || _areaConstraint == augmented ) {
	_energy += 0.5 * _penaltyArea * sqr(dA/_prescribedArea);
      }
      if( _areaConstraint == multiplier || _areaConstraint == augmented ) {
	_energy += _fixedTension * dA;
      }

    }
    
    double pressure = 0.0;
    double tension = 0.0;
    if( _volumeConstraint == multiplier || _volumeConstraint ==augmented ) {
      pressure = _fixedPressure;
    }

    if( _volumeConstraint == penalty ||_volumeConstraint == augmented ) {
      pressure += 
	-_penaltyVolume * (_volume/_prescribedVolume-1.0)/_prescribedVolume;
    }
    _pressureNode->setPoint(pressure);

    if( _areaConstraint == multiplier || _areaConstraint == augmented ) {
      tension = _fixedTension;
    }
    if( _areaConstraint == penalty || _areaConstraint == augmented ) {
      tension += 
	_penaltyArea * (_area/_prescribedArea-1.0)/_prescribedArea;
    }
    _tensionNode->setPoint(tension);


    // Model now initializes forces and stiffness
//     if(f1) {
//       for(NodeIterator n=_nodes.begin(); n!=_nodes.end(); n++) 
// 	for(int i=0; i<(*n)->dof(); i++)
// 	  (*n)->setForce(i,0.0);
//     }
    //
    // Need to zero out stiffness too!!!!!!!!!!
    //
    
    //
    // compute energy, forces and stiffness matrix in each element
    // loop through all elements
#ifdef _OPENMP	
#pragma omp parallel for 		\
  schedule(static) default(shared) 
#endif	
    for(int ei=eBegin; ei<eEnd; ei++) {
//       Element* e=_elements[ei];
      _elements[ei]->compute( f0, f1, f2 );
    }

    if(f0) { 
// #ifdef _OPENMP	
// #pragma omp parallel for
// #endif	
      for(int ei=eBegin; ei<eEnd; ei++) {
	Element* e=_elements[ei];
	_energy += e->energy();
      }
// #ifdef WITH_MPI
//       double myEnergy=_energy;
//       MPI_Allreduce(&myEnergy, &_energy, 1, MPI_DOUBLE, 
// 		    MPI_SUM, MPI_COMM_WORLD);
// #endif 
    }


// #ifdef WITH_MPI
//     if(f1) {
//       for(NodeIterator n=_nodes.begin(); n!=_nodes.end(); n++) 
// 	for(int i=0; i<(*n)->dof(); i++) {
// 	  double force = 0.0;
// 	  double myForce = (*n)->getForce(i); 
// 	  MPI_Allreduce(&myForce, &force, 1, MPI_DOUBLE, 
// 			MPI_SUM, MPI_COMM_WORLD);
// 	  (*n)->setForce(i,force);
// 	}
//     }
// #endif

    // Predictor/corrector approach for contact
    for(ContactIterator c=_contacts.begin(); c!=_contacts.end(); c++) {
      (*c)->correct();
    }
    return;
  }

  //! create input file used by Paraview, a 3D viewer
  template < class Material_t >
  void GLBody< Material_t >::printParaview(const std::string name) const
  {

    int numShells = _shells.size();
    
    blitz::Array<double,1> energy(blitz::shape(numShells));
    blitz::Array<double,1> bendingEnergy(blitz::shape(numShells));
    blitz::Array<double,1> inplaneEnergy(blitz::shape(numShells));
    blitz::Array<double,1> curvature(blitz::shape(numShells));

    blitz::Array<Tensor3D,1> strain(numShells);
    blitz::Array<Tensor3D,1> stress(numShells);

    energy = 0.0;
    curvature = 0.0;
    bendingEnergy = 0.0;
    inplaneEnergy = 0.0;
    strain = Tensor3D(0.0);
    stress = Tensor3D(0.0);

#ifdef WITH_MPI
    int eBegin=0, eEnd=0;
    MPE_Decomp1d( _shells.size(), _nProcessors, _processorRank, 
		  &eBegin, &eEnd );
    eBegin--;
#endif

#ifdef WITH_MPI
    for(int e=eBegin; e<eEnd && e<_shells.size(); e++) {
	const FeElement_t*const*const pe=&(_shells[e]);
#else
    ConstFeElementIterator pe = _shells.begin();
    for ( int e = 0; pe != _shells.end(); pe++, e++) {
#endif
      int npts=0;
      typename FeElement_t::ConstQuadPointIterator 
	p = (*pe)->quadraturePoints().begin();
      for( ; p != (*pe)->quadraturePoints().end(); p++){
	energy(e) += p->material.energyDensity();
	bendingEnergy(e) += p->material.bendingEnergy();
	inplaneEnergy(e) += p->material.stretchingEnergy();
	curvature(e) += p->material.meanCurvature();

	// compute global Cartesian strain and stress
	Tensor3D E(0.0),S(0.0);
	const Tensor2D & Ecov = p->material.strain();
	const Tensor2D & Scon = p->material.stress();

	for(int alpha=0; alpha<2; alpha++) {
	  const Vector3D & base_alpha = 
	    p->material.refShellGeometry().a()[alpha];
	  const Vector3D & dual_alpha = 
	    p->material.refShellGeometry().aDual()[alpha];
	  for(int beta=0; beta<2; beta++) {
	    const Vector3D & base_beta = 
	      p->material.refShellGeometry().a()[beta];
	    const Vector3D & dual_beta = 
	      p->material.refShellGeometry().aDual()[beta];
	    for(int I=0; I<3; I++) {
	      for(int J=0; J<3; J++) {
		E(I,J) += Ecov(alpha,beta)*dual_alpha(I)*dual_beta(J);
		S(I,J) += Scov(alpha,beta)*base_alpha(I)*base_beta(J);
	      }
	    }
	  }
	}
	strain(e) += E;
	stress(e) += S;

	npts++;
      }
      energy(e) /= (double)( npts );
      bendingEnergy(e) /= (double)( npts );
      inplaneEnergy(e) /= (double)( npts );
      curvature(e) /= (double)( npts );
      strain(e) /= (double)( npts );
      stress(e) /= (double)( npts );
    }
//     char pid[30]; 
//     sprintf(pid, "%d", _processorRank);
//     std::string dataName = name + ".printData." + pid;
//     std::ofstream dataStream(dataName.c_str());
//     dataStream << energy << std::endl
// 	       << curvature << std::endl;
//     dataStream.close();

#ifdef WITH_MPI
    blitz::Array<double,1> globalsum(blitz::shape(numShells));
    globalsum = 0.0;
    MPI_Reduce(energy.data(), globalsum.data(), energy.size(), 
	       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    blitz::cycleArrays(energy,globalsum);
    MPI_Reduce(bendingEnergy.data(), globalsum.data(), bendingEnergy.size(), 
	       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    blitz::cycleArrays(bendingEnergy,globalsum);
    MPI_Reduce(inplaneEnergy.data(), globalsum.data(), inplaneEnergy.size(), 
	       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    blitz::cycleArrays(inplaneEnergy,globalsum);
    MPI_Reduce(curvature.data(), globalsum.data(), curvature.size(), 
	       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    blitz::cycleArrays(curvature,globalsum);
#endif

#ifdef WITH_MPI
    if(_processorRank!=0) {
      return;
    }
#endif

//     dataName = name + ".printData.sum";
//     dataStream.open(dataName.c_str());
//     dataStream << energy << std::endl
// 	       << curvature << std::endl;
//     dataStream.close();

    std::string fileName = name + ".vtk";
    std::ofstream ofs(fileName.c_str());
//     ofs.open(fileName.c_str(), std::ios::out);
    if (!ofs) {
      std::cout << "can not open output ("
		<< fileName
		<< ") file." << std::endl;
      exit(0);
    }
    
    ////////////////////////////////////////////////////////////////////
    //
    //    Node Section
    //
    ofs << "# vtk DataFile Version 2.0\n"
	<< "Test example" << std::endl
	<< "ASCII" << std::endl
	<< "DATASET POLYDATA" << std::endl
	<< "POINTS  " << _shellNodes.size() << "  double" << std::endl;
    
    //
    // output nodal postions
    ConstFeNodeIterator pn = _shellNodes.begin();
    for ( ; pn!= _shellNodes.end(); pn ++) {
      const Vector3D & nodalPos =  (*pn)->position();
      ofs << std::setprecision(16) 
	  << nodalPos(0) << "  "
	  << nodalPos(1) << "  "
	  << nodalPos(2) << std::endl;
    }

    /////////////////////////////////////////////////////////////////////
    //
    //    Element Section
    //
    ofs << "POLYGONS  " << _shells.size() << "  "
	<< 4*_shells.size() << std::endl;
    for (ConstFeElementIterator pe = _shells.begin(); pe != _shells.end(); pe++) {
      const typename FeElement_t::NodeContainer & pnc = (*pe)->nodes();
      ofs << 3 << "  "
	  << std::setw(10) << pnc[0] -> id()
	  << std::setw(10) << pnc[1] -> id()
	  << std::setw(10) << pnc[2] -> id()
	  << std::endl;
    }
     	
    //////////////////////////////////////////////////////////////////////////
    //
    //  output color for each element ( corresponding to mean
    //  curvature, energy...)
    //
    ofs << "CELL_DATA    " << _shells.size() << std::endl;
    //
    // output for strain energy
    ofs << "SCALARS    strainEnergy    double    1" << std::endl;
    ofs << "LOOKUP_TABLE default" << std::endl;
    for ( int e = 0; e<numShells; e++) 
      ofs << energy(e) << std::endl;
    ofs << std::endl;

    //
    // output for strain energy
    ofs << "SCALARS    bendingEnergy    double    1" << std::endl;
    ofs << "LOOKUP_TABLE default" << std::endl;
    for ( int e = 0; e<numShells; e++) 
      ofs << bendingEnergy(e) << std::endl;
    ofs << std::endl;

    //
    // output for strain energy
    ofs << "SCALARS    inplaneEnergy    double    1" << std::endl;
    ofs << "LOOKUP_TABLE default" << std::endl;
    for ( int e = 0; e<numShells; e++) 
      ofs << inplaneEnergy(e) << std::endl;
    ofs << std::endl;

    //
    // output for mean curvature
    // since one color is mapping into one element, for
    // several guass points integration, we compute their average value
    ofs << "SCALARS    meanCurvature    double    1" << std::endl;
    ofs << "LOOKUP_TABLE default    " << std::endl;		
    for ( int e = 0; e<numShells; e++) 
      ofs << curvature(e) << std::endl;
    ofs << std::endl;

    //
    // output for strain 
    ofs << "TENSORS    strain    double" << std::endl;
    ofs << "LOOKUP_TABLE default" << std::endl;
    for ( int e = 0; e<numShells; e++) 
      for( int I=0; I<3; I++)
	for( int J=0; J<3; J++)
	  ofs << strain(e)(I,J) << " "; 
    ofs << std::endl;

    ofs << endl << "POINT_DATA " << _shellNodes.size() << endl
	<< "VECTORS displacements double" << endl;
    /*     << "LOOKUP_TABLE displacements" << endl; */
    //
    // output nodal postions
    pn = _shellNodes.begin();
    for ( ; pn!= _shellNodes.end(); pn ++) {
      Vector3D nodalDisp;
      nodalDisp = (*pn)->point() - (*pn)->position();
      ofs << std::setprecision(16) 
	  << nodalDisp(0)
	  << '\t' <<nodalDisp(1)
	  << '\t' <<nodalDisp(2) << std::endl;
    }
    
    ofs << endl << "VECTORS forces double" << endl;
    /*     << "LOOKUP_TABLE displacements" << endl; */
    //
    // output nodal postions
    pn = _shellNodes.begin();
    for ( ; pn!= _shellNodes.end(); pn ++) {
      const tvmet::Vector<double, 4> nodalForce = (*pn)->force();
      ofs << std::setprecision(16) 
	  << nodalForce(0)
	  << '\t' <<nodalForce(1)
	  << '\t' <<nodalForce(2) << std::endl;
    }
    //
    // output for reaction coordinate
    ofs << "SCALARS    reacCoord    double    1" << std::endl;
    ofs << "LOOKUP_TABLE default" << std::endl;

    pn = _shellNodes.begin();
    for ( ; pn!= _shellNodes.end(); pn ++) {
      double rc;
      rc = (*pn)->getPoint(3); 
      ofs << rc << std::endl;
    }

    ofs.close();

    return;
  }

  //! create Alias-Wavefront .obj file
  template < class Material_t >
  void GLBody< Material_t >::printObj(const std::string name) const
  {
#ifdef WITH_MPI
    if(_processorRank!=0) return;
#endif

    std::string fileName = name + ".obj";
    std::ofstream ofs(fileName.c_str(), std::ios::out);
    if (!ofs) {
      std::cout << "can not open output ("
		<< fileName
		<< ") file." << std::endl;
      exit(0);
    }

    // output nodal postions
    ConstFeNodeIterator pn = _shellNodes.begin();
    for ( ; pn!= _shellNodes.end(); pn ++) {
      const Vector3D & nodalPos =  (*pn)->position();
      ofs << "v " << std::setprecision(16) 
	  << std::setw(24) << nodalPos(0) 
	  << std::setw(24) << nodalPos(1) 
	  << std::setw(24) << nodalPos(2) << std::endl;
    }
    
    // output face connectivities
    for(ConstFeElementIterator pe=_shells.begin(); pe!=_shells.end(); pe++) {
      const typename FeElement_t::NodeContainer & pnc = (*pe)->nodes();
      ofs << "f "
	  << std::setw(10) << pnc[0] -> id()
	  << std::setw(10) << pnc[1] -> id()
	  << std::setw(10) << pnc[2] -> id()
	  << std::endl;
    }
    ofs.close();		

    return;
  }
  
} // namespace voom
