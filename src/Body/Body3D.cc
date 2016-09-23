// -*- C++ -*-
//----------------------------------------------------------------------
//
//                 Melissa M. Gibbons & Luigi Perotti
//                University of California Los Angeles
//                   (C) 2006 All Rights Reserved
//
//----------------------------------------------------------------------

#include <string>
#include <fstream>
#include <blitz/array-impl.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

#ifdef WITH_MPI
#include <mpe.h>
#endif

namespace voom
{

  Body3D::Body3D(vector<Material *> Mat, 
		 const ConnectivityContainer & Connectivities,
		 const DefNodeContainer & DefNodes, 
		 Quadrature<3> & Quad,
		 Shape<3> & Sh,
		 double k)
  {

#ifdef WITH_MPI
  MPI_Comm_size( MPI_COMM_WORLD, &_nProcessors );
  MPI_Comm_rank( MPI_COMM_WORLD, &_processorRank );
#endif

    // Initialize _dof
    _dof = 0;
    for(ConstDefNodeContainerIterator n = DefNodes.begin(); n != DefNodes.end(); n++)
    {
      _nodes.push_back(*n); 
      _dof += (*n)->dof();
    }
        
    std::cout << "Total dof = "<< _dof  << std::endl;
 
    // Create  and store elements
    _elements.reserve( Connectivities.size() );
    assert(Mat.size() == Connectivities.size());

    int el = 0;
    for(ConstConnectivityIterator c = Connectivities.begin(); c != Connectivities.end(); c++, el++) 
    {
      // Build elemental node container
      DefNodeContainer nds;
      
      for(ElementConnectivity::const_iterator e = c->begin(); e != c->end(); e++) {
	nds.push_back(DefNodes[*e] );
      }
      // cout << DefNodes.size() << endl;
      // Create element
      _elements.push_back(new Element3D(nds, Mat[el], Quad, Sh, k) );      
    }
    std::cout << "Total elements = "<< _elements.size()  << std::endl;

  }
  

  void Body3D::reset()
  {
     for(ConstElementIterator e = _elements.begin(); e != _elements.end(); e++)
     {
       dynamic_cast<Element3D* >(*e)->reset();
     }	
  }


  //! compute
  void Body3D::compute( bool f0, bool f1, bool f2 )
  { 
#ifdef WITH_MPI
    MPE_Decomp1d( _elements.size(), _nProcessors, _processorRank, 
		  &eBegin, &eEnd );
    eBegin--;
    MPE_Decomp1d( _capsidElements.size(), _nProcessors, _processorRank, 
		  &sBegin, &sEnd );
    sBegin--;
#endif

    // Initialize energy and forces
    if(f0) _energy = 0.0;
    if(f1) {
      // bug: should allow model to initialize forces in case multiple
      // bodies share some nodes
      /*for(NodeIterator n=_nodes.begin(); n!=_nodes.end(); n++) 
	for(int i=0; i<(*n)->dof(); i++)
	(*n)->setForce(i,0.0);*/
    }
    // Need to zero out stiffness too!!!!!!!!!!

    // compute energy, forces and stiffness matrix in each element
    // loop through all elements
#ifdef _OPENMP	
#pragma omp parallel for 		\
  schedule(static) default(shared) 
#endif	

    for(int e = 0; e < _elements.size(); e++)
    {
      _elements[e]->compute( f0, f1, f2);      
    }

    if(f0)
    { 
      for(int e = 0; e < _elements.size(); e++)
      {
	_energy += _elements[e]->energy();
      }
    }

    return;
  }

  //! Create input file used by Paraview, a 3D viewer
  void Body3D::printParaviewLinearTet(const string name) const
  {
    // will output data to a .vtk file at every increment of "plate" displacement
    int numEls = _elements.size();

#ifdef WITH_MPI
    blitz::Array<double,1> globalsum(blitz::shape(numEls));
    globalsum = 0.0;
    MPI_Reduce(energy.data(), globalsum.data(), energy.size(), 
	       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    blitz::cycleArrays(energy,globalsum);
#endif

#ifdef WITH_MPI
    if(_processorRank!=0) {
      return;
    }
#endif

    std::string fileName = name + ".vtk";
    std::ofstream ofs(fileName.c_str());
    if (!ofs) {
      std::cout << "can not open output ("
		<< fileName
		<< ") file." << std::endl;
      exit(0);
    }
    /*
    ////////////////////////////////////////////////////////////////////
    //
    //    Node Section
    //
    ofs << "# vtk DataFile Version 2.0\n"
	<< "Test example" << std::endl
	<< "ASCII" << std::endl
	<< "DATASET UNSTRUCTURED_GRID" << std::endl
	<< "POINTS  " << _nodes.size() << "  double" << std::endl;
    
    //
    // output nodal postions
    for (ConstNodeIterator pn = _nodes.begin(); pn != _nodes.end(); pn ++) {
      const Vector3D & nodalPos =  dynamic_cast<DeformationNode<3>* >(*pn)->position();
      ofs << std::setprecision(16) 
	  << nodalPos(0) << "  "
	  << nodalPos(1) << "  "
	  << nodalPos(2) << std::endl;
    }

    /////////////////////////////////////////////////////////////////////
    //
    //    Element Section
    //
    ofs << "CELLS  " << _elements.size() << "  "
	<< 5*_elements.size() << std::endl;
    for (ConstElementIterator pe = _elements.begin(); pe != _elements.end(); pe++)
    {
      const Element3D::NodeContainer & pnc = dynamic_cast<Element3D* >(*pe)->nodes();
      ofs << 4 << "  "
	  << std::setw(10) << pnc[0] -> id()
	  << std::setw(10) << pnc[1] -> id()
	  << std::setw(10) << pnc[2] -> id()
	  << std::setw(10) << pnc[3] -> id()
	  << std::endl;
    }
     
    ofs << "CELL_TYPES " << _elements.size() << std::endl;
    for( int i=0; i<_elements.size(); i++ ) {
      ofs << 10 << std::endl;
    }
	
    ofs << endl << "POINT_DATA " << _nodes.size() << endl
	<< "VECTORS displacements float" << endl;
   
    // output nodal postions
  
    for (ConstNodeIterator pn = _nodes.begin(); pn != _nodes.end(); pn ++) {
      Vector3D nodalPoint =  dynamic_cast<DeformationNode<3>* >(*pn)->point(); 
      Vector3D nodalPosition =  dynamic_cast<DeformationNode<3>* >(*pn)->position(); 
      Vector3D nodalDisp = nodalPoint - nodalPosition;

      ofs << std::setprecision(16) 
	  << nodalDisp(0)
	  << '\t' <<nodalDisp(1)
	  << '\t' <<nodalDisp(2) << std::endl;
    }
  */
    ofs.close();

    return;
  }

 //! create input file of output data used by Paraview, a 3D viewer
  void Body3D::printParaviewPostProcess(const string name) const
  {
    // will output data to a .vtk file at every increment of "plate" displacement
    /*
    int numEls = _capsidElements.size();
    

#ifdef WITH_MPI
    blitz::Array<double,1> globalsum(blitz::shape(numEls));
    globalsum = 0.0;
    MPI_Reduce(energy.data(), globalsum.data(), energy.size(), 
	       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    blitz::cycleArrays(energy,globalsum);
#endif

#ifdef WITH_MPI
    if(_processorRank!=0) {
      return;
    }
#endif

    std::string fileName = name + ".vtk";
    std::ofstream ofs(fileName.c_str());
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
	<< "DATASET UNSTRUCTURED_GRID" << std::endl
	<< "POINTS  " << _capsidNodes.size() << "  double" << std::endl;
    
    //
    // output nodal postions
    ConstCapsidNodeIterator pn = _capsidNodes.begin();
    for ( ; pn!= _capsidNodes.end(); pn ++) {
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
    ofs << "CELLS  " << _capsidElements.size() << "  "
	<< 5*_capsidElements.size() << std::endl;
    for (ConstCapsidElementIterator pe = _capsidElements.begin(); pe != _capsidElements.end(); pe++) {
      const typename CapsidElement_t::NodeContainer & pnc = (*pe)->nodes();
      ofs << 4 << "  "
	  << std::setw(10) << pnc[0] -> id()
	  << std::setw(10) << pnc[1] -> id()
	  << std::setw(10) << pnc[2] -> id()
	  << std::setw(10) << pnc[3] -> id()
	  << std::endl;
    }
     
    ofs << "CELL_TYPES " << _capsidElements.size() << std::endl;
    for( int i=0; i<_capsidElements.size(); i++ ) {
      ofs << 10 << std::endl;
    }
	
    //////////////////////////////////////////////////////////////////////////
    //
    //  output color for each element ( corresponding to energy...)
    //
    ofs << "CELL_DATA    " << _capsidElements.size() << std::endl;
   
    // output for von Mises stress
    ofs << "SCALARS    MisesStress    float    1" << std::endl;
    ofs << "LOOKUP_TABLE default" << std::endl;
    for (ConstCapsidElementIterator pe = _capsidElements.begin(); pe != _capsidElements.end(); pe++ ) 
      ofs << (*pe)->CalcMisesStress() << std::endl;
    ofs << std::endl;

    // output for principal strains in each element
    ofs << "VECTORS PrincipalStrains floats" << endl;
    for (ConstCapsidElementIterator pe = _capsidElements.begin(); pe != _capsidElements.end(); pe++ ) {
      Vector3D prinStrain;
      prinStrain = (*pe)->CalcPrincipalStrains();
      ofs << std::setprecision(16)
	  << prinStrain(0)
	  << '\t' << prinStrain(1)
	  << '\t' << prinStrain(2) << std::endl;
    }
    ofs << std::endl;

    ofs << endl << "POINT_DATA " << _capsidNodes.size() << endl
	<< "VECTORS displacements float" << endl;
   
    // output nodal postions
    pn = _capsidNodes.begin();
    for ( ; pn!= _capsidNodes.end(); pn ++) {
      Vector3D nodalDisp;
      nodalDisp = (*pn)->point() - (*pn)->position();
      ofs << std::setprecision(16) 
	  << nodalDisp(0)
	  << '\t' <<nodalDisp(1)
	  << '\t' <<nodalDisp(2) << std::endl;
    }
    
    ofs << endl << "VECTORS forces float" << endl;
   
    // output nodal forces
    pn = _capsidNodes.begin();
    for ( ; pn!= _capsidNodes.end(); pn ++) {
      const Vector3D & nodalForce = (*pn)->force();
      ofs << std::setprecision(16) 
	  << nodalForce(0)
	  << '\t' <<nodalForce(1)
	  << '\t' <<nodalForce(2) << std::endl;
    }
    ofs.close();
    */

    return;
  }

    //! create input file used by Paraview, a 3D viewer
  void Body3D::printParaviewQuadTet(const string name) const
  {
    /*
    // will output data to a .vtk file at every increment of "plate" displacement

    int numEls = _capsidElements.size();
    
    // blitz::Array<double,1> energy(blitz::shape(numEls));

//     energy = 0.0;

// #ifdef WITH_MPI
//     int eBegin=0, eEnd=0;
//     MPE_Decomp1d( _capsidElements.size(), _nProcessors, _processorRank, 
// 		  &eBegin, &eEnd );
//     eBegin--;
// #endif

// #ifdef WITH_MPI
//     for(int e=eBegin; e<eEnd && e<_capsidElements.size(); e++) {
// 	const CapsidElement_t*const*const pe=&(_capsidElements[e]);
// #else
//     ConstCapsidElementIterator pe = _capsidElements.begin();
//     for ( int e = 0; pe != _capsidElements.end(); pe++, e++) {
// #endif
//       int npts=0;
//       typename CapsidElement_t::ConstQuadPointIterator 
// 	p = (*pe)->quadraturePoints().begin();
//       for( ; p != (*pe)->quadraturePoints().end(); p++){
// 	energy(e) += p->material.energyDensity();
// 	npts++;
//       }
//       energy(e) /= (double)( npts );
//     }

#ifdef WITH_MPI
    blitz::Array<double,1> globalsum(blitz::shape(numEls));
    globalsum = 0.0;
    MPI_Reduce(energy.data(), globalsum.data(), energy.size(), 
	       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    blitz::cycleArrays(energy,globalsum);
#endif

#ifdef WITH_MPI
    if(_processorRank!=0) {
      return;
    }
#endif

    std::string fileName = name + ".vtk";
    std::ofstream ofs(fileName.c_str());
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
	<< "DATASET UNSTRUCTURED_GRID" << std::endl
	<< "POINTS  " << _capsidNodes.size() << "  double" << std::endl;
    
    //
    // output nodal postions
    ConstCapsidNodeIterator pn = _capsidNodes.begin();
    for ( ; pn!= _capsidNodes.end(); pn ++) {
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
    ofs << "CELLS  " << _capsidElements.size() << "  "
	<< 11*_capsidElements.size() << std::endl;
    for (ConstCapsidElementIterator pe = _capsidElements.begin(); pe != _capsidElements.end(); pe++) {
      const typename CapsidElement_t::NodeContainer & pnc = (*pe)->nodes();
      ofs << 10 << "  "
	  << std::setw(10) << pnc[0] -> id()
	  << std::setw(10) << pnc[1] -> id()
	  << std::setw(10) << pnc[2] -> id()
	  << std::setw(10) << pnc[3] -> id()
	  << std::setw(10) << pnc[4] -> id()
	  << std::setw(10) << pnc[5] -> id()
	  << std::setw(10) << pnc[6] -> id()
	  << std::setw(10) << pnc[7] -> id()
	  << std::setw(10) << pnc[8] -> id()
	  << std::setw(10) << pnc[9] -> id()
	  << std::endl;
    }
     
    ofs << "CELL_TYPES " << _capsidElements.size() << std::endl;
    for( int i=0; i<_capsidElements.size(); i++ ) {
      ofs << 24 << std::endl;
    }
	
    //////////////////////////////////////////////////////////////////////////
    //
    //  output color for each element ( corresponding to energy...)
    //
    ofs << endl << "POINT_DATA " << _capsidNodes.size() << endl
	<< "VECTORS displacements float" << endl;
    
    // output nodal postions
    pn = _capsidNodes.begin();
    for ( ; pn!= _capsidNodes.end(); pn ++) {
      Vector3D nodalDisp;
      nodalDisp = (*pn)->point() - (*pn)->position();
      ofs << std::setprecision(16) 
	  << nodalDisp(0)
	  << '\t' <<nodalDisp(1)
	  << '\t' <<nodalDisp(2) << std::endl;
    }
    
   //  ofs << endl << "VECTORS forces float" << endl;
   
//     // output nodal forces
//     pn = _capsidNodes.begin();
//     for ( ; pn!= _capsidNodes.end(); pn ++) {
//       const Vector3D & nodalForce = (*pn)->force();
//       ofs << std::setprecision(16) 
// 	  << nodalForce(0)
// 	  << '\t' <<nodalForce(1)
// 	  << '\t' <<nodalForce(2) << std::endl;
//     }


//     ofs << "CELL_DATA    " << _capsidElements.size() << std::endl;
//     //
//     // output for strain energy
//     //ofs << "SCALARS    strainEnergy    float    1" << std::endl;
//     //ofs << "LOOKUP_TABLE default" << std::endl;
//     //for ( int e = 0; e<numEls; e++) 
//     //  ofs << energy(e) << std::endl;
//     //ofs << std::endl;
    
//     // output for von Mises stress
//     ofs << "SCALARS    MisesStress    float    1" << std::endl;
//     ofs << "LOOKUP_TABLE default" << std::endl;
//     for (ConstCapsidElementIterator pe = _capsidElements.begin(); pe != _capsidElements.end(); pe++ ) 
//       ofs << (*pe)->CalcMisesStress() << std::endl;
//     ofs << std::endl;

    // output for stress tensors
    //ofs << "TENSORS CauchyStress float" << std::endl;
    //for (ConstCapsidElementIterator pe = _capsidElements.begin(); pe != _capsidElements.end(); pe++) {
    // Matrix3D cauchy = (*pe)->cauchyStress();
    //ofs << cauchy(0,0) << " " << cauchy(0,1) << " " << cauchy(0,2) << std::endl
    //  << cauchy(1,0) << " " << cauchy(1,1) << " " << cauchy(1,2) << std::endl
    //  << cauchy(2,0) << " " << cauchy(2,1) << " " << cauchy(2,2) << std::endl
    //  << std::endl;
    //}
    
    ofs.close();
    */
    return;
  }

  //! set I1mostProb (most probable invariants) etc. for the experimental material which helps to find the "true" reference state
  void Body3D::setMPinv(std::vector<double> & I1mp, std::vector<double> & I2mp, std::vector<double> & Jmp )
  {	
    int el = 0;
    for(ElementIterator pe = _elements.begin(); pe != _elements.end(); pe++, el++)
    {
      Element3D::QuadPointContainer quad = dynamic_cast<Element3D *>(*pe)->quadPoints(); 
      // assert(dynamic_cast<Element3D *>(*pe) != NULL);
      for (Element3D::ConstQuadPointIterator pq = quad.begin(); pq != quad.end(); pq++)
      {
	dynamic_cast<HomogMP*>(pq->material)->setMPinv(I1mp[el],I2mp[el],Jmp[el]);
	// assert(dynamic_cast<HomogMP*>(pq->material) != NULL);
      }
    }
     
    return;
  }


} // namespace voom
