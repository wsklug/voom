//----------------------------------------------------------------------
//
//                    William S. Klug, Feng Feng
//                University of California Los Angeles
//                   (C) 2004 All Rights Reserved
//
//----------------------------------------------------------------------


/*! 
  \file TwoPhaseBody.cc

*/

#include <string>
#include <fstream>
#include <blitz/array-impl.h>

namespace voom
{
  using std::endl;
  using std::cout;

  template<class Material_t>
  void TwoPhaseBody<Material_t>::initializeBody
  (Material_t material,
   ConnectivityContainer & connectivities,
   const NodeContainer & nodes,
   const unsigned quadOrder,
   const int numberOfBoundaries,
   const double pressure,
   const double tension,
   const double chemicalTension,
   const double penaltyVolume,
   const double penaltyArea,
   const double penaltyAreaOne,
   const double viscosity,
   GlobalConstraint volumeConstraint,
   GlobalConstraint areaConstraint,
   GlobalConstraint areaOneConstraint   ) {

    _penaltyArea = penaltyArea;
    _penaltyVolume = penaltyVolume;
    _penaltyAreaOne = penaltyAreaOne;
    _viscosity = viscosity;
    _areaConstraint = areaConstraint;
    _volumeConstraint = volumeConstraint;
    _areaOneConstraint = areaOneConstraint;
    _volume = 0.0;
    _area = 0.0;
    _areaOne = 0.0;
	  
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

    _tensionNode = new MultiplierNode(_nodes.size(), idx, tension );
    
    _chemicalTensionNode = new MultiplierNode(_nodes.size(), idx, chemicalTension);

    std::cout << "TwoPhaseBody::initializeBody(): dof = "<< _dof 
	      << std::endl;

    std::cout << "beginning to create HDS." << std::endl;
    
    _createHDS(connectivities, numberOfBoundaries);
    //_initNodeNeighbors();
    std::cout << "HDS has been created." << std::endl;
	  
    //
    // create elements
    //	std::cout << "begin to create elements..." << std::endl;
    Face_handle fh = _hds.faces_begin();
    for ( ; fh != _hds.faces_end(); fh ++){
      _createElement(material, fh, quadOrder);
    }
    _elements.reserve(_shells.size());
    for(ConstFeElementIterator e=_shells.begin(); e!=_shells.end(); e++) {
      _elements.push_back(*e);
    }
    //	std::cout << "finished creating elements." << std::endl;

    // compute mechanics
    compute(false, false, false);
    _constraintVolume = _volume;
    std::cout << "Set constraint volume: " << _constraintVolume << std::endl;
    _constraintArea = _area;
    std::cout << "Set constraint area: " << _constraintArea << std::endl;
    _constraintAreaOne = _areaOne;
    std::cout << "Set constraint area of the first component: " << _constraintAreaOne << std::endl;
  }

  //! compute
  template< class Material_t >
  void TwoPhaseBody<Material_t>::compute( bool f0, bool f1, bool f2 )
  {


    //
    // compute element areas and volumes
    if( _volumeConstraint != FeElement_t::none || 
	_areaConstraint != FeElement_t::none ||
	_areaOneConstraint != FeElement_t::none) {
      _volume = 0.0;
      _area = 0.0;
      _areaOne = 0.0;

#ifdef _OPENMP	
#pragma omp parallel for
#endif	      
      for(FeElementIterator e=_shells.begin(); e!=_shells.end(); e++) {
        (*e)->compute( false,false,false );
	_volume += (*e)->volume();     
	_area += (*e)->area();
	_areaOne += (*e)->areaOne();
      }

      //std::cout << "Body volume = " << _volume << std::endl;

      // if total volume is less than 0, the element faces may be
      // incorrectly oriented.
      if (_volume <= 0.0){ 
	std::cout << "The current volume is less than Zero and it = "
		  << _volume
		  << std::endl
		  << "Element faces may be incorrectly oriented or mesh may be inverted." 
		  << std::endl;
	print("NonPositiveVolume");

	//exit(0);
      }

      double pressure = -_penaltyVolume * (_volume/_constraintVolume-1.0)/_constraintVolume;
      _pressureNode->setPoint(pressure);
      
      double tension = _penaltyArea * (_area/_constraintArea-1.0)/_constraintArea;
      _tensionNode->setPoint(tension);

      double chemicalTension = _penaltyAreaOne * (_areaOne/_constraintAreaOne-1.0)/_constraintAreaOne;
      _chemicalTensionNode->setPoint(chemicalTension);
    }

    if(f0) {
      _energy = 0.0;

      if( _volumeConstraint == FeElement_t::penalty ) {
	double dv = _volume/_constraintVolume - 1.0;	
	_energy += 0.5 * _penaltyVolume * dv * dv;
	_tempWork1 = 0.5 * _penaltyVolume * dv * dv;
      }

      if( _areaConstraint == FeElement_t::penalty ) {
	double da = _area/_constraintArea - 1.0;
	_energy += 0.5 * _penaltyArea * da * da;
	_tempWork2 = 0.5 * _penaltyArea * da * da;
      }

      if( _areaOneConstraint == FeElement_t::penalty ) {
	double da1 = _areaOne/_constraintAreaOne - 1.0;
	_energy += 0.5 * _penaltyAreaOne * da1 * da1;
      }
    }

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
#pragma omp parallel for
#endif	
    for(Body::ElementIterator e=_elements.begin(); e!=_elements.end(); e++) {
      (*e)->compute( f0, f1, f2 );
    }
		
    if(f0) { 
#ifdef _OPENMP	
#pragma omp parallel for
#endif	
      for(Body::ElementIterator e=_elements.begin(); e!=_elements.end(); e++) {
	_energy += (*e)->energy();
      }
    }

    // viscous regularization part
    if(f0) {
      _viscousEnergy = 0.0;
      for(FeNodeIterator n=_shellNodes.begin(); n!=_shellNodes.end(); n++) {
	typename FeNode_t::Point dx;
	dx = (*n)->point() - (*n)->position();
	_viscousEnergy += 0.5*_viscosity*tvmet::dot(dx,dx);	
      } 
      _energy += _viscousEnergy;
    }
    if(f1) {
      for(FeNodeIterator n=_shellNodes.begin(); n!=_shellNodes.end(); n++) {
	typename FeNode_t::Point f;
	f = _viscosity * ( (*n)->point() - (*n)->position() );
	(*n)->updateForce( f );
      } 
    }
    return;
  }
	

  //! output HDS info.
  template< class Material_t >
  void TwoPhaseBody< Material_t >::printByHDS()
  {
    assert(_hds.is_valid());
    Halfedge_handle hh;  
    Face_handle fh = _hds.faces_begin();
    std::cout << std::endl << std::endl;
    std::cout << " --------- Halfedge Data Structure information -----------" 
	      << std::endl;
    for ( ; fh != _hds.faces_end(); fh ++) {
      std::cout << _hds.index(fh) << std::endl;
      hh = fh -> halfedge();
      std::cout << std::setw(6) 
		<< hh -> vertex() -> getNodePointer() -> id()
		<< std::setw(6) 
		<< hh -> next() -> vertex() -> getNodePointer() -> id()
		<< std::setw(6) 
		<< hh -> prev() -> vertex() -> getNodePointer() -> id()
		<< std::endl;
      }
  }
	
  //! create input file used by Paraview, a 3D viewer
  template < class Material_t >
  void TwoPhaseBody< Material_t >::printParaview(const std::string name) const
  {
    std::cout << "Constraint Energy = " << constraintEnergy() << endl
	      << "total Strain Energy = " << totalStrainEnergy() << endl;
    //<< "work = "    << work() <<endl;
    std::string fileName = 
      //(std::string)(getenv("HOME")) +
      //"/scratch/VTKTwoPhaseI2/" +
      name + ".vtk";				
    std::cout << name << "  " << fileName << std::endl;
    std::ofstream ofs(fileName.c_str(), std::ios::out);
    if (!ofs) {
      std::cout << "can not open output ("
		<< fileName
		<< ") file." << std::endl;
      exit(0);
    }

    //////////////////////////////////////////////////////////////////////////
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

    //////////////////////////////////////////////////////////////////////////
    //
    //    Normal Section     :::::::: Currently, use nodal internal force
    //
    //////////////////////////////////////////////////////////////////////////
    //
    //    Eelement Section
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
    //  curvature, energy, concentration...)
    //
    ofs << "CELL_DATA    " << _shells.size() << std::endl;
    //
    // output for strain energy
    ofs << "SCALARS    strainEnergy    float    1" << std::endl;
    ofs << "LOOKUP_TABLE default" << std::endl;
    ConstFeElementIterator pe = _shells.begin();
    for ( int i = 0; pe != _shells.end(); pe ++, i++) {
      int npts=0;
      double energy=0.0;
      typename FeElement_t::ConstQuadPointIterator 
	p = (*pe)->quadraturePoints().begin();
      for( ; p != (*pe)->quadraturePoints().end(); p++){
	energy += p->material.energyDensity();
	npts++;
      }
      energy /= (double)( npts );
      ofs << energy << std::endl;
    }
    ofs << std::endl;


    //
    // output for mean curvature
    // since one color is mapping into one element, for
    // several guass points integration, we compute their average value
    ofs << "SCALARS    meanCurvature    float    1" << std::endl;
    ofs << "LOOKUP_TABLE meanCurvature    " << /* _shells.size() <<*/ std::endl;		
    pe = _shells.begin();
    for ( int i = 0; pe != _shells.end(); pe ++, i++) {
      double meanCurvature = 0.0;
      int npts=0;
      typename FeElement_t::ConstQuadPointIterator 
	p = (*pe)->quadraturePoints().begin();
      for( ; p != (*pe)->quadraturePoints().end(); p++){
	meanCurvature += p->material.meanCurvature();
	npts++;
      }
	meanCurvature /= (double)( npts );
	ofs << meanCurvature << std::endl;
    }


    ofs << endl << "POINT_DATA " << _shellNodes.size() << endl
	<< "VECTORS displacements float" << endl;
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

    ofs << endl << "VECTORS forces float" << endl;
    /*     << "LOOKUP_TABLE displacements" << endl; */
    //
    // output nodal postions
    pn = _shellNodes.begin();
    for ( ; pn!= _shellNodes.end(); pn ++) {
      const tvmet::Vector<double, 4>  nodalForce = (*pn)->force();
      ofs << std::setprecision(16) 
	  << nodalForce(0)
	  << '\t' <<nodalForce(1)
	  << '\t' <<nodalForce(2) << std::endl;
    }
    ofs << std::endl;

    //
    // output for concentration
    ofs << "SCALARS    concentration    float    1" << std::endl;
    ofs << "LOOKUP_TABLE default" << std::endl;

    pn = _shellNodes.begin();
    for ( ; pn!= _shellNodes.end(); pn ++) {
      double concentration;
      concentration = (*pn)->getPoint(3); 
      ofs << concentration << std::endl;
    }


    /*    pe = _shells.begin();
    for ( int i = 0; pe != _shells.end(); pe ++, i++) {
      int npts=0;
      double concentration=0.0;
      typename FeElement_t::ConstQuadPointIterator 
	p = (*pe)->quadraturePoints().begin();
      for( ; p != (*pe)->quadraturePoints().end(); p++){
	concentration += p->material.concentration();
	npts++;
      }
      concentration /= (double)( npts );
      ofs << concentration << std::endl;
      }*/
    

    ofs.close();
		
  }

  
  
} // namespace voom
