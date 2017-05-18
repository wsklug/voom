// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          William S. Klug
//                University of California Los Angeles
//                   (C) 2009 All Rights Reserved
//
//----------------------------------------------------------------------

#include "C0MembraneGLoutput.h"

namespace voom {

  void C0MembraneGLoutput::operator()
    (ElementContainer & elements, 
     DefNodeContainer & defNodes,
     GLNodeContainer & glNodes,
     std::string name)
  {

    
    int numShells = elements.size();
    
    blitz::Array<double,1> energy(blitz::shape(numShells));
    blitz::Array<double,1> stretch(blitz::shape(numShells));
    blitz::Array<Tensor3D,1> strain(numShells);
    blitz::Array<Tensor3D,1> stress(numShells);


    energy = 0.0;
    stretch = 0.0;
    strain = Tensor3D(0.0);
    stress = Tensor3D(0.0);


    int e=0;
    for(ConstElementIterator pe=elements.begin(); pe!=elements.end(); pe++){
      int npts=0;
      C0MembraneGL::ConstQuadPointIterator 
	p = (*pe)->quadPoints().begin();
      for( ; p != (*pe)->quadPoints().end(); p++){
	energy(e) += p->material.energyDensity();
	stretch(e) += p->material.stretchingEnergy();

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
		S(I,J) += Scon(alpha,beta)*base_alpha(I)*base_beta(J);
	      }
	    }
	  }
	}
	strain(e) += E;
	stress(e) += S;

	npts++;
      }
      energy(e) /= (double)( npts );
      stretch(e) /= (double)( npts );
      strain(e) /= (double)( npts );
      stress(e) /= (double)( npts );

      e++;
    }


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
	<< "DATASET UNSTRUCTURED_GRID" << std::endl
	<< "POINTS  " << defNodes.size() << "  double" << std::endl;
    
    //
    // output nodal postions
    ConstDefNodeIterator pn = defNodes.begin();
    for ( ; pn!= defNodes.end(); pn ++) {
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

    int nodes_per_element = elements[0]->defNodes().size();

    ofs << "CELLS  " << elements.size() << "  "
	<< (nodes_per_element+1)*(elements.size()) << std::endl;
    for(int e=0; e<elements.size(); e++) {
      const C0MembraneGL * pe = elements[e];
      const C0MembraneGL::DefNodeContainer & pnc = (pe)->defNodes();
      ofs << nodes_per_element << "  ";
      for(int a=0; a<nodes_per_element; a++) {
	ofs << std::setw(10) << pnc[a] -> id();
      }
      ofs << std::endl;
    }
    
    int vtk_cell_type=0;
    if( nodes_per_element == 3 ) {
      // 3-node triangle
      vtk_cell_type = 5;
    } else if( nodes_per_element == 6 ) {
      // 6-node triangle
      vtk_cell_type = 22;
    } 
    
    ofs << "CELL_TYPES  " << elements.size() << std::endl;
    for(int f=0; f<elements.size(); f++) {
      ofs << vtk_cell_type << std::endl;
    }
 	
    //////////////////////////////////////////////////////////////////////////
    //
    //  output color for each element ( corresponding to mean
    //  curvature, energy...)
    //
    ofs << "CELL_DATA    " << elements.size() << std::endl;
    //
    // output for strain energy
    ofs << "SCALARS    energy    double    1" << std::endl;
    ofs << "LOOKUP_TABLE default" << std::endl;
    for ( int e = 0; e<elements.size(); e++) {
      ofs << energy(e) << std::endl;
    }
    ofs << std::endl;

    //
    // output for strain energy
    ofs << "SCALARS    stretchingEnergy    double    1" << std::endl;
    ofs << "LOOKUP_TABLE default" << std::endl;
    for ( int e = 0; e<elements.size(); e++) {
      ofs << stretch(e) << std::endl;
    }
    ofs << std::endl;

//     //
//     // output for strain 
//     ofs << "TENSORS    strain    double" << std::endl;
//     ofs << "LOOKUP_TABLE default" << std::endl;
//     for ( int e = 0; e<numShells; e++) {
//       for( int I=0; I<3; I++)
// 	for( int J=0; J<3; J++)
// 	  ofs << strain(e)(I,J) << " "; 
//       ofs << std::endl;
//     }
//     ofs << std::endl;

    ofs << std::endl << "POINT_DATA " << defNodes.size() << std::endl
	<< "VECTORS displacements double" << std::endl;
    /*     << "LOOKUP_TABLE displacements" << endl; */
    //
    // output nodal postions
    pn = defNodes.begin();
    for ( ; pn!= defNodes.end(); pn ++) {
      Vector3D nodalDisp;
      nodalDisp = (*pn)->point() - (*pn)->position();
      ofs << std::setprecision(16) 
	  << nodalDisp(0)
	  << '\t' <<nodalDisp(1)
	  << '\t' <<nodalDisp(2) << std::endl;
    }
    
    ofs << std::endl << "VECTORS forces double" << std::endl;
    /*     << "LOOKUP_TABLE displacements" << std::endl; */
    //
    // output nodal postions
    pn = defNodes.begin();
    for ( ; pn!= defNodes.end(); pn ++) {
      const Vector3D & nodalForce = (*pn)->force();
      ofs << std::setprecision(16) 
	  << nodalForce(0)
	  << '\t' <<nodalForce(1)
	  << '\t' <<nodalForce(2) << std::endl;
    }

    ofs << std::endl << "SCALARS field double" 
	<< std::endl << "LOOKUP_TABLE default"
	<< std::endl;
    //
    // output nodal postions
    for( ConstGLNodeIterator sn= glNodes.begin(); sn!=glNodes.end(); sn++) {
      ofs << std::setprecision(16) 
	  << (*sn)->point()
	  << std::endl;
    }

    ofs << std::endl << "SCALARS potential double" 
	<< std::endl << "LOOKUP_TABLE default"
	<< std::endl;
    //
    // output nodal postions
    for( ConstGLNodeIterator sn= glNodes.begin(); sn!=glNodes.end(); sn++) {
      ofs << std::setprecision(16) 
	  << (*sn)->force()
	  << std::endl;
    }

    ofs.close();


    return;
  } 
  

};
