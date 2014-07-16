// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          William S. Klug
//                University of California Los Angeles
//                   (C) 2009 All Rights Reserved
//
//----------------------------------------------------------------------

#include "GelOutput.h"

namespace voom {

  template<>
  void GelOutput<2>::operator()(Gel * gel, std::string name) {

    //
    // open file
    //
    std::string fileName = name + ".vtp";
    std::ofstream ofs(fileName.c_str());
    if (!ofs) {
      std::cout << "Error: can not open paraview output file "
		<< fileName
		<< std::endl;
      return;
    }
    
    //
    //    Node position data
    //

    ofs << "<?xml version=\"1.0\"?>" << endl;

    ofs << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << endl;

    ofs << "<PolyData>" << endl;

    // count up all filament nodes
    int nNodes = 0;
    for( Gel::ConstFilamentIterator f=gel->filaments().begin(); f!=gel->filaments().end(); f++ ) {
      nNodes += (*f)->nodes.size();
    }
    

//     // nodes on crosslinks
    for( Gel::ConstCrosslinkIterator c=gel->crosslinks().begin(); c!=gel->crosslinks().end(); c++ ) {
       nNodes += 2;
    }

    std::vector< Vector2D > points;
    points.reserve(2*nNodes);

    std::vector< std::vector<int> > lines;
    lines.reserve(2*nNodes);

    // loop through filaments and crosslinks, adding points and line
    // connectivities to lists

    for(Gel::ConstFilamentIterator f=gel->filaments().begin(); f!=gel->filaments().end(); f++) {
      for(int a=0; a < (*f)->bonds.size(); a++) {
	
	const Vector2D & x1 = (*f)->nodes[a  ]->point() ;
	const Vector2D & x2 = (*f)->nodes[a+1]->point();
	Vector2D y1(x1), y2(x2);

	(gel->box())->mapPoint(y1);
	(gel->box())->mapPoint(y2);

	Vector2D u(x2-x1);
	gel->box()->mapDistance(u);

	if( norm2(u) == norm2(y2-y1) ) {  
	  // points are in the same periodic domain, so just plot their
	  // principle images
	  points.push_back( y1 );
	  points.push_back( y2 );
	  std::vector< int > l(2);
	  l[0] = points.size()-2; 
	  l[1] = l[0]+1;
	  lines.push_back( l );
	} else {
	  // if bond is split across box boundary, so print two copies
	  // of it, one on each side of the boundary.

	  // version of the bond starting at y1
	  points.push_back( y1 );
	  y1+=u;
	  points.push_back( y1 );
	  std::vector< int > l(2);
	  l[0] = points.size()-2; 
	  l[1] = l[0]+1;
	  lines.push_back( l );

	  // version of the bond ending at y2
	  points.push_back( y2 );
	  y2-=u;
	  points.push_back( y2 );
	  l[0] = points.size()-2; 
	  l[1] = l[0]+1;
	  lines.push_back( l );
	  
	}
      }
    }
    
    for(Gel::ConstCrosslinkIterator c=gel->crosslinks().begin(); c!=gel->crosslinks().end(); c++) {
      const Vector2D & x1 = (*c)->getNode(0)->point() ;
      const Vector2D & x2 = (*c)->getNode(1)->point();
      Vector2D y1(x1), y2(x2);
      
      gel->box()->mapPoint(y1);
      gel->box()->mapPoint(y2);
      
      Vector2D u(x2-x1);
      gel->box()->mapDistance(u);

      if( norm2(u) == norm2(y2-y1) ) {  
	// points are in the same periodic domain, so just plot their
	// principle images
	points.push_back( y1 );
	points.push_back( y2 );
	
	std::vector< int > l(2);
	l[0] = points.size()-2; 
	l[1] = l[0]+1;
	lines.push_back( l );
	
      } else {
	  // if link is split across box boundary, so print two copies
	  // of it, one on each side of the boundary.

	  // version of the bond starting at y1
	  points.push_back( y1 );
	  y1+=u;
	  points.push_back( y1 );
	  std::vector< int > l(2);
	  l[0] = points.size()-2; 
	  l[1] = l[0]+1;
	  lines.push_back( l );

	  // version of the bond ending at y2
	  points.push_back( y2 );
	  y2-=u;
	  points.push_back( y2 );
	  l[0] = points.size()-2; 
	  l[1] = l[0]+1;
	  lines.push_back( l );

      }
    }
    
    // print stuff
    ofs << "<Piece NumberOfPoints=\"" << points.size() << "\" "
	<< "NumberOfVerts=\"0\" " 
	<< "NumberOfLines=\"" << lines.size() << "\" "
	<< "NumberOfStrips=\"0\" NumberOfPolys=\"0\">"
	<< endl;

    ofs << "<Points>" << endl
	<< "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\" >" << endl;
    for(int a=0; a<points.size(); a++) {
      ofs << std::setprecision(16)
	  << points[a](0) << " 0 " << points[a](1) << endl;
    }
    ofs << "</DataArray>" << endl;

    ofs << "</Points>" << endl;

    ofs << "<Lines>" << endl
	<< "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"1e+299\" RangeMax=\"-1e+299\">"
	<< endl;
    for(int e=0; e<lines.size(); e++) {
      ofs << lines[e][0] << " " << lines[e][1] << endl;
    }
    ofs << "</DataArray>" << endl
	<< "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\" RangeMin=\"1e+299\" RangeMax=\"-1e+299\">"
	<< endl;
    for(int e=0; e<lines.size(); e++) {
      ofs << 2*(e+1) << " " << endl;
    }

    ofs << "</DataArray>" << endl
	<< "</Lines>" << endl
	<< "</Piece>" << endl
	<< "</PolyData>" << endl
	<< "</VTKFile>" << endl;

    ofs.close();

    return;

  }

  template<>
  void GelOutput<2>::printGel(Gel* gel, std::string filename) {
    //     std::cout << "SemiflexibleGel<N>::printParaview()." << std::endl;
    //
    // open file
    //
    std::string fileName = filename + ".vtk";
    std::ofstream ofs(fileName.c_str());
    if (!ofs) {
      std::cout << "Error: can not open paraview output file "
		<< fileName
		<< std::endl;
      return;
    }
    
    //
    //    Node position data
    //
    
    std::map<DefNode*,int> nodeIDMap;
    std::vector<DefNode*> printNodes;
    std::set<DefNode*> imageNodes;
    std::vector< std::pair<DefNode*,DefNode*> > mainPairs;
    std::vector< std::pair<DefNode*,DefNode*> > imgPairs;
    int nodeCount = 0;

    // std::map<DefNode*,VectorND> & finalPosMap = gel->finalPosMap();

    //double shear = _box->shear();
    LeesEdwards* box = (LeesEdwards*)(gel->box());
    double shear = box->shear();
    box->setShear(0.0);
    int nodeID = 0;
    for(ConstFilamentIterator fi=gel->filaments().begin(); fi!=gel->filaments().end(); fi++) {
      Filament* f1 = *fi;
      int nfNodes = f1->nodes.size();
      for(int nn=0; nn<nfNodes-1; nn++) {
	VectorND mapPos1;
	mapPos1 = f1->nodes[nn]->position();
	box->mapPoint(mapPos1);
// 	VectorND def1;
// 	def1 = f1->nodes[nn]->point() - mapPos1;
// 	_box->setShear(shear);
// 	_box->mapDistance(def1);
// 	_box->setShear(0.0);
// 	def1 += mapPos1;
// 	DefNode* newNode1 = new BrownianNode<2>(nodeID,f1->nodes[nn]->index(),mapPoint,def1);
	nodeIDMap.insert(pair<DefNode*,int>(f1->nodes[nn],nodeID));
	nodeID++;
	nodeCount++;
	printNodes.push_back(f1->nodes[nn]);
	VectorND mapPos2;
	mapPos2 = f1->nodes[nn+1]->position();
	box->mapPoint(mapPos2);
	VectorND diff;
	diff = mapPos2 - mapPos1;
	double candist = norm2(diff);
	box->mapDistance(diff);
	double dist = norm2(diff);
	// if bond crosses boundary, make image nodes //
	if(fabs(dist-candist) > 1.0e-6) {
	  VectorND img1pos;
	  img1pos = mapPos2 - diff;
	  VectorND disp;
	  disp = f1->nodes[nn]->point() - img1pos;
	  box->setShear(shear);
	  box->mapDistance(disp);
	  box->setShear(0.0);
	  disp += img1pos;
	  DefNode* imgNode1 = new BrownianNode<2>(nodeID,f1->nodes[nn]->index(),img1pos,disp);
	  nodeIDMap.insert(pair<DefNode*,int>(imgNode1,nodeID));
	  nodeID++;
	  nodeCount++;
	  printNodes.push_back(imgNode1);
	  imageNodes.insert(imgNode1);
	  imgPairs.push_back(std::pair<DefNode*,DefNode*>(imgNode1,f1->nodes[nn+1]));
	  VectorND img2pos;
	  img2pos = mapPos1 + diff;
	  VectorND disp2;
	  disp2 = f1->nodes[nn+1]->point() - img2pos;
	  box->setShear(shear);
	  box->mapDistance(disp2);
	  box->setShear(0.0);
	  disp2 += img2pos;
	  DefNode* imgNode2 = new BrownianNode<2>(f1->nodes[nn+1]->id(),f1->nodes[nn+1]->index(),img2pos,disp2);
	  nodeIDMap.insert(pair<DefNode*,int>(imgNode2,nodeID));
	  nodeID++;
	  nodeCount++;
	  printNodes.push_back(imgNode2);
	  imageNodes.insert(imgNode2);
	  mainPairs.push_back(pair<DefNode*,DefNode*>(f1->nodes[nn],imgNode2));
	}
	// if not, just add node pair //
	else {
	  mainPairs.push_back(pair<DefNode*,DefNode*>(f1->nodes[nn],f1->nodes[nn+1]));
	}
      }
      // now add last node //
      nodeIDMap.insert(pair<DefNode*,int>(f1->nodes[nfNodes-1],nodeID));
      nodeID++;
      nodeCount++;
      printNodes.push_back(f1->nodes[nfNodes-1]);
    }

    // count up all filament nodes
    int nNodes = nodeCount;

    // now do same image thingy with motors/pinches... //
    //nNodes += _motors.size()*2;
    //nNodes += _pinches.size()*2;

    //Mo
//     for( ConstCrosslinkIterator c=gel->crosslinks().begin(); c!=_crosslinks.end(); c++ ) {
//       const DefNodeContainer & crosslinknodes = (*c)-> getNodes();
//       nNodes += crosslinknodes.size();
//     }

    assert(nNodes == printNodes.size());
    
    ofs << "# vtk DataFile Version 2.0\n"
	<< "Test example" << std::endl
	<< "ASCII" << std::endl
	<< "DATASET POLYDATA" << std::endl
	<< "POINTS  " << nNodes << "  double" << std::endl;
    

    // output nodal postions
    for(DefNodeIterator nid=printNodes.begin(); nid!=printNodes.end(); nid++ ) {
      VectorND nodalPos;
      nodalPos = (*nid)->position();
      // if the node is NOT an image node, map its position //
      if(imageNodes.find(*nid) == imageNodes.end()) {
	box->setShear(0.0);
	box->mapPoint(nodalPos);       
      }
      ofs << std::setprecision(16) 
	  << nodalPos(0) << "  "
	  << 0.0 << "  "
	  << nodalPos(1) << std::endl;
    }

//     for (ConstMotorIterator m=_motors.begin(); m!=_motors.end(); m++) {
//       const Vector2D & nodalPos1 = (*m)->getStartPoint();
//       const Vector2D & nodalPos2 = (*m)->getStartPoint();
//       ofs << std::setprecision(16)
// 	  << nodalPos1(0) << "  "
// 	  << 0.0 << "  "
//           << nodalPos1(1) <<std::endl;
//       ofs << std::setprecision(16)
// 	  << nodalPos2(0) << "  "
// 	  << 0.0 << "  "
//           << nodalPos2(1) <<std::endl;
//     }

//     //Mo
//     for( ConstCrosslinkIterator c=_crosslinks.begin(); c!=_crosslinks.end(); c++ ) {
//       const DefNodeContainer & crosslinknodes = (*c)-> getNodes();
//       ConstDefNodeIterator pn = crosslinknodes.begin();
//       for ( ; pn!= crosslinknodes.end(); pn ++) {
//         const Vector2D & nodalPos = (*pn)->position();
//         ofs << std::setprecision(16) 
// 	    << nodalPos(0) << "  "
// 	    << 0.0 << "  "
// 	    << nodalPos(1) << std::endl;
//       }
//     }

    //
    // segment connectivity data
    //

    int nSegments = mainPairs.size() + imgPairs.size();

    int nSegsMeas = 0;
    for( ConstFilamentIterator f=gel->filaments().begin(); f!=gel->filaments().end(); f++ ) {
      nSegsMeas += (*f)->bonds.size();
    }
    assert(nSegsMeas == mainPairs.size());

//     nSegments += _motors.size();

//     //Mo
//     for( ConstCrosslinkIterator c=_crosslinks.begin(); c!=_crosslinks.end(); c++ ) {
//       nSegments++;
//     }

    ofs << "LINES  " << nSegments << "  "
	<< 3*nSegments << std::endl;
    for(std::vector< std::pair<DefNode*,DefNode*> >::iterator np=mainPairs.begin(); np!=mainPairs.end(); np++ ) {
      DefNode* n1 = np->first;
      DefNode* n2 = np->second;
      int indx1 = nodeIDMap[n1];
      int indx2 = nodeIDMap[n2];

      ofs << 2 << "  "
	  << std::setw(10) << indx1
	  << std::setw(10) << indx2
	  << std::endl;
     
    }

    for(std::vector< std::pair<DefNode*,DefNode*> >::iterator np=imgPairs.begin(); np!=imgPairs.end(); np++ ) {
      DefNode* n1 = np->first;
      DefNode* n2 = np->second;
      int indx1 = nodeIDMap[n1];
      int indx2 = nodeIDMap[n2];

      ofs << 2 << "  "
	  << std::setw(10) << indx1
	  << std::setw(10) << indx2
	  << std::endl;
     
    }
//     for( ConstMotorIterator m=_motors.begin(); m!=_motors.end(); m++ ) {
//       ofs << 2 << "  "
// 	  << std::setw(10) << nodeID
// 	  << std::setw(10) << nodeID+1
// 	  << std::endl;
//       nodeID += 2;
//     }
     	
//     //Mo
//     for( ConstCrosslinkIterator c=_crosslinks.begin(); c!=_crosslinks.end(); c++ ) {
//       ofs << 2 << "  "
// 	  << std::setw(10) << nodeID
// 	  << std::setw(10) << nodeID+1
// 	  << std::endl;
//       nodeID+=2;
//     }

    //
    //  output energy for each segment
    //
    ofs << std::endl;
    ofs << "CELL_DATA    " << nSegments << std::endl;

//     ofs << "SCALARS    Energy    double    1" << std::endl;
//     ofs << "LOOKUP_TABLE default" << std::endl;
//     for( ConstFilamentIterator f=_filaments.begin(); f!=_filaments.end(); f++ ) {
//       int firstAngle=0;
//       int lastAngle=(*f)->angles.size()-1;   //What happened to those two variables???  --Mo
//       for(int a=0; a < (*f)->bonds.size(); a++) {
// 	double energy = (*f)->bonds[a]->energy();
// 	// 	energy += 0.5*(*f)->bonds[std::max(a-1,firstAngle)]->energy();
// 	// 	energy += 0.5*(*f)->bonds[std::min(a,lastAngle)]->energy();
// 	ofs << energy << std::endl;
//       }
//     }
//     for( ConstMotorIterator m=_motors.begin(); m!=_motors.end(); m++ ) {
//       ofs << (*m)->energy() << std::endl;
//     }
//     ofs << std::endl;

//     //Mo
//     for( ConstCrosslinkIterator c=_crosslinks.begin(); c!=_crosslinks.end(); c++ ) {
//       double energy = (*c)->energy();
//       ofs << energy << std::endl;
//     }
      
//     ofs << std::endl;

    ofs << "SCALARS    ImageorReal    int    1" << std::endl;
    ofs << "LOOKUP_TABLE default" << std::endl;

    for(int imp=0; imp<mainPairs.size(); imp++) ofs << 1 << std::endl;

    for(int iip=0; iip<imgPairs.size(); iip++) ofs << 0 << std::endl;

    ofs << std::endl;

    ofs << "SCALARS    Energy    double    1" << std::endl;
    ofs << "LOOKUP_TABLE default" << std::endl;
    for( ConstFilamentIterator f=gel->filaments().begin(); f!=gel->filaments().end(); f++ ) {
      int nAngles = (*f)->angles.size();   //What happened to those two variables???  --Mo
      for(int a=0; a < (*f)->bonds.size(); a++) {
	double energy = (*f)->bonds[a]->energy();
	if(a==0 && nAngles>0) energy += .5*(*f)->angles[0]->energy();
	else if(a==(*f)->bonds.size()-1 && nAngles>0) energy += .5*(*f)->angles[nAngles-1]->energy();
	else if(nAngles>0) {
	  energy += 0.5*(*f)->angles[a]->energy();
	  energy += 0.5*(*f)->angles[a-1]->energy();
	}
	ofs << energy << std::endl;
      }
    }

    for(int iip=0; iip<imgPairs.size(); iip++) ofs << -1.0 << std::endl;
 
//     for( ConstMotorIterator m=_motors.begin(); m!=_motors.end(); m++ ) {
//       ofs << (*m)->energy() << std::endl;
//     }
//     ofs << std::endl;

//     //Mo
//     for( ConstCrosslinkIterator c=_crosslinks.begin(); c!=_crosslinks.end(); c++ ) {
//       double energy = (*c)->energy();
//       ofs << energy << std::endl;
//     }
      
    ofs << std::endl;

    ofs << "SCALARS    BendEnergy    double    1" << std::endl;
    ofs << "LOOKUP_TABLE default" << std::endl;
    for( ConstFilamentIterator f=gel->filaments().begin(); f!=gel->filaments().end(); f++ ) {
      int nAngles = (*f)->angles.size();   //What happened to those two variables???  --Mo
      for(int a=0; a < (*f)->bonds.size(); a++) {
	double energy = 0.0;
	if(a==0 && nAngles>0) energy += .5*((*f)->angles[0]->energy());
	else if(a==(*f)->bonds.size()-1 && nAngles>0) energy += 0.5*((*f)->angles[nAngles-1]->energy());
	else if(nAngles>0) {
	  energy += 0.5*((*f)->angles[a]->energy());
	  energy += 0.5*((*f)->angles[a-1]->energy());
	}
	ofs << energy << std::endl;
      }
    }

    for(int iip=0; iip<imgPairs.size(); iip++) ofs << -1.0 << std::endl;

//     for( ConstMotorIterator m=_motors.begin(); m!=_motors.end(); m++ ) {
//       ofs << 0.0 << std::endl;
//     }
//    ofs << std::endl;

    //Mo
//     for( ConstCrosslinkIterator c=_crosslinks.begin(); c!=_crosslinks.end(); c++ ) {
//       double energy = (*c)->energy();
//       ofs << 0.0 << std::endl;
//     }
      
    ofs << std::endl;

    ofs << "SCALARS    StretchEnergy    double    1" << std::endl;
    ofs << "LOOKUP_TABLE default" << std::endl;
    for( ConstFilamentIterator f=gel->filaments().begin(); f!=gel->filaments().end(); f++ ) {
      for(int a=0; a < (*f)->bonds.size(); a++) {
	double energy = (*f)->bonds[a]->energy();
	ofs << energy << std::endl;
      }
    }
    for(int iip=0; iip<imgPairs.size(); iip++) ofs << -1.0 << std::endl;

//     for( ConstMotorIterator m=_motors.begin(); m!=_motors.end(); m++ ) {
//       ofs << (*m)->energy() << std::endl;
//     }
//    ofs << std::endl;

    //Mo
//     for( ConstCrosslinkIterator c=_crosslinks.begin(); c!=_crosslinks.end(); c++ ) {
//       double energy = (*c)->energy();
//       ofs << energy << std::endl;
//     }

    ofs << std::endl;

    ofs << "SCALARS    Strain    double    1" << std::endl;
    ofs << "LOOKUP_TABLE default" << std::endl;
    for( ConstFilamentIterator f=gel->filaments().begin(); f!=gel->filaments().end(); f++ ) {
      for(int a=0; a < (*f)->bonds.size(); a++) {
	double lenchange = (*f)->bonds[a]->strain();
	ofs << lenchange << std::endl;
      }
    }
    for(int iip=0; iip<imgPairs.size(); iip++) ofs << -1.0 << std::endl;

    ofs << std::endl;

    ofs << "SCALARS    Stiffness    double    1" << std::endl;
    ofs << "LOOKUP_TABLE default" << std::endl;
    for( ConstFilamentIterator f=gel->filaments().begin(); f!=gel->filaments().end(); f++ ) {
      for(int a=0; a < (*f)->bonds.size(); a++) {
	double stiff = (*f)->bonds[a]->stiffness();
	ofs << stiff << std::endl;
      }
    }
    for(int iip=0; iip<imgPairs.size(); iip++) ofs << -1.0 << std::endl;

    ofs << std::endl;

    ofs << "SCALARS    StiffnessChange    double    1" << std::endl;
    ofs << "LOOKUP_TABLE default" << std::endl;
    for( ConstFilamentIterator f=gel->filaments().begin(); f!=gel->filaments().end(); f++ ) {
      for(int a=0; a < (*f)->bonds.size(); a++) {
	double stiffChange = (*f)->bonds[a]->stiffnessChange();
	ofs << stiffChange << std::endl;
      }
    }
    for(int iip=0; iip<imgPairs.size(); iip++) ofs << -1.0 << std::endl;

    ofs << std::endl;

    ofs << "SCALARS    Angle    double    1" << std::endl;
    ofs << "LOOKUP_TABLE default" << std::endl;
    for( ConstFilamentIterator f=gel->filaments().begin(); f!=gel->filaments().end(); f++ ) {
      for(int a=0; a < (*f)->bonds.size(); a++) {
	VectorND sep;
	sep = (*f)->nodes[a+1]->position() - (*f)->nodes[a]->position();
	double ang = atan2(sep[1],sep[0]);
	if(ang < -M_PI/2.0) ang += M_PI;
	else if(ang >= M_PI/2.0) ang -= M_PI;
	ofs << ang << std::endl;
      }
    }

    for(int iip=0; iip<imgPairs.size(); iip++) ofs << -1.0 << std::endl;

//     for( ConstMotorIterator m=_motors.begin(); m!=_motors.end(); m++ ) {
//       ofs << 0.0 << std::endl;
//     }
//    ofs << std::endl;

    //Mo
//     for( ConstCrosslinkIterator c=_crosslinks.begin(); c!=_crosslinks.end(); c++ ) {
//       double energy = (*c)->energy();
//       ofs << 0.0 << std::endl;
//     }
      
    ofs << std::endl;

    ofs << "SCALARS    FilLen    double    1" << std::endl;
    ofs << "LOOKUP_TABLE default" << std::endl;
    for( ConstFilamentIterator f=gel->filaments().begin(); f!=gel->filaments().end(); f++ ) {
      VectorND e2e;
      int nNodes = (*f)->nodes.size();
      e2e = (*f)->nodes[nNodes-1]->position() - (*f)->nodes[0]->position();
      double fl = norm2(e2e);
      for(int a=0; a < (*f)->bonds.size(); a++) {
	ofs << fl << std::endl;
      }
    }

    for(int iip=0; iip<imgPairs.size(); iip++) ofs << -1.0 << std::endl;

//     for( ConstMotorIterator m=_motors.begin(); m!=_motors.end(); m++ ) {
//       ofs << 0.0 << std::endl;
//     }
//    ofs << std::endl;

    //Mo
//     for( ConstCrosslinkIterator c=_crosslinks.begin(); c!=_crosslinks.end(); c++ ) {
//       double energy = (*c)->energy();
//       ofs << 0.0 << std::endl;
//     }
      
//    ofs << std::endl;

    ofs << std::endl << "POINT_DATA " << nNodes << std::endl
	<< "VECTORS displacements double" << std::endl;

    // output nodal displacements
    for(DefNodeIterator nid=printNodes.begin(); nid!=printNodes.end(); nid++ ) {
      VectorND nodalPos = (*nid)->position();
      VectorND nodalDisp;
      // if the node is NOT an image node, map its point //
      if(imageNodes.find(*nid) == imageNodes.end()) {
	box->setShear(0.0);
	box->mapPoint(nodalPos);
	nodalDisp = (*nid)->point() - nodalPos;
	box->setShear(shear);
	box->mapDistance(nodalDisp);
	box->setShear(0.0);
      }
      else nodalDisp = (*nid)->point() - (*nid)->position();
      
      ofs << std::setprecision(16) 
	  << nodalDisp(0) << '\t' 
	  << 0.0 << '\t' 
	  << nodalDisp(1)  << std::endl;
      
    }

//     for (ConstMotorIterator m=_motors.begin(); m!=_motors.end(); m++) {
//       Vector2D nodalPos1; 
//       nodalPos1 = (*m)->getEndPoint1()-(*m)->getStartPoint();
//       Vector2D nodalPos2;
//       nodalPos2 = (*m)->getEndPoint2()-(*m)->getStartPoint();
//       ofs << std::setprecision(16)
// 	  << nodalPos1(0) << "  "
// 	  << 0.0 << "  "
//           << nodalPos1(1) <<std::endl;
//       ofs << std::setprecision(16)
// 	  << nodalPos2(0) << "  "
// 	  << 0.0 << "  "
//           << nodalPos2(1) <<std::endl;
//     }
    
    //Mo
//     for( ConstCrosslinkIterator c=_crosslinks.begin(); c!=_crosslinks.end(); c++ ) {
//       const DefNodeContainer & crosslinknodes = (*c)-> getNodes();
//       ConstDefNodeIterator pn = crosslinknodes.begin();
//       for ( ; pn!= crosslinknodes.end(); pn ++) {
//         Vector2D nodalDisp;
//         nodalDisp = (*pn)->point() - (*pn)->position();
//         ofs << std::setprecision(16) 
// 	    << nodalDisp(0) << '\t'
// 	    << 0.0 << '\t'
// 	    << nodalDisp(1) << std::endl;
// 	//cout << " semiflexible " << (*pn) <<' ' << (*pn)->point() <<endl;
//       }
      
//     }


//     ofs << std::endl << "VECTORS finaldisplacements double" << std::endl;
//     // output final nodal displacements
//     for(DefNodeIterator nid=printNodes.begin(); nid!=printNodes.end(); nid++ ) {
//       VectorND nodalPos = (*nid)->position();
//       VectorND nodalDisp;
//       // if the node is NOT an image node, map its point //
//       if(imageNodes.find(*nid) == imageNodes.end()) {
// 	box->setShear(0.0);
// 	box->mapPoint(nodalPos);
// 	nodalDisp = finalPosMap[*nid] - nodalPos;
// 	box->setShear(shear);
// 	box->mapDistance(nodalDisp);
// 	box->setShear(0.0);
//       }
//       else nodalDisp = 0.0,0.0;
      
//       ofs << std::setprecision(16) 
// 	  << nodalDisp(0) << '\t' 
// 	  << 0.0 << '\t' 
// 	  << nodalDisp(1)  << std::endl;
      
//     }

    ofs << std::endl << "VECTORS forces double" << std::endl;
    // output nodal forces
    for(DefNodeIterator nid=printNodes.begin(); nid!=printNodes.end(); nid++ ) {
//       Vector2D & nodalPos = (*nid)->position();
//       Vector2D nodalDisp;
//       // if the node is NOT an image node, map its point //
//       if(imageNodes.find(*nid) == imageNodes.end()) {
// 	_box->setShear(0.0);
// 	_box->mapPoint(nodalPos);
// 	nodalDisp = (*nid)->point() - nodalPos;
// 	_box->setShear(shear);
// 	_box->mapDistance(nodalDisp);
// 	_box->setShear(0.0);
//       }
//       else nodalDisp = (*nid)->point() - (*nid)->position();
      
      Vector2D nodalForce = (*nid)->force();
      if(imageNodes.find(*nid) != imageNodes.end()) {
	nodalForce = 0.0;
      }
      ofs << std::setprecision(16) 
	  << nodalForce(0) << '\t' 
	  << 0.0 << '\t' 
	  << nodalForce(1)  << std::endl;
      
    }

  //   for (ConstMotorIterator m=_motors.begin(); m!=_motors.end(); m++) {
//       ofs << std::setprecision(16)
// 	  << 0.0 << '\t'
// 	  << 0.0 << '\t'
// 	  << 0.0 << std::endl;
//       ofs << std::setprecision(16)
// 	  << 0.0 << '\t'
// 	  << 0.0 << '\t'
// 	  << 0.0 << std::endl;
//     }


//     //Mo
//     for( ConstCrosslinkIterator c=_crosslinks.begin(); c!=_crosslinks.end(); c++ ) {
//       const DefNodeContainer & crosslinknodes = (*c)-> getNodes();
//       ConstDefNodeIterator pn = crosslinknodes.begin();
//       for ( ; pn!= crosslinknodes.end(); pn ++) {
//         const Vector2D & nodalForce = (*pn)->force();
//         ofs << std::setprecision(16) 
// 	    << nodalForce(0) << '\t'
// 	    << 0.0 << '\t'
// 	    << nodalForce(1) << std::endl;
//       }
//     }

    ofs << std::endl;

    ofs << "SCALARS    ImageorReal    int    1" << std::endl;
    ofs << "LOOKUP_TABLE default" << std::endl;

    for(DefNodeIterator dni=printNodes.begin(); dni!=printNodes.end(); dni++) {
      if(imageNodes.find(*dni) == imageNodes.end()) ofs << 1 << std::endl;
      else ofs << 0 << std::endl;
    }

    ofs << std::endl;

    ofs << "SCALARS    Pinch    int    1" << std::endl;
    ofs << "LOOKUP_TABLE default" << std::endl;

    for(DefNodeIterator dni=printNodes.begin(); dni!=printNodes.end(); dni++) {
      if(gel->pinchNodes().find(*dni) != gel->pinchNodes().end()) ofs << 1 << std::endl;
      else ofs << 0 << std::endl;
    }

    ofs << std::endl;
    
    ofs.close();

    printNodes.clear();
    imageNodes.clear();
    mainPairs.clear();
    imgPairs.clear();
    nodeIDMap.clear();

    return;
}

  template<>
  void GelOutput<2>::printParamHeader(std::string enFileName, std::string paramHeader) {
    std::ofstream enFile(enFileName.c_str(),ios::trunc);
    enFile << "#" << paramHeader << std::endl;
    enFile.close();
  }
  
  template<>
  void GelOutput<2>::printFieldLabels(std::string enFileName, std::string fieldLabels) {
    std::ofstream enFile(enFileName.c_str(),ios::app);
    enFile << "#" << fieldLabels << std::endl;
    enFile.close();
  }

  template<>
  void GelOutput<2>::printEnergies(Gel * gel, std::string enFileName, double shear) {
    std::ofstream enFile(enFileName.c_str(),ios::app);
    double vRegEn;
    if(gel->viscReg() == 0) vRegEn = 0.0;
    else vRegEn = gel->viscReg()->energy();
    enFile << shear << "\t" << gel->energy() << "\t" << gel->filenergy() << "\t" << gel->bendingenergy() << "\t" << gel->stretchingenergy()<< "\t" << gel->crosslinkenergy() << "\t" << gel->motorenergy() << "\t" << gel->pinchenergy() << "\t" << vRegEn << std::endl;
    enFile.close();
  }
  
  template<>
  void GelOutput<2>::printCrosslinkData(Gel * gel, std::string fileName) {
    std::ofstream clFile(fileName.c_str());
    clFile << "#Mean crosslink separation = " << gel->getMeanCLsep() << std::endl << "#Dist\tFreq" << std::endl;
    SemiflexibleGel<2>::CrosslinkDistFreq & cldfreqs = gel->getCrossDistro();
    SemiflexibleGel<2>::ClDFIter cldi;
    for(cldi=cldfreqs.begin(); cldi!=cldfreqs.end(); cldi++) {
      clFile << cldi->first << "\t" << cldi->second << std::endl;
    }
    clFile.close();
  }

  template<>
  void GelOutput<2>::printFilLengthData(Gel * gel, std::string fileName) {
    std::ofstream filFile(fileName.c_str());
    filFile << "#Mean filament length = " << gel->getMeanFilLen() << std::endl << "#L\tFreq" << std::endl;
    for(map< double, int >::iterator mapIt=gel->getLengthDistro().begin(); mapIt!=gel->getLengthDistro().end(); mapIt++) {
      filFile << mapIt->first << "\t" << mapIt->second << std::endl;
    }
    filFile.close();
  }
 
  template<>
  void GelOutput<2>::printNematicData(Gel *gel, std::string fileName) {
    std::ofstream nemFile(fileName.c_str());
    nemFile << "#Nematic order parameter = " << gel->getNematicOP() << std::endl << "#theta\tFreq" << std::endl;
    for(map<double,int>::iterator mapNem=gel->getNematicDistro().begin(); mapNem!=gel->getNematicDistro().end(); mapNem++) {
      nemFile << mapNem->first << "\t" << mapNem->second << std::endl;
    }
    nemFile.close();
  }

//   template<>
//   void GelOutput<2>::printSolverData(Solver * solver, std::string fileName) {
//     std::ofstream solFile(fileName.c_str());
//     // get maximum force //
//     double maxForce = 0.0;
//     for(int i=0; i<solver->size(); i++) {
//       double tmpForce = abs(solver->gradient(i));
//       if(tmpForce > maxForce) maxForce = tmpForce;
//     }
//     solFile << solver->function() << "\t" << maxForce << "\t" << solver->iterationNo() << std::endl;
//     solFile.close();
//   }
};
