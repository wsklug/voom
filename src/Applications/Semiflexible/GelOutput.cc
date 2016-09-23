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
    enFile << shear << "\t" << gel->energy() << "\t" << gel->crosslinkenergy() << "\t" << gel->bendingenergy() << "\t" << gel->stretchingenergy() << std::endl;
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
};
