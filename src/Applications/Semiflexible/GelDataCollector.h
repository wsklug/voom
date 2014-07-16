// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          Andrew R. Missel
//                University of California Los Angeles
//                   (C) 2009 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__GelDataCollector_h__)
#define __GelDataCollector_h__

#include <blitz/array.h>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include "SemiflexibleGel.h"
#include "Grid.h"
#include "voom.h"
#include "Node.h"

#ifdef WITH_MPI
#include <mpi.h>
#endif

namespace voom {

  /*! 
    GelDataCollector includes methods to collect and print gel data
  */

  template<int N>
  class GelDataCollector {
  
  public:
    
    typedef SemiflexibleGel<N> Gel;
    typedef Gel::Filament Fil;
    typedef typename tvmet::Vector<double,N> VectorND; 

    //! Default constructor
    GelDataCollector() {
      _gel = 0;
      
    }

    //! Constructor with pointer to gel
    GelDataCollector(Gel* gel) {
      _gel = gel;
    }
    
    //! virtual destructor
    virtual ~GelDataCollector() {;}

    std::vector<std::pair< VectorND, double > > nematicOrderParameterField(double boxSize) {

      // create the grid //
      VectorND gridSize;
      for(int i=0; i<N; i++) gridSize[i] = boxSize;
      Grid<Fil,Fil,2> * FilGrid = new Grid<Fil,Fil,2>();
      FilGrid->setBox(_gel->box());
      FilGrid->setPosFunc(&Fil::point);
      FilGrid->setGridSpace(gridSize);
      FilGrid->setComputeNeighbors(false);
      
      gridSize = FilGrid->gridSpace();

      // add all filaments to the grid //
      int nFils = _gel->filaments.size();
      for(int j=0; j<nFils; j++) {
	FilGrid->addElem(_gel->filament(j));
      }
      
      // compute the nematic order parameter for each square //
      std::vector< std::set<Fil *> > & fils = FilGrid->elemBoxes();
      int nGridSquares = fils.size();
      for(int k=0; k<nGridSquares; k++) {
	double meanang = 0.0;
	int numf = fils[k].size();
	
      }
      

    }
    


    void printEnergies(std::string enFileName, std::string paramHeader) {
      if(_gel == 0) {
	std::cerr << "Error: gel data collector not initialized with gel; please use member function setGel(Gel * gel)." << std::endl;
	exit(1);
      }
      else {
	std::ofstream 
      }
    }

    void printEnergies(Gel* gel, std::string enFileName, std::string paramHeader) {
      setGel(gel);
      printEnergies(enFileName,paramHeader);
      setGel(0);
    }



  private:
	
    void setGel(Gel * gel) {
      _gel = gel;
    }
    
    Gel * _gel;
    
  };  
  
  
} // namespace voom

#endif // __GelDataCollector_h__
