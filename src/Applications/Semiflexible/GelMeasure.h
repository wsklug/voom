// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          Andrew R. Missel
//                University of California Los Angeles
//                   (C) 2010 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__GelMeasure_h__)
#define __GelMeasure_h__

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
#include "AffinityElement.h"
#include "AffinityMeasure.h"
#include "voom.h"
#include "Node.h"
#include "Lbfgsb.h"
#include "Grid.h"

#ifdef WITH_MPI
#include <mpi.h>
#endif

namespace voom
{

  template<int N> class GelMeasure
  {
  public:
    
    typedef  SemiflexibleGel<N> Gel;

    typedef BrownianNode<N> DefNode;
    typedef typename std::vector< DefNode* > DefNodeContainer;
    typedef typename DefNodeContainer::iterator DefNodeIterator;
    typedef typename DefNodeContainer::const_iterator ConstDefNodeIterator;

    typedef DeformationNode<N> BaseDefNode;

    typedef typename Gel::Filament Filament;
    typedef std::vector< Filament* > FilamentContainer;
    typedef typename FilamentContainer::iterator FilamentIterator;
    typedef typename FilamentContainer::const_iterator ConstFilamentIterator;

    typedef tvmet::Vector<double,N> VectorND;

    struct DataPoint {
      VectorND pos;

      double dat1;
      double dat2;

      const VectorND & pos { return pos; }
      
    };

    typedef Grid<Filament,Filament,N> FilGrid;
    typedef Grid<DefNode,BaseDefNode,N> NodeGrid;
    typedef Grid<AffinityElement,AffinityElement,N> AffElementGrid;
    typedef Grid<DataPoint,DataPoint,N> DataGrid;

    typedef typename std::vector< std::pair<double,double> > CorrelationData;

    //! Construct from stuff
    GelMeasure(Gel * gel) {
      _gel = gel;
    }
   
    ~GelMeasure() {}

    // this function creates the grid and computes quantities //
    void setScale(double scale, double headAffScale = -1.0) {
      if(headAffScale > 0.0) {
	
      }

      else {

      }
    }

    void computeCorrelationFunction(double scale, std::string quant1, std::string quant2, double affScale = -1.0, std::string fileName = "") {
      DataGrid data;

      CorrelationData corrData;

      bool firstQuant = false;

      if(quant1.find("nonaffinity (Head)")!=string::npos) {
	
      }

      else if(quant1.find("nonaffinity (rot)")!string::npos) {

      }
      
      else if(quant1.find("long filament length density")!=string::npos) {
	
	
      }

      else if(quant1.find("length density")!=string::npos) {

      }

      else if(quant1.find("bending energy")!=string::npos) {

      }

      else if(quant1.find("stretching energy")!=string::npos) {

      }

      else if(quant1.find("energy")!=string::npos) {

      }

      else {
	std::cout << "
      }
    }

    void affinityMeasure() {

    }
    
    
    void getCorrelations(CorrelationData & data, DataGrid & data, double max) {
      // given a set of data points, compute the correlation function //
    }
			 



  private:
    
    Gel * _gel;
								

  };  


} // namespace voom

#endif // __GelMeasure_h__
