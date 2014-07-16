// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          William S. Klug
//                University of California Los Angeles
//                   (C) 2009 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__GelOutput_h__)
#define __GelOutput_h__

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
#include "voom.h"
#include "Node.h"
#include "Lbfgsb.h"

#ifdef WITH_MPI
#include <mpi.h>
#endif

namespace voom
{

  /*! 
    GelOutput prints out VTK files for a mesh of C0MembraneGL elements
  */
  template<int N>
  class GelOutput
  {
  public:
    
    typedef  SemiflexibleGel<N> Gel;

    typedef BrownianNode<N> DefNode;
    typedef typename std::vector< DefNode* > DefNodeContainer;
    typedef typename DefNodeContainer::iterator DefNodeIterator;
    typedef typename DefNodeContainer::const_iterator ConstDefNodeIterator;

    typedef typename Gel::Filament Filament;
    typedef std::vector< Filament* > FilamentContainer;
    typedef typename FilamentContainer::iterator FilamentIterator;
    typedef typename FilamentContainer::const_iterator ConstFilamentIterator;

    typedef tvmet::Vector<double,N> VectorND;

    //! Construct from stuff
    GelOutput() {}
    
    //! virtual destructor
    virtual ~GelOutput() {;}

    virtual void operator()(Gel * gel, std::string filename);
    
    virtual void printGel(Gel* gel, std::string filename);

    virtual void printParamHeader(std::string enFileName, std::string paramHeader);

    virtual void printFieldLabels(std::string enFileName, std::string paramHeader);

    virtual void printEnergies(Gel * gel, std::string enFileName, double shear);

    virtual void printCrosslinkData(Gel * gel, std::string fileName);
    
    virtual void printFilLengthData(Gel * gel, std::string fileName);

    virtual void printNematicData(Gel * gel, std::string fileName);

    //virtual void printSolverData(Gel * gel, Solver * solver, std::string fileName);

  };  


} // namespace voom

#endif // __GelOutput_h__
