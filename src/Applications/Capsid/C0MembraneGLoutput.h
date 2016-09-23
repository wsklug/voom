// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          William S. Klug
//                University of California Los Angeles
//                   (C) 2009 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__C0MembraneGLoutput_h__)
#define __C0MembraneGLoutput_h__

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
#include "C0MembraneGL.h"
#include "voom.h"
#include "Node.h"

#ifdef WITH_MPI
#include <mpi.h>
#endif

namespace voom
{

  /*! 
    C0MembraneGLoutput prints out VTK files for a mesh of C0MembraneGL elements
  */
  class C0MembraneGLoutput
  {
  public:

    typedef std::vector< C0MembraneGL* > 	ElementContainer;
    typedef ElementContainer::iterator 		ElementIterator;
    typedef ElementContainer::const_iterator 	ConstElementIterator;

    typedef std::vector< DeformationNode<3>* > 	DefNodeContainer;
    typedef DefNodeContainer::iterator 		DefNodeIterator;
    typedef DefNodeContainer::const_iterator 	ConstDefNodeIterator;

    typedef std::vector< ScalarFieldNode<3>* > 	GLNodeContainer;
    typedef GLNodeContainer::iterator 		GLNodeIterator;
    typedef GLNodeContainer::const_iterator 	ConstGLNodeIterator;

    //! Construct from stuff
    C0MembraneGLoutput() {}
    
    //! virtual destructor
    virtual ~C0MembraneGLoutput() {;}

    virtual void operator()(ElementContainer & elements, 
			    DefNodeContainer & defNodes,
			    GLNodeContainer & glNodes,
			    std::string filename); 
  };  
} // namespace voom

#endif // __C0MembraneGLoutput_h__
