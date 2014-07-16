// -*- C++ -*-
//----------------------------------------------------------------------
//
//                    William S. Klug & Feng Feng
//                University of California Los Angeles
//                   (C) 2004 All Rights Reserved
//
//----------------------------------------------------------------------
//
// $Log$
// Revision 1.18  2005/10/22 19:16:46  klug
// Separated Node and NodeBase.
//
// Revision 1.17  2005/08/22 22:22:18  klug
// Assembly shifted from elements to nodes.  Bodies now compute energy.
// Model::setField renamed less ambiguously to putField.
//
// Revision 1.16  2005/05/25 02:14:37  klug
// Removed Constraint.h Potential.h includes.
//
// Revision 1.15  2005/05/23 17:54:15  klug
// New interface with data storage pushed to solver.
//
//----------------------------------------------------------------------

/*! 
  \file Model.h

  \brief Model is a class for a Finite Element model consisting of
  elements, potentials and constraints.

*/

#if !defined(__Model_h__)
#define __Model_h__

#include<blitz/array.h>
#include<vector>
#include "voom.h"
#include "NodeBase.h"
#include "Body.h"
#include "Element.h"
#include "Constraint.h"

#ifdef WITH_MPI
#include <mpi.h>
#endif

namespace voom
{

  /*!  Class for a Finite Element model consisting of elements,
    potentials and constraints.
  */
  class Model
  {
 
  public:

    typedef std::vector< Body* > 			BodyContainer;
    typedef std::vector< Body* >::iterator 		BodyIterator;
    typedef std::vector< Body* >::const_iterator 	ConstBodyIterator;
    
    typedef std::vector< NodeBase* > 			NodeContainer;
    typedef std::vector< NodeBase* >::iterator 		NodeIterator;
    typedef std::vector< NodeBase* >::const_iterator 	ConstNodeIterator;

    typedef std::vector< Constraint* > ConstraintContainer;
    typedef ConstraintContainer::iterator ConstraintIterator;
    typedef ConstraintContainer::const_iterator ConstConstraintIterator;
    
    //! Default Constructor
    Model() {};

    Model( const BodyContainer & bodies, const NodeContainer & nodes );

    Model( const NodeContainer & nodes );

    template<class Solver_t>
    void getField(Solver_t & solver) const;

    template<class Solver_t>
    void putField(const Solver_t & solver);

    template<class Solver_t>
    void addField(const Solver_t & solver);

    //! Do mechanics, send data to solver
    template<class Solver_t>
    void computeAndAssemble( Solver_t & solver, bool f0, bool f1, bool f2 );

    //! Get the number of degrees of freedom in the model
    const int dof() const {return _dof;} 

    const BodyContainer & bodies() const {return _bodies;}

    void pushBackBody( Body * bd ) { _bodies.push_back( bd ); }

    //! check consistency
    bool checkConsistency(bool f1=true, bool f2=true);

    //! check consistency
    bool checkRank(const int rank, bool numerical=false);

    void print(std::string name) {      
//       for(ConstBodyIterator b=_bodies.begin(); b!=_bodies.end(); b++)
// 	(*b)->print(name);
      char str[10];
      if(_bodies.size() == 1) {
	_bodies[0]->print(name);
      } else {
	for(int i=0; i<_bodies.size(); i++) {
	  std::string bdname = name;
	  sprintf(str,"-bd%d",i);
	  bdname += str;
	  _bodies[i]->print(bdname);
	}
      }
    }
    
    // Como?  El printQuad no es bueno.  --WSK
//     void printQuad(std::string name) {
//       for(ConstBodyIterator b=_bodies.begin(); b!=_bodies.end(); b++)
// 	(*b)->printQuad(name);
//     }

    void pushBackConstraint( Constraint * c ) { _constraints.push_back( c ); }

  private:

    //! total degree of freedom in the Model
    int _dof;

    //! container of bodies which contribute to the total model energy
    BodyContainer _bodies;

    //! container of the nodes that represent all of the dof in the model
    NodeContainer _nodes;
    // Note: the model used to access nodes via bodies.  Now the
    // master list of active nodes is stored here.  This should
    // eliminate some redundancy in the bodies and elements, and makes
    // much of the model code cleaner.
    //
    // -- WSK 28 Aug 2007

    ConstraintContainer _constraints;

#ifdef WITH_MPI
    int _nProcessors;
    int _processorRank;
#endif

  };

}; // namespace voom

#include "Model.icc"

#endif // __Model_h__
