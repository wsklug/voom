// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2009 All Rights Reserved
//
//----------------------------------------------------------------------

/*! 
  \file GenericBody.h

  \brief GenericBody is a generic implemention for the Body interface,
  designed to call the compute functions of a generic list of Elements
  connected to a list of shared nodes.

*/

#if !defined(__GenericBody_h__)
#define __GenericBody_h__

#include<blitz/array.h>
#include<vector>
#include "NodeBase.h"
#include "Element.h"
#include "Body.h"

namespace voom
{

  /*!  GenericBody is a generic implemention for the Body interface,
  designed to call the compute functions of a generic list of Elements
  connected to a list of shared nodes.
  */
  class GenericBody : public Body
  {
    
  public:

    //! Default Constructor
    GenericBody();
    
    //! Element and Node Constructor
    GenericBody(ElementContainer & elems, NodeContainer & nodes);

    //! Default Destructor
    virtual ~GenericBody() {};
    
    //! Do mechanics on GenericBody
    void compute( bool f0, bool f1, bool f2 );
    
    //! Print input file for Paraview
    void printParaview(std::string name) const {};

    
  };

} // namespace voom
#endif // __GenericBody_h__
