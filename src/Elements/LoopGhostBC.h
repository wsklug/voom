// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          William S. Klug
//                University of California Los Angeles
//                    (C) 2008 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__LoopGhostBC_h__)
#define __LoopGhostBC_h__

#include "Node.h"
#include "Constraint.h"

namespace voom
{
  //! Free boundary condition for a ghost boundary node in a LoopShellBody

  /*! Consider a free boundary edge connecting nodes 0 and 2 with
      opposite interior node 1.  Following Cirak, et al. (IJNME,
      2000), we modify the subdivision scheme for the edge by
      introducing a ghost node 3 outside the boundary.  For the edge
      to be traction free, the position of 3 is constrained to be

      \f[
      \bm{x}_3 = \bm{x}_0 + \bm{x}_2 - \bm{x}_1
      \f]

            3
           / \
         /     \
        0-------2
         \     /
           \ /  
            1

   */
  class LoopGhostBC : public Constraint
  {
  public:

    typedef DeformationNode<3> Node_t;
    
    LoopGhostBC(Node_t * N0, Node_t * N1, Node_t * N2, Node_t *N3);

    virtual void predict();
    virtual void correct();

  private:

    Node_t * _N0;
    Node_t * _N1;
    Node_t * _N2;
    Node_t * _N3;
  };

} // namespace voom

#endif // __Constraint_h__
