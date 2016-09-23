// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2004 All Rights Reserved
//
//----------------------------------------------------------------------
//
// $Log$
// Revision 1.6  2005/10/21 00:28:12  klug
// Removed some old commented-out stuff.  Made updateState virtual to allow
// derived classes.
//
// Revision 1.5  2005/05/23 17:38:39  klug
// Added cvs logging.
//
//----------------------------------------------------------------------

/*! 
  \file SCElastic.h

  \brief Interface for simple Helfrich spontaneous curvature model for
  an elastic lipid bilayer.

*/

#if !defined(__SCElastic_h__)
#define __SCElastic_h__

#include "ShellMaterial.h"
#include "VoomMath.h"


namespace voom
{

  /*!  Class implementing the simple Helfrich spontaneous curvature
    model for an elastic lipid bilayer.
  */
  
  class SCElastic : public ShellMaterial<void*>
  {
  protected:

    double _kC;
    double _kG;
		
    double _C0;
    double _H; // mean curvature
    double _K; // gaussian curvature

  public:

    SCElastic(const double kC, const double kG, const double C0=0.0)
      : _kC(kC), _kG(kG), _C0(C0) {}

    virtual void updateState(bool f0, bool f1, bool f2 );

    double meanCurvature() const {return _H;}
    double gaussianCurvature() const {return _K;}

    double spontaneousCurvature() const { return _C0; }
    double bendingModulus() const { return _kC; }
    double gaussianModulus() const { return _kG; }

    void setBendingModulus( double kC ) { _kC = kC; };
    void setGaussianModulus( double kG ) { _kG = kG; };

    virtual double bendingEnergy() const {return _W;}
    virtual double stretchingEnergy() const {return 0;}
  };
  
} //namespace voom

#endif //  !defined(__SCElastic_h__)
