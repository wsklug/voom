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
// Revision 1.5  2005/05/23 17:35:44  klug
// Added cvs logging.
//
//----------------------------------------------------------------------

/*! 
  \file Material.h

  \brief Base Class for a material object.  Typical materials should be
  derived from this.

*/

#if !defined(__Material_h__)
#define __Material_h__

#include "voom.h"
#include "VoomMath.h"

namespace voom
{

/*!  Base Class for a material object.  Typical materials should be
  derived from this.
*/

class Material
{
 public:
  void setDeformationGradient( const Tensor3D & F ){
    if( determinant(F) > 0 ) {
      _F=F;
    } else{
      //std::cout << "The determinant of deformation gradient tensor is not greater than zero!" << " J = " << determinant(F) << std::endl;
    }
  };
  virtual void updateState(bool f0, bool f1, bool f2) = 0;
  virtual double energyDensity() const { return _W; }
  virtual const Tensor3D & deformationGradient() const { return _F; }
  virtual const Tensor3D & piolaStress() const { return _P; }
  virtual const Tensor3D & cauchyStress() { return _cauchy; }
  virtual double vonMisesStress() const { return _vMises; }

  virtual Material * copy() = 0;

 protected:
  double _W; // Energy Density
  Tensor3D _F; // Deformation Gradient
  Tensor3D _P; // 1st Piola-Kirchhoff Stress Tensor
  Tensor3D _cauchy; // Cauchy Stress Tensor
  double _vMises; // von Mises stress
  blitz::Array<double,1> _Fp; // Internal (history) variables

};

} //namespace voom

#endif //  !defined(__Material_h__)
