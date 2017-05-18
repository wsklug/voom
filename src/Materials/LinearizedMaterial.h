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
// Revision 1.3  2005/05/23 17:37:31  klug
// Added cvs logging.
//
//----------------------------------------------------------------------

/*! 
  \file LinearizedMaterial.h

  \brief Templated base class for a small strain material object.
  Specific materials can be derived from this.

*/

#if !defined(__LinearizedMaterial_h__)
#define __LinearizedMaterial_h__


#include "voom.h"
#include <iostream> 

namespace voom
{

//!  Templated base class for 3-D or 2-D small strain material objects.

/*!  Specific small-strain material objects can be derived from this
  class which represents the abstract notion of the point-wise
  response of a material under small-strains.  It is designed to be
  used for materials formulated for <i>variational updates</i>, i.e.,
  for which a scalar function \f$w\f$ exists which serves as a
  potential for the cauchy stress:
  \f[
  \sigma_{ij} = \frac{\partial w}{\partial\epsilon_{ij}} 
  \f]
  where \f$\epsilon\f$ is, in the case of elasticity the linearized
  infinitesimal strain \f$\epsilon_{ij}=\frac{1}{2}(u_{i,j}+u_{j,i})\f$,
  and in the case of inelasticity an <i>incremental</i> strain.

  This template defines the interface for working with such material
  objects, through which the user may 
  <ol>
  <li> specify the strain using the <tt>setStrain()</tt> method, 
  <li> access the strain using the <tt>setStrain()</tt> method, 
  <li> request computational updates of the energy density and its 
  first and second derivatives using the <tt>updateState()</tt> method, 
  and
  <li> access the most recently computed values of the energy density 
  and its first and second derivatives using the <tt>energyDensity()</tt>,
  <tt>stress()</tt>, and <tt>moduli()</tt> methods, respectively.
  </ol>
  Data arrays are stored in Voigt form using Vector and Matrix
  containers from the open-source
  <a href=http://tvmet.sourceforge.net>tvmet</a> library.

  In addition, this class provides implementations of two verification
  tests to be used with concrete derived class objects: 
  <ol>
 
  <li> <tt>checkConsistency()</tt> which computes the derivatives of
  the energy density numerically, and checks that they match the
  directly computed values, and 

  <li> <tt>checkIsotropy()</tt> which verifies that a material is
  isotropic by subjecting it to a random rotation and checking the
  response.
  </ol>

  <b>All</b> materials derived from this template should pass the
  consistency check.  Derived isotropic materials should pass the
  isotropy check as well.
*/

template <unsigned int STRAIN_DIMENSION>
class LinearizedMaterial
{
public:
  //! typename for "vectors" in Voigt notation
  typedef typename
  tvmet::Vector<double,STRAIN_DIMENSION> VoigtVector; 
  //! typename for "matrices" in Voigt notation
  typedef typename 
  tvmet::Matrix<double,STRAIN_DIMENSION,STRAIN_DIMENSION> VoigtMatrix;

  //! Default constructor does nothing but verify that the template parameter is either 3 or 6
  LinearizedMaterial() {
    if( STRAIN_DIMENSION != 3 && STRAIN_DIMENSION != 6 ) {
      std::cout << "LinearizedMaterial: Template parameter STRAIN_DIMENSION must be either 3 (planar problems) or 6 (3-D problems)." << std::endl;
      exit(0);
    }
  }

  virtual ~LinearizedMaterial() {}

  //! Modify the strain.
  void setStrain( const VoigtVector & eps ){ _epsilon = eps; };

  //! Method for updating the response of the material to the strain value.
  /*! 
    If <tt>f0==true</tt> compute the energy density.<br> 
    If <tt>f1==true</tt> compute the first derivative of the energy
    density.<br> 
    If <tt>f2==true</tt> compute the second derivative of the energy
    density.<br>
    
    This method is purely virtual, so it must be implemented by
    derived classes.
  */
  virtual void updateState(bool f0, bool f1, bool f2)=0;

  //! Access the energy density.
  virtual double energyDensity() const { return _w; }
  //! Access the strain.
  virtual const VoigtVector & strain() const { return _epsilon; }
  //! Access the stress (1st derivative of w).
  virtual const VoigtVector & stress() const { return _sigma; }
  //! Access the moduli (2nd derivative of w).
  virtual const VoigtMatrix & moduli() const { return _c; }
  
  //! Check that the 1st and 2nd derivatives of the energy density are consistent
  /*!  This method computes the 1st and 2nd derivatives of the energy
    density twice, the first time directly using the
    <tt>updateState()</tt> method, and the second time by numerical
    differentiation with a 2-point formula.  The results should be
    identical to within a tolerance; if so, return <tt>true</tt>,
    otherwise return <tt>false</tt>.
  */
  bool checkConsistency();

  //! Check that a material is isotropic.
  /*!  Generate a random rotation, and compare the responses of the
    unrotated material state with the rotated state.  If they are
    equivalent, return <tt>true</tt>, otherwise <tt>false</tt>.
  */
  bool checkIsotropy();

protected:
  //! Strain energy density.
  double _w; 		  
  //! Small strain components in Voigt form.
  VoigtVector _epsilon;
  //! Cauchy stress components in Voigt form.
  VoigtVector _sigma;
  //! Matrix of elastic (or effective) moduli in Voigt form.
  VoigtMatrix _c;      
};

} //namespace voom

#include "LinearizedMaterial.icc"


#endif //  !defined(__LinearizedMaterial_h__)
