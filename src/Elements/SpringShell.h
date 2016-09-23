// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2010 All Rights Reserved
//
//----------------------------------------------------------------------

/*! 
  \file SpringShell.h

  \brief SpringShell is derived from element, implementing an angle
         spring between two triangles for shell bending.

*/

#if !defined(__SpringShell_h__)
#define __SpringShell_h__

#include <vector>
#include <cstdio>
#include <ctime>

#include "Element.h"
#include "Node.h"

namespace voom
{

  //! Element implementing an angle spring between two triangles for shell bending.  
  /*!  Nodal arrangement is defined as follows:

    \verbatim

   2-------1
   | (0) / |
   |   /   |
   | / (1) |
   0-------3               \endverbatim

    Energy of the spring is quadratic in the angle
    \f[
         E = \frac{k}{2}(\theta-\theta_0)^2
    \f]

    where 

    \f[
         \cos\theta = x \equiv n_0\cdot n_1 .
    \f]

    The normal for a triangle connecting nodes \f$(a,b,c)\f$ are
    computed as \f$ \vec{n} = \vec{A}/A, \f$ where 

    \f[
    2\vec{A} = \vec{r}_{ab}\times\vec{r}_{ac} = 
    \vec{r}_{a}\times\vec{r}_{b} +
    \vec{r}_{b}\times\vec{r}_{c} +
    \vec{r}_{c}\times\vec{r}_{a} .
    \f]
    
    and \f$A\equiv|\vec{A}|\f$ is the area of the triangle.

    The derivatives of the normal are calculated as
    \f[
       \frac{\partial\vec{n}}{\partial\vec{r}_I} = 
               P \frac{\partial (2\vec{A})}{\partial\vec{r}_I}
    \f]
    where projection operator \f$P\f$ is 
    \f[
        P = \frac{1}{2A}[ I - \vec{n}\otimes\vec{n} ] , 
    \f]

    with

    \f[
     \frac{\partial (2\vec{A})}{\partial\vec{r}_a} = [(\vec{r}_b-\vec{r}_c)\times]^T \quad
     \frac{\partial (2\vec{A})}{\partial\vec{r}_b} = [(\vec{r}_c-\vec{r}_a)\times]^T \quad
     \frac{\partial (2\vec{A})}{\partial\vec{r}_c} = [(\vec{r}_a-\vec{r}_b)\times]^T 
    \f]

    where \f$[\vec{v}\times]\f$ is the dual of \f$\vec{v}\f$, i.e.,
    the skew symmetric tensor with the property
    \f$[\vec{v}\times]\vec{w} = \vec{v}\times\vec{w}\f$, for any
    vector \f$\vec{w}\f$.

    Accordingly the internal forces are computed as 

    \f[
    \frac{\partial E}{\partial \vec{r}_I} = 
      \frac{dE}{d\theta}\frac{d\theta}{dx}\Big( 
        \vec{n}^{(0)}\cdot\frac{\partial\vec{n}^{(1)}}{\partial\vec{r}_I} +
	\vec{n}^{(1)}\cdot\frac{\partial\vec{n}^{(0)}}{\partial\vec{r}_I} 
      \Big) 
    \f]

    where \f$ \theta=\arccos(x) \f$ and 
    \f$ x =\vec{n}^{(0)}\cdot\vec{n}^{(1)}\f$ give

    \f[
    \frac{d\theta}{dx} = \frac{d}{dx}\arccos(x) = \frac{-1}{\sqrt{1-x^2}} .
    \f]

    Explicitly then, we have the internal forces as

    \f[
    \vec{f}_{0} = \frac{\partial E}{\partial \vec{r}_0} = 
    \frac{dE}{d\theta}\frac{d\theta}{dx}\big[ 
        (\vec{r}_1-\vec{r}_2)\times(P^{(0)}\vec{n}^{(1)}) +
        (\vec{r}_3-\vec{r}_1)\times(P^{(1)}\vec{n}^{(0)}) 
      \big] 
    \f]

    \f[
    \vec{f}_{1} = \frac{\partial E}{\partial \vec{r}_1} = 
    \frac{dE}{d\theta}\frac{d\theta}{dx}\big[ 
        (\vec{r}_2-\vec{r}_0)\times(P^{(0)}\vec{n}^{(1)}) +
        (\vec{r}_0-\vec{r}_3)\times(P^{(1)}\vec{n}^{(0)}) 
      \big] 
    \f]

    \f[
    \vec{f}_{2} = \frac{\partial E}{\partial \vec{r}_2} = 
    \frac{dE}{d\theta}\frac{d\theta}{dx}\big[ 
        (\vec{r}_0-\vec{r}_1)\times(P^{(0)}\vec{n}^{(1)}) 
      \big] 
    \f]

    \f[
    \vec{f}_{3} = \frac{\partial E}{\partial \vec{r}_3} = 
    \frac{dE}{d\theta}\frac{d\theta}{dx}\big[ 
        (\vec{r}_1-\vec{r}_0)\times(P^{(1)}\vec{n}^{(0)}) 
      \big] 
    \f]

  */
  class SpringShell : public Element
  {

  public:

    typedef DeformationNode<3> DefNode; // nickname for mechanics nodes

    typedef std::vector< DefNode* > DefNodeContainer;
    typedef std::vector< DefNode* >::iterator DefNodeIterator;
    typedef std::vector< DefNode* >::const_iterator ConstDefNodeIterator;

    //! virtual destructor
    virtual ~SpringShell() {;}

    //! Constructor
    SpringShell( const DefNodeContainer & defNodes, 
		 double k, double theta0 = 0.0 );

    const DefNodeContainer & defNodes() const {return _defNodes;}

    //! Do mechanics on element; compute energy, forces, and/or stiffness.
    virtual void compute(bool f0, bool f1, bool f2);
    
    double strainEnergy() const { return _strainEnergy; }

    //! Accessor for element work done by pressure
    double work() const { return 0.0; }
		
    //! Set the reference angle
    void setReferenceAngle(double theta0) {_theta0 = theta0;}

    //! Access reference angle
    double referenceAngle() const {return _theta;}

    //
    //   data
    //
  private:

    DefNodeContainer _defNodes;

    //! strain energy
    double _strainEnergy;

    //! angle between face normals
    double _theta;

    //! reference angle between face normals
    double _theta0;

    //! spring constant
    double _k;

  };  // end of class
} // namespace voom

#endif // __SpringShell_h__
