// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                 (C) 2004-2005 All Rights Reserved
//
//----------------------------------------------------------------------
//
// $Log$
// Revision 1.1  2005/05/23 18:40:05  klug
// Initial checkin.
//
//----------------------------------------------------------------------

#if !defined(__Constraint_h__)
#define __Constraint_h__

#include <vector>
#include "Node.h"
#include "Element.h"

namespace voom
{

  template<int dim_n>
  class MultiplierNodeConstraint : public ElementBase 
  {
  public:
    
    typedef tvmet::Vector<double, dim_n> Vector;
    MultiplierNodeConstraint( const DeformationNode<dim_n> * defNode, 
			      const MultiplierNode * multNode, 
			      const Vector & direction,
			      double pointValue ) {
      _defNode = defNode;
      _multNode = multNode;
      _direction = direction;
      _value = pointValue;

      double norm=tvmet::norm2(_direction);
      if( std::abs(norm-1.0) > 1.0e-12 ) _direction /= norm;

      _index.resize(dim_n+1);
      _index[0] = _multNode->index()[0];
      for(int i=0; i<dim_n; i++) 
	_index[i+1] = _defNode->index()[i];

      _energy = 0.0;
      _forces.resize(dim_n+1);
      _forces = 0.0;
      _stiffness.resize(dim_n+1,dim_n+1);
      _stiffness = 0.0;
      for(int i=0; i<dim_n; i++) {
	_stiffness(0,i+1) = _stiffness(i+1,0) = -_direction(i);
      }
    }

    void compute(bool f0, bool f1, bool f2, bool updateDof=true) {
      const tvmet::Vector<double,dim_n> & point = _defNode->point(); 
      double g = tvmet::dot( _direction, point ) - _value;
      double lambda = _multNode->point(); 
      if(f0) {
	_energy = - lambda*g;
      }

      if(f1) {
	_forces(0) = -g;
	for(int i=0; i<dim_n; i++) _forces(i+1) = -lambda*_direction(i);
      }
      
    }

    bool checkConsistency();

  private:
    const DeformationNode<dim_n> * _defNode;
    const MultiplierNode * _multNode; 
    Vector _direction;
    double _value;

  };

  template<int dim_n>
  bool MultiplierNodeConstraint<dim_n>::checkConsistency() {

    ElementVector forces_n(_forces.shape());
    ElementMatrix stiffness_n(_stiffness.shape());
    forces_n = 0.0;
    stiffness_n = 0.0;
    
    MultiplierNode * mn = const_cast<MultiplierNode*>(_multNode);
    DeformationNode<dim_n> * dn = const_cast<DeformationNode<dim_n>*>(_defNode);

    srand(time(0));
    mn->addPoint( 0, 0.05*(double)rand()/RAND_MAX );
    for(int i=0; i<dim_n; i++) 
      dn->addPoint(i, 0.05*(double)rand()/RAND_MAX );
    
    double h=1.0e-6;
    mn->addPoint(0,h);
    compute(true,true,false);
    forces_n(0) = _energy;
    blitz::Range all = blitz::Range::all();
    stiffness_n(0,all) = _forces;

    mn->addPoint(0,-2.0*h);
    compute(true,true,false);
    forces_n(0) -= _energy;
    stiffness_n(0,all) -= _forces;
      
    mn->addPoint(0,h);

    for(int i=0; i<dim_n; i++) {
      dn->addPoint(i,h);
      compute(true,true,false);
      forces_n(i+1) = _energy;
      stiffness_n(i+1,all) = _forces;

      dn->addPoint(i,-2.0*h);
      compute(true,true,false);
      forces_n(i+1) -= _energy;
      stiffness_n(i+1,all) -= _forces;
      
      dn->addPoint(i,h);
    } 
    forces_n /= 2.0*h;
    stiffness_n /= 2.0*h;
   
    compute(false,true,true);

    double Ferror=0.0;
    double Fnorm =0.0;
    double Kerror=0.0;
    double Knorm =0.0;
    double tol=1.0e-10;
    const int nDOF = dim_n+1;
    for(int i=0; i<nDOF; i++){
      const double & f = _forces(i);
      const double & fn = forces_n(i);
      Ferror = std::max(std::abs(f-fn),Ferror);
      Fnorm += (f)*(f);
      for(int j=0; j<nDOF; j++){
	const double & k = _stiffness(i,j);
	const double & kn = stiffness_n(i,j);
	Kerror = std::max(std::abs(k-kn),Kerror);
	Knorm += (k)*(k);
      }
    }
    Fnorm = sqrt(Fnorm);
    Knorm = sqrt(Knorm);
    
    std::cout <<std::setw(18)<<"Ferror ="<<std::setw(12)<<Ferror 
	      <<std::setw(18)<<"tol*Fnorm =" <<std::setw(12)<<tol*Fnorm
	      <<std::endl
	      <<std::setw(18)<<"Kerror ="<<std::setw(12)<<Kerror 
	      <<std::setw(18)<<"tol*Knorm =" <<std::setw(12)<<tol*Knorm
	      <<std::endl;
    if( Ferror < tol*Fnorm ||  Kerror < tol*Knorm) {
      std::cout << "Constraint consistency check PASSED!"<<std::endl;
      return true;
    }
    std::cout << "Constraint consistency check FAILED!"<<std::endl
	      <<"forces = "<<_forces<<std::endl
	      <<"forces_n = "<<forces_n<<std::endl
	      <<"stiffness = "<<_stiffness<<std::endl
	      <<"stiffness_n = "<<stiffness_n<<std::endl;
    return false;
  }

} // end namespace voom
#endif // __Constraint_h__
