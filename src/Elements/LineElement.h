/*!
  \file LineElement.h

  \LineElement is derived from element, designed to calculate the line tesion 

*/
#if !defined(__LineElement_h__)
#define __LineElement_h__

#include "Element.h"
#include "Node.h"
#include "VoomMath.h"


namespace voom
{

  class LineElement : public Element
  {
  public:

    //! virtual destructor
    virtual ~LineElement() {;}

    //! constructor
    LineElement(const double chi,
		DeformationNode<3>* node1,
		DeformationNode<3>* node2)
      {
	_chi = chi;
	_node1 = node1;
	_node2 = node2;
      };

    virtual void compute(bool f0, bool f1, bool f2){
      const DeformationNode<3>::Point & x1 = _node1->point();
      const DeformationNode<3>::Point & x2 = _node2->point(); 

      double d = norm2(x1-x2);

      if (f0) {
	_energy = _chi * d;
      }

      //DeformationNode<3>::Point f;

      if (f1) {
	tvmet::Vector< double, 3 > f;
	f = _chi*(x1-x2)/d;

	for (int i=0; i<3; i++){
	  _node1->addForce(i,  f[i]);
	  _node2->addForce(i, -f[i]);
	}	
      }

    }


    double length() { 
      const DeformationNode<3>::Point & x1 = _node1->point();
      const DeformationNode<3>::Point & x2 = _node2->point(); 

      double d = norm2(x1-x2); 
      return d;
    }
      
    const  DeformationNode<3>* node1() { return _node1; }
    const  DeformationNode<3>* node2() { return _node2; }  

  private:

    double _chi; //lineTension
    DeformationNode<3>* _node1;
    DeformationNode<3>* _node2;

  };//end of class

} // namespace voom


#endif // __LineElement_h__
