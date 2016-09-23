#include "../EntropicSpring.h"

using namespace voom;
using std::cout;
using std::setw;
using std::endl;

int main() {
  
  tvmet::Matrix<double, 2, 2> xn;
  xn = 
    0.0, 0.5,
    0.5, 0.0;
  
  typedef DeformationNode<2> Node_t;
  typedef std::vector<Node_t*>::iterator NodeIterator;
  std::vector< Node_t * > nodes;
  for( int a=0; a<2; a++ ) {

    unsigned id = a;
    Node_t::Point X;
    Node_t::Point x;
    X = xn(a,0), xn(a,1);

    NodeBase::DofIndexMap idx(2);
    idx[0] = 2*a; idx[1] = 2*a+1; 
    
    Node_t * nd = new Node_t(id,idx,X,X);
    
    nd -> setId(id);
    
    cout << "Read "<<a<<"th node from input: " << endl
	 << nd->id() << endl
	 << "X = " 
	 << std::setw(16) << nd->position()(0)
	 << std::setw(16) << nd->position()(1) 
	 << endl;
    nodes.push_back( nd );
  }

  cout << "Creating spring." << endl;
  
  double kT = 1.0;
  double kC = 15.0;
  double Larc = (1.0/sqrt(2.0))+5.0e-3;
  bool coeffsKnown = true;

  for(  int fitOrder = 3; fitOrder <=5; fitOrder+=2 ) {
    cout << " fitOrder = " << fitOrder << endl;

    voom::EntropicSpring<2> spring( nodes[0], nodes[1], kT,  kC,  fitOrder,  Larc,  coeffsKnown ); 

//   cout << "Checking consistency of order = 3." << endl;
//   spring.checkConsistency();
    for( NodeIterator n=nodes.begin(); n!=nodes.end(); n++) {
      for(int i=0; i<(*n)->dof(); i++) (*n)->setForce(i,0.0);
    }
    
    spring.compute(true,true,false);
    cout << setw(6) << "Node"
	 << setw(12) << "x"
	 << setw(12) << "y"
	 << setw(12) << "fx"
	 << setw(12) << "fy" 
	 << setw(12) << "|f|"
         << setw(12) << "U"
	 << endl;
    
    for( int a=0; a<nodes.size(); a++) {
      cout << setw(6) << a 
	   << setw(12) << nodes[a]->getPoint(0)
	   << setw(12) << nodes[a]->getPoint(1)
	   << setw(12) << nodes[a]->getForce(0)
	   << setw(12) << nodes[a]->getForce(1) 
	   << setw(12) << tvmet::norm2( nodes[a]->force() )
           << setw(12) << spring.energy()<<endl;
    }

//   fitOrder = 5;
//   voom::EntropicSpring<2> spring5( nodes[0], nodes[1], double kT, double kC, int fitOrder, double Larc, bool coeffsKnown ); 
  }
  cout << "Bye bye." << endl;
  return 0;

}
