#include "../EntropicSpring.h"
#include "../Body/SemiflexibleGel.h"
#include "../PeriodicBox.h"
#include <string>
#include <iostream>
#include <fstream>

using namespace voom;
using std::cout;
using std::setw;
using std::endl;

int main() {
  
//   tvmet::Matrix<double, 2, 2> xn;
//   xn = 
//     0.0, 0.785,
//     0.785, 0.0;

  ranlib::Normal<double> rngang(0,M_PI/2.0);
  rngang.seed((unsigned int)time(0));
  
  double ang = rngang.random();

  ranlib::Normal<double> rnglen(.1,.02);
  rnglen.seed((unsigned int)time(0));

  typedef DeformationNode<2> Node_t;
  typedef std::vector<Node_t*>::iterator NodeIterator;
  std::vector< Node_t * > nodes;
  Node_t::Point Xold = 50.0,50.0;
  for( int a=0; a<10; a++ ) {

    unsigned id = a;
    Node_t::Point X;
    double len = rnglen.random();
    assert(len>0.0);
    X(0) = Xold + cos(ang)*len;
    X(1) = Xold + sin(ang)*len;
    Xold = X;

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

  SemiflexibleGel<2>* gel = new SemiflexibleGel<2>();
  PeriodicBox* box = new PeriodicBox(100.0,100.0);
  gel->setBox(box);
  //cout << "Creating spring." << endl;
  
  double kT = .0041;
  double kC = 5.0*.0041;

  gel->addFilament();

 //  tvmet::Vector<double,2> sep;
//   sep = nodes[1]->position() - nodes[0]->position();
//   double dist = norm2(sep);

//   voom::EntropicSpring<2> spring( nodes[0], nodes[1], kC,  5.0, dist , 1000.0 ); 

//   spring.reportSpringProperties();

//   std::cout << "Please enter the name of a file in which to store the entropic spring test data: ";
//   char fn[256];
//   std::cin >> fn;
  
//   std::ofstream outFile(fn);
//   outFile << "#L\tE\tF\tk" << std::endl;
  
//   double minlen = dist*.97;
//   double maxlen = dist*1.03;
//   double step = dist/1000.0;
//   int nSteps = (int)((maxlen-minlen)/step);
//   for(int ns=0; ns<nSteps; ns++) {
//     double L = minlen + ns*step;
//     std::vector<double> res = spring.check(L);
//     outFile << L << "\t" << res[0] << "\t" << res[1] << "\t" << res[2] << std::endl;
//   }

//   outFile.close();

//   cout << "Checking consistency of order = 3." << endl;
//   spring.checkConsistency();
//   for( NodeIterator n=nodes.begin(); n!=nodes.end(); n++) {
//     for(int i=0; i<(*n)->dof(); i++) (*n)->setForce(i,0.0);
//   }
    
//     spring.compute(false,true,false);
//     cout << setw(6) << "Node"
//         << setw(12) << "x"
//         << setw(12) << "y"
//         << setw(12) << "fx"
//         << setw(12) << "fy" 
//         << setw(12) << "|f|"
//         << setw(12) << "U"
//         << endl;
    
//     for( int a=0; a<nodes.size(); a++) {
//       cout << setw(6) << a 
//           << setw(12) << nodes[a]->getPoint(0)
//           << setw(12) << nodes[a]->getPoint(1)
//           << setw(12) << nodes[a]->getForce(0)
//           << setw(12) << nodes[a]->getForce(1) 
//           << setw(12) << tvmet::norm2( nodes[a]->force() )
//           << setw(12) << spring.energy()<<endl;
//     }

//   fitOrder = 5;
//   voom::EntropicSpring<2> spring5( nodes[0], nodes[1], double kT, double kC, int fitOrder, double Larc, bool coeffsKnown ); 
  cout << "Bye bye." << endl;
  return 0;

}
