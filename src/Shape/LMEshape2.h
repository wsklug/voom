#if !defined(__LMEshape2__)
#define __LMEshape2__

#include <vector>
#include "voom.h"
#include<tvmet/Vector.h>
#include<tvmet/Matrix.h>
#include "Node.h"
#include "NodeBase.h"
#include "VoomMath.h"
#include "Shape.h"

using namespace std;

namespace voom
{
  // Class for Linear Maximum Entropy shape functions (LME)
  
  class LMEshape2 : public Shape<2> {
  public:
    
    virtual ~LMEshape2() {}
    
    //Temporary constructor
    LMEshape2(double beta, double searchR, double tol, int nItMax);
   
    //! Compute the shape functions and derivatives
    void compute(const CoordinateArray & s) {;};
    bool compute(const CoordinateArray & s, const std::vector<DeformationNode<2>* >& Nodes, bool FindNeigh);

    double checkPartitionOfUnity();
    PositionContainer nodalCoordinates() {
	return _positions;
    }
   
    vector<double>  LME_r(const CoordinateArray & s);
    double LME_pa(const CoordinateArray & s, int a);
    vector<vector<double > > LME_J(const CoordinateArray & s);
   
    vector<vector<double > > DyadicProduct(const vector<double > & a, const vector<double > & b);
    vector<vector<double > > ScalarProduct(const double a, const vector<vector<double > > A);
    void MatrixSum (vector<vector<double > > & A, const vector<vector<double > > & B);
    void MatrixSub (vector<vector<double > > & A, const vector<vector<double > > & B);

  protected:
    // Data

    // Search radius, parameter beta, minimizer lambda, tolerance and maximum iterations
    double _beta;
    double _searchR;
    double _tol;
    int _nItMax;
    vector<double> _lambda;

  };
}
#endif
