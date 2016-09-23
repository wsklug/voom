#if !defined(__LMEshape__)
#define __LMEshape__

#include <vector>
#include "voom.h"
#include<tvmet/Vector.h>
#include<tvmet/Matrix.h>
#include "Node.h"
#include "NodeBase.h"
#include "VoomMath.h"

using namespace std;

// External LAPACK routine to invert matrices 
extern "C" void dpotrf_(char *uplo, int *n, double *a, int *lda, int *info);
extern "C" void dpotrs_(char *uplo, int *n, int *nrhs, double *a, int *lda,
                        double *b, int *ldb, int *info);
extern "C" void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv,
                        double *b, int *ldb, int *info);

namespace voom
{
  // Class for Linear Maximum Entropy shape functions (LME)
  
  class LMEshape {
  public:
    typedef std::vector<double> FunctionContainer;
    typedef std::vector<int> NodeNContainer;
    typedef tvmet::Vector<double,3> CoordinateArray;
    
    virtual ~LMEshape() {}
    
    //Temporary constructor
    LMEshape(double beta, double searchR, double tol, int nItMax);
    
    //! Return shape functions and their derivatives
    const FunctionContainer & functions() const {return _functions;}
    
    const FunctionContainer & xderivative() const {return _xderivative;}
    
    const FunctionContainer & yderivative() const {return _yderivative;}
    
    const FunctionContainer & zderivative() const {return _zderivative;}
    
    const NodeNContainer & neighb() const {return _nodes;}
    NodeNContainer & neighb() {return _nodes;}


    //! Compute the shape functions and derivatives
    bool compute(const CoordinateArray & s, const std::vector<DeformationNode<3>* >& nodes, bool FindNeigh, bool f1);

    double checkPartitionOfUnity();

    vector<double>  LME_r(const CoordinateArray & s, const std::vector<DeformationNode<3>* >& nodes);
    double LME_pa(const CoordinateArray & s, const std::vector<DeformationNode<3>* >& nodes, int);
    vector<vector<double > > LME_J(const CoordinateArray & s, const std::vector<DeformationNode<3>* >& nodes);
   

    vector<vector<double > > DyadicProduct(const vector<double > & a, const vector<double > & b);
    vector<vector<double > > ScalarProduct(const double a, const vector<vector<double > > A);
    void MatrixSum (vector<vector<double > > & A, const vector<vector<double > > & B);
    void MatrixSub (vector<vector<double > > & A, const vector<vector<double > > & B);



  protected:
    // Data
    
    //! Shape functions and their derivatives
    FunctionContainer _functions, _xderivative, _yderivative, _zderivative;
    
    //! Node number corresponding to the shape functions calculated
    NodeNContainer _nodes;

    // Search radius, parameter beta, minimizer lambda
    double _beta;
    double _searchR;
    double _tol;
    int _nItMax;
    vector<double> _lambda;

  };
}
#endif
