#include "LMEshape2.h"
#include <math.h>
namespace voom {

  LMEshape2::LMEshape2(double beta, double searchR, double tol, int nItMax): 
                     _beta(beta), _searchR(searchR), _tol(tol), _nItMax(nItMax)
  {
    for (int i = 0; i < 3; i++) {
      _lambda.push_back(0.0);
    }
  };

  bool LMEshape2::compute(const CoordinateArray & s, const std::vector<DeformationNode<2>* >& Nodes, bool FindNeigh)
  {
    // Input data are: the point at which to evaluate the shape function, the nodes in the domain.
    bool Success = true;
    int nNodes = Nodes.size(), i = 0;
    
    if (FindNeigh) 
    {
      // Find Neighborhood
      for(i = 0; i < nNodes; i++) 
	{
	  CoordinateArray x(Nodes[i]->getPosition(0), Nodes[i]->getPosition(1));
	  CoordinateArray temp(x-s);
	  if(tvmet::norm2(temp) <= _searchR)
	  {
	    _positions.push_back(x);
	  }
	}
    }

    // Initialize functions and derivatives
    nNodes = _positions.size();
    _functions.resize(nNodes);
    _derivatives.resize(nNodes);

    vector<double > r = LME_r(s);
    double res = sqrt(r[0]*r[0] + r[1]*r[1]), det = 0.0; 
    i = 0;
    while(res > _tol && i < _nItMax)
    {
      vector<vector<double > > J = LME_J(s);
      r = LME_r(s);

      det = J[0][0]*J[1][1]-J[1][0]*J[0][1];
    
      _lambda[0] -= (r[0]*J[1][1] - r[1]*J[0][1])/det;
      _lambda[1] -= (r[1]*J[0][0] - r[0]*J[1][0])/det;

      r = LME_r(s);
      i++;
      res = sqrt(r[0]*r[0] + r[1]*r[1]);

      cout << _lambda[0] << " " << _lambda[1] << endl;
    }

    if (isnan(res) != 0)
    {
      Success = false;
      cout << "Error at iteration = " << i << " , position = " << s <<  " res = " << res << endl;
    }


    // compute shape functions and shape functions derivatives
    vector<vector<double > > J = LME_J(s);
    // Compute determinant of J
    double detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];

    vector<vector<double > > Jinv(2,vector<double >(2,0.0));
    Jinv[0][0] = J[1][1]/detJ;
    Jinv[0][1] =-J[0][1]/detJ;
    Jinv[1][0] =-J[1][0]/detJ;
    Jinv[1][1] = J[0][0]/detJ;

    // compute shape functions
    vector<double > temp(2,0.0);
    for(unsigned int i = 0; i < nNodes; i++)
    {
      // compute shape functions
      _functions[i] = LME_pa(s, i);

      temp[0] = s[0]-_positions[i][0];
      temp[1] = s[1]-_positions[i][1];

      _derivatives[i][0] = -LME_pa(s, i)*(temp[0]*Jinv[0][0] + temp[1]*Jinv[0][1]);
      _derivatives[i][1] = -LME_pa(s, i)*(temp[0]*Jinv[1][0] + temp[1]*Jinv[1][1]);
    }
    
    return Success;
   
  };



  double LMEshape2::checkPartitionOfUnity()
  {
    double Sum = 0.0;
    for(unsigned int i = 0; i < _functions.size(); i++)
    {
      Sum += _functions[i];
    }
    return Sum;
  };



  vector<double > LMEshape2::LME_r(const CoordinateArray & s)
  {
    vector<double > r(2, 0.0), temp(2, 0.0);

    for (int i = 0; i < _positions.size(); i++)
    {
      temp[0] = s[0]-_positions[i][0];
      temp[1] = s[1]-_positions[i][1];
   
      r[0] += LME_pa(s, i)*temp[0];
      r[1] += LME_pa(s, i)*temp[1];
    }
    
    return r;
  }



  double LMEshape2::LME_pa(const CoordinateArray & s, int j)
  {
    double Z = 0.0, pa = 0.0;
    vector<double> temp(2, 0.0);

    for (int i = 0; i < _positions.size(); i++)
    {
      temp[0] = s[0]-_positions[i][0];
      temp[1] = s[1]-_positions[i][1];

      Z +=  exp(-_beta*(temp[0]*temp[0] + temp[1]*temp[1]) + _lambda[0]*temp[0] + _lambda[1]*temp[1]);
    }

    temp[0] = s[0]-_positions[j][0];
    temp[1] = s[1]-_positions[j][1];

    pa = exp(-_beta*(temp[0]*temp[0] + temp[1]*temp[1]) + _lambda[0]*temp[0] + _lambda[1]*temp[1])/Z;
    
    return pa;
    
  }



  vector<vector<double > > LMEshape2::LME_J(const CoordinateArray & s)
  {
    vector<vector<double > > J(2, vector<double > (2,0.0));
    vector<double > temp(2, 0.0);
    for (int i = 0; i < _positions.size(); i++)
    {
      temp[0] = s[0]-_positions[i][0];
      temp[1] = s[1]-_positions[i][1];
 
      MatrixSum(J, ScalarProduct(LME_pa(s, i), DyadicProduct(temp,temp) )); 
    }

    MatrixSub(J, DyadicProduct(LME_r(s), LME_r(s) ) ); 

    return J;
  }



  vector<vector<double > > LMEshape2::DyadicProduct(const vector<double > & a, const vector<double > & b)
  {
    vector<vector<double > > D(2, vector<double > (2,0.0));
    D[0][0] = a[0]*b[0];   D[0][1] = a[0]*b[1];
    D[1][0] = a[1]*b[0];   D[1][1] = a[1]*b[1];

    return D;
  }



  vector<vector<double > > LMEshape2::ScalarProduct(const double a, const vector<vector<double > > A)
  {
    vector<vector<double > > D(2, vector<double > (2,0.0));

    D[0][0] = A[0][0]*a;   D[0][1] = A[0][1]*a;  
    D[1][0] = A[1][0]*a;   D[1][1] = A[1][1]*a;
 
    return D;
  }



  void LMEshape2::MatrixSum (vector<vector<double > > & A, const vector<vector<double > > & B)
  {
    A[0][0] += B[0][0];   A[0][1] += B[0][1];
    A[1][0] += B[1][0];   A[1][1] += B[1][1];
   }



  void LMEshape2::MatrixSub (vector<vector<double > > & A, const vector<vector<double > > & B)
  {
    A[0][0] -= B[0][0];   A[0][1] -= B[0][1];
    A[1][0] -= B[1][0];   A[1][1] -= B[1][1];
  }


}
