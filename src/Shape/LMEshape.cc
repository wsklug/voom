#include "LMEshape.h"
#include <math.h>
namespace voom {

  LMEshape::LMEshape(double beta, double searchR, double tol, int nItMax): 
                     _beta(beta), _searchR(searchR), _tol(tol), _nItMax(nItMax)
  {
    for (int i = 0; i < 3; i++) {
      _lambda.push_back(0.0);
    }
  };

  bool LMEshape::compute(const CoordinateArray & s, const std::vector<DeformationNode<3>* >& Nodes, bool FindNeigh, bool f1)
  {
    // Input data are: the point at which to evaluate the shape function, the nodes in the domain.
    bool Success = true;
    int nNodes = Nodes.size(), i = 0;
    
    if (FindNeigh) 
    {
      // Find Neighborhood
      for(i = 0; i < nNodes; i++) 
	{
	  CoordinateArray x(Nodes[i]->getPosition(0), Nodes[i]->getPosition(1), Nodes[i]->getPosition(2));
	  CoordinateArray temp(x-s);
	  if(tvmet::norm2(temp) <= _searchR)
	  {
	    _nodes.push_back(i);
	  }
	}
    }

    // Initialize functions and derivatives
    int _Ns = _nodes.size();
    // cout << _Ns << endl;
    _functions.resize(_Ns);
    _xderivative.resize(_Ns);
    _yderivative.resize(_Ns);
    _zderivative.resize(_Ns);

    

    char uplo = 'U';
    int  n = 3;
    int  nrhs = 1;
    int  lda  = n;
    int  ldb  = n;
    int  info;
    double a[9], b[3];
    int ipiv[3];
    
    vector<double > r = LME_r(s, Nodes);
    double res = r[0]*r[0] + r[1]*r[1] + r[2]*r[2]; 
    i = 0;
    while(res > _tol && i < _nItMax)
    {
      vector<vector<double > > J = LME_J(s, Nodes);
      r = LME_r(s, Nodes);
      
      a[0] = J[0][0];
      a[1] = J[1][0];
      a[2] = J[2][0];
      a[3] = J[0][1];
      a[4] = J[1][1];
      a[5] = J[2][1];
      a[6] = J[0][2];
      a[7] = J[1][2];
      a[8] = J[2][2];

      b[0] = -r[0];
      b[1] = -r[1];
      b[2] = -r[2];

      // Solve for lambda
      // dpotrf_(&uplo, &n, a, &lda, &info);
      // dpotrs_(&uplo, &n, &nrhs, a, &lda, b, &ldb, &info);
      dgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    
      _lambda[0] += b[0];
      _lambda[1] += b[1];
      _lambda[2] += b[2];

      r = LME_r(s, Nodes);
      i++;
      res = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
    }

    if (isnan(res) != 0 || i > _nItMax)
    {
      Success = false;
      cout << "i = " << i << " res = " << res << endl;
    }


    // compute shape functions and shape functions derivatives
    vector<vector<double > > J = LME_J(s, Nodes);
    // Compute determinant of J
    double detJ = J[0][0]*J[1][1]*J[2][2] +  J[0][1]*J[1][2]*J[2][0] +  J[0][2]*J[1][0]*J[2][1] -
                  J[0][0]*J[1][2]*J[2][1] -  J[0][1]*J[1][0]*J[2][2] -  J[0][2]*J[2][0]*J[1][1];
    
    // if (f1 == true) cout << detJ << endl;
    // if (detJ < _tol && f1 == true)
    //{
      /*      cout << s << endl;
      cout << detJ << endl;
for (int m=0; m<3; m++)
      {
      for (int n=0; n<3; n++)
	{
	  cout << J[m][n] << " ";
	}
      cout << endl;
      }
      */
      // return false; // Likely boundary point !
      //}

    
    // cout << "detJ = " << detJ << endl;
    
    // assert(0);
    

    vector<vector<double > > Jinv(3,vector<double >(3,0.0));
    Jinv[0][0] = (J[1][1]*J[2][2] -  J[1][2]*J[2][1])/detJ;
    Jinv[0][1] = (J[0][2]*J[2][1] -  J[0][1]*J[2][2])/detJ;
    Jinv[0][2] = (J[0][1]*J[1][2] -  J[1][1]*J[0][2])/detJ;
    Jinv[1][0] = (J[1][2]*J[2][0] -  J[1][0]*J[2][2])/detJ;
    Jinv[1][1] = (J[0][0]*J[2][2] -  J[0][2]*J[2][0])/detJ;
    Jinv[1][2] = (J[0][2]*J[1][0] -  J[0][0]*J[1][2])/detJ;
    Jinv[2][0] = (J[1][0]*J[2][1] -  J[2][0]*J[1][1])/detJ;
    Jinv[2][1] = (J[0][1]*J[2][0] -  J[0][0]*J[2][1])/detJ;
    Jinv[2][2] = (J[0][0]*J[1][1] -  J[1][0]*J[0][1])/detJ;

    // compute shape functions
    vector<double > temp(3,0.0);
    for(unsigned int i = 0; i < _Ns; i++)
    {
      // compute shape functions
      _functions[i] = LME_pa(s, Nodes, i);

      temp[0] = s[0]-Nodes[_nodes[i]]->getPosition(0);
      temp[1] = s[1]-Nodes[_nodes[i]]->getPosition(1);
      temp[2] = s[2]-Nodes[_nodes[i]]->getPosition(2);

      _xderivative[i] = -LME_pa(s, Nodes, i)*(temp[0]*Jinv[0][0] + temp[1]*Jinv[0][1] + temp[2]*Jinv[0][2]);
      _yderivative[i] = -LME_pa(s, Nodes, i)*(temp[0]*Jinv[1][0] + temp[1]*Jinv[1][1] + temp[2]*Jinv[1][2]);
      _zderivative[i] = -LME_pa(s, Nodes, i)*(temp[0]*Jinv[2][0] + temp[1]*Jinv[2][1] + temp[2]*Jinv[2][2]);
    }
    
    return Success;
   
  };



  double LMEshape::checkPartitionOfUnity()
  {
    double Sum = 0.0;
    for(unsigned int i = 0; i < _nodes.size(); i++)
    {
      Sum += _functions[i];
    }
    return Sum;
  };



  vector<double > LMEshape::LME_r(const CoordinateArray & s, const std::vector<DeformationNode<3>* >& Nodes)
  {
    vector<double > r(3, 0.0), temp(3, 0.0);

    for (int i = 0; i < _nodes.size(); i++)
    {
      temp[0] = s[0]-Nodes[_nodes[i]]->getPosition(0);
      temp[1] = s[1]-Nodes[_nodes[i]]->getPosition(1);
      temp[2] = s[2]-Nodes[_nodes[i]]->getPosition(2);
   
      r[0] += LME_pa(s, Nodes, i)*temp[0];
      r[1] += LME_pa(s, Nodes, i)*temp[1];
      r[2] += LME_pa(s, Nodes, i)*temp[2];
    }
    
    return r;
  }



  double LMEshape::LME_pa(const CoordinateArray & s, const std::vector<DeformationNode<3>* >& Nodes, int j)
  {
    double Z = 0.0, pa = 0.0;
    vector<double> temp(3, 0.0);

    for (int i = 0; i < _nodes.size(); i++)
    {
      temp[0] = s[0]-Nodes[_nodes[i]]->getPosition(0);
      temp[1] = s[1]-Nodes[_nodes[i]]->getPosition(1);
      temp[2] = s[2]-Nodes[_nodes[i]]->getPosition(2);

      Z +=  exp(-_beta*(temp[0]*temp[0] + temp[1]*temp[1] + temp[2]*temp[2]) + _lambda[0]*temp[0] + _lambda[1]*temp[1] +  _lambda[2]*temp[2]);
    }

    temp[0] = s[0]-Nodes[_nodes[j]]->getPosition(0);
    temp[1] = s[1]-Nodes[_nodes[j]]->getPosition(1);
    temp[2] = s[2]-Nodes[_nodes[j]]->getPosition(2);

    pa = exp(-_beta*(temp[0]*temp[0] + temp[1]*temp[1] + temp[2]*temp[2]) + _lambda[0]*temp[0] + _lambda[1]*temp[1] + _lambda[2]*temp[2])/Z;
    
    return pa;
    
  }



  vector<vector<double > > LMEshape::LME_J(const CoordinateArray & s, const std::vector<DeformationNode<3>* >& Nodes)
  {
    vector<vector<double > > J(3, vector<double > (3,0.0));
    vector<double > temp(3, 0.0);
    for (int i = 0; i < _nodes.size(); i++)
    {
      temp[0] = s[0]-Nodes[i]->getPosition(0);
      temp[1] = s[1]-Nodes[i]->getPosition(1);
      temp[2] = s[2]-Nodes[i]->getPosition(2);

      MatrixSum(J, ScalarProduct(LME_pa(s, Nodes, i), DyadicProduct(temp,temp) )); 
    }

    MatrixSub(J, DyadicProduct(LME_r(s, Nodes), LME_r(s, Nodes) ) ); 

    return J;
  }



  vector<vector<double > > LMEshape::DyadicProduct(const vector<double > & a, const vector<double > & b)
  {
    vector<vector<double > > D(3, vector<double > (3,0.0));
    D[0][0] = a[0]*b[0];   D[0][1] = a[0]*b[1];   D[0][2] = a[0]*b[2];
    D[1][0] = a[1]*b[0];   D[1][1] = a[1]*b[1];   D[1][2] = a[1]*b[2];
    D[2][0] = a[2]*b[0];   D[2][1] = a[2]*b[1];   D[2][2] = a[2]*b[2];

    return D;
  }



  vector<vector<double > > LMEshape::ScalarProduct(const double a, const vector<vector<double > > A)
  {
    vector<vector<double > > D(3, vector<double > (3,0.0));

    D[0][0] = A[0][0]*a;   D[0][1] = A[0][1]*a;   D[0][2] = A[0][2]*a;
    D[1][0] = A[1][0]*a;   D[1][1] = A[1][1]*a;   D[1][2] = A[1][2]*a;
    D[2][0] = A[2][0]*a;   D[2][1] = A[2][1]*a;   D[2][2] = A[2][2]*a;

    return D;

  }



  void LMEshape::MatrixSum (vector<vector<double > > & A, const vector<vector<double > > & B)
  {
    A[0][0] += B[0][0];   A[0][1] += B[0][1];   A[0][2] += B[0][2];
    A[1][0] += B[1][0];   A[1][1] += B[1][1];   A[1][2] += B[1][2];
    A[2][0] += B[2][0];   A[2][1] += B[2][1];   A[2][2] += B[2][2];
  }



  void LMEshape::MatrixSub (vector<vector<double > > & A, const vector<vector<double > > & B)
  {
    A[0][0] -= B[0][0];   A[0][1] -= B[0][1];   A[0][2] -= B[0][2];
    A[1][0] -= B[1][0];   A[1][1] -= B[1][1];   A[1][2] -= B[1][2];
    A[2][0] -= B[2][0];   A[2][1] -= B[2][1];   A[2][2] -= B[2][2];
  }


}
