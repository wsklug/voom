#include "LMEtet.h"
#include <math.h>
namespace voom {

  LMEtet::LMEtet(double beta, double searchR, double tol, int nItMax, vector<int> NodesInd, vector<tvmet::Vector<double, 3> > EvaluationPoints, vector<double > Weights): 
    _beta(beta), _searchR(searchR), _tol(tol), _nItMax(nItMax), _nodesInd(NodesInd), _evaluationPoints(EvaluationPoints), _weights(Weights)
  {
    for (int i = 0; i < 3; i++) {
      _lambda.push_back(0.0);
    }
  };

  bool LMEtet::compute(const CoordinateArray & s, const std::vector<DeformationNode<3>* >& Nodes)
  {
    _lambda[0] = 0.0;
    _lambda[1] = 0.0;
    _lambda[2] = 0.0;
    // Input data are: the point at which to evaluate the shape function, the nodes in the domain.
    bool Success = true;

    // Initialize functions and derivatives
    int _Ns = _nodesInd.size();
    _functions.resize(_Ns);
    _xderivative.resize(_Ns);
    _yderivative.resize(_Ns);
    _zderivative.resize(_Ns);

  

    // Compute lambda
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
    unsigned int i = 0;
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

    if (std::isnan(res) != 0 || i > _nItMax)
    {
      Success = false;
      cout << "i = " << i << " res = " << res << endl;
      return Success;
    }



    // Compute shape functions and shape functions derivatives
    vector<vector<double > > J = LME_J(s, Nodes);
    // Compute determinant of J
    double detJ = J[0][0]*J[1][1]*J[2][2] +  J[0][1]*J[1][2]*J[2][0] +  J[0][2]*J[1][0]*J[2][1] -
                  J[0][0]*J[1][2]*J[2][1] -  J[0][1]*J[1][0]*J[2][2] -  J[0][2]*J[2][0]*J[1][1];    

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

    // Compute shape functions
    vector<double > temp(3,0.0);
    for(unsigned int i = 0; i < _Ns; i++)
    {
      // compute shape functions
      _functions[i] = LME_pa(s, Nodes, i);

      temp[0] = s[0]-Nodes[_nodesInd[i]]->getPosition(0);
      temp[1] = s[1]-Nodes[_nodesInd[i]]->getPosition(1);
      temp[2] = s[2]-Nodes[_nodesInd[i]]->getPosition(2);

      _xderivative[i] = -LME_pa(s, Nodes, i)*(temp[0]*Jinv[0][0] + temp[1]*Jinv[0][1] + temp[2]*Jinv[0][2]);
      _yderivative[i] = -LME_pa(s, Nodes, i)*(temp[0]*Jinv[1][0] + temp[1]*Jinv[1][1] + temp[2]*Jinv[1][2]);
      _zderivative[i] = -LME_pa(s, Nodes, i)*(temp[0]*Jinv[2][0] + temp[1]*Jinv[2][1] + temp[2]*Jinv[2][2]);
    }
    
    return Success;
  };



  void LMEtet::LocalStiffness(const std::vector<DeformationNode<3>* >& Nodes, const double Kmat[3][3][3][3], double *LMEK)
  {
    bool LMEok = false;
    unsigned int LMEfailed_E = 0, LMEfailed_N = 0;

    unsigned int Ns = _nodesInd.size();          // Number of nodes for which we compute the shape function
    double LMEKloc[Ns*3*Ns*3];

    unsigned int Ne = _evaluationPoints.size();  // Number of quadrature points
    unsigned int ne = 0, a = 0, b = 0, i = 0, j = 0;
    // Compute stiffness matrix
    for (i = 0; i < Ns*Ns*9 ;i++) {
      LMEKloc[i] = 0.0;
    }
    
    for (ne = 0; ne < Ne; ne++)
    {
      LMEok = this->compute(_evaluationPoints[ne], Nodes);

      cout << ne << endl;

      /*for (a = 0; a < Ns; a++)
      {
	cout << _functions[a] << "  " << _xderivative[a] << "  " << _yderivative[a] << "  " << _zderivative[a] << endl;
	}*/

      if (LMEok)
      {
	for (a = 0; a < Ns; a++)
	{
	  // cout << endl;
	  for (i = 0; i < 3; i++)
	  {
	    for (b = 0; b < Ns; b++)
	    {
	      for (j = 0; j < 3; j++)
	      {
		LMEKloc[a*3+i + (b*3+j)*3*Ns] +=
		  (Kmat[i][0][j][0]*_xderivative[a]*_xderivative[b]+
		   Kmat[i][1][j][0]*_yderivative[a]*_xderivative[b]+
		   Kmat[i][2][j][0]*_zderivative[a]*_xderivative[b]+
		   Kmat[i][0][j][1]*_xderivative[a]*_yderivative[b]+
		   Kmat[i][1][j][1]*_yderivative[a]*_yderivative[b]+
		   Kmat[i][2][j][1]*_zderivative[a]*_yderivative[b]+
		   Kmat[i][0][j][2]*_xderivative[a]*_zderivative[b]+
		   Kmat[i][1][j][2]*_yderivative[a]*_zderivative[b]+
		   Kmat[i][2][j][2]*_zderivative[a]*_zderivative[b])*_weights[ne];
		if (ne == Ne-1) { cout << LMEKloc[a*3+i + (b*3+j)*3*Ns] << "  "; }
	      }
	    }
	  }
	}
      }
      else
      {
	LMEfailed_E++;
      }
    }

    double *A = new double[3*Ns*3*Ns];
    CoordinateArray s;
    // Compute inversion matrix
    unsigned int ind = 0;
    for (a = 0; a < Ns; a++)
    {
      s(0) = Nodes[_nodesInd[a]]->getPosition(0);
      s(1) = Nodes[_nodesInd[a]]->getPosition(1);
      s(2) = Nodes[_nodesInd[a]]->getPosition(2);
      LMEok = this->compute(s, Nodes);
      cout << endl;
      if (LMEok)
      {
	for (b = 0; b < Ns; b++)
	{
	  A[ind] = _functions[b];
	  A[ind + Ns*3 + 1] = _functions[b];
	  A[ind + Ns*6 + 2] = _functions[b];
	  ind += 3;
	}
	ind += 6*Ns;
      }
      else
      {
	LMEfailed_N++;
      }
    }

    cout << "Evaluation points: shape functions failed = " << LMEfailed_E << " out of " << Ne << endl;
    cout << "Nodes: shape functions failed             = " << LMEfailed_N << " out of " << Ns << endl;
    
    // Invert matrix A
    inverse(A, Ns*3);
    // Transform LME stiffness matrix using A
    multiplyAKAT(A, Ns*3, LMEKloc, LMEK);
    
    
  }



  double LMEtet::checkPartitionOfUnity()
  {
    double Sum = 0.0;
    for(unsigned int i = 0; i < _nodesInd.size(); i++)
    {
      Sum += _functions[i];
    }
    return Sum;
  };



  vector<double > LMEtet::LME_r(const CoordinateArray & s, const std::vector<DeformationNode<3>* >& Nodes)
  {
    vector<double > r(3, 0.0), temp(3, 0.0);

    for (int i = 0; i < _nodesInd.size(); i++)
    {
      temp[0] = s[0]-Nodes[_nodesInd[i]]->getPosition(0);
      temp[1] = s[1]-Nodes[_nodesInd[i]]->getPosition(1);
      temp[2] = s[2]-Nodes[_nodesInd[i]]->getPosition(2);
   
      r[0] += LME_pa(s, Nodes, i)*temp[0];
      r[1] += LME_pa(s, Nodes, i)*temp[1];
      r[2] += LME_pa(s, Nodes, i)*temp[2];
    }
    
    return r;
  }



  double LMEtet::LME_pa(const CoordinateArray & s, const std::vector<DeformationNode<3>* >& Nodes, int j)
  {
    double Z = 0.0, pa = 0.0;
    vector<double> temp(3, 0.0);

    for (int i = 0; i < _nodesInd.size(); i++)
    {
      temp[0] = s[0]-Nodes[_nodesInd[i]]->getPosition(0);
      temp[1] = s[1]-Nodes[_nodesInd[i]]->getPosition(1);
      temp[2] = s[2]-Nodes[_nodesInd[i]]->getPosition(2);

      Z +=  exp(-_beta*(temp[0]*temp[0] + temp[1]*temp[1] + temp[2]*temp[2]) + _lambda[0]*temp[0] + _lambda[1]*temp[1] +  _lambda[2]*temp[2]);
    }

    temp[0] = s[0]-Nodes[_nodesInd[j]]->getPosition(0);
    temp[1] = s[1]-Nodes[_nodesInd[j]]->getPosition(1);
    temp[2] = s[2]-Nodes[_nodesInd[j]]->getPosition(2);

    pa = exp(-_beta*(temp[0]*temp[0] + temp[1]*temp[1] + temp[2]*temp[2]) + _lambda[0]*temp[0] + _lambda[1]*temp[1] + _lambda[2]*temp[2])/Z;
    
    return pa;
    
  }



  vector<vector<double > > LMEtet::LME_J(const CoordinateArray & s, const std::vector<DeformationNode<3>* >& Nodes)
  {
    vector<vector<double > > J(3, vector<double > (3,0.0));
    vector<double > temp(3, 0.0);
    for (int i = 0; i < _nodesInd.size(); i++)
    {
      temp[0] = s[0]-Nodes[_nodesInd[i]]->getPosition(0);
      temp[1] = s[1]-Nodes[_nodesInd[i]]->getPosition(1);
      temp[2] = s[2]-Nodes[_nodesInd[i]]->getPosition(2);

      MatrixSum(J, ScalarProduct(LME_pa(s, Nodes, i), DyadicProduct(temp,temp) )); 
    }

    MatrixSub(J, DyadicProduct(LME_r(s, Nodes), LME_r(s, Nodes) ) ); 

    return J;
  }



  vector<vector<double > > LMEtet::DyadicProduct(const vector<double > & a, const vector<double > & b)
  {
    vector<vector<double > > D(3, vector<double > (3,0.0));
    D[0][0] = a[0]*b[0];   D[0][1] = a[0]*b[1];   D[0][2] = a[0]*b[2];
    D[1][0] = a[1]*b[0];   D[1][1] = a[1]*b[1];   D[1][2] = a[1]*b[2];
    D[2][0] = a[2]*b[0];   D[2][1] = a[2]*b[1];   D[2][2] = a[2]*b[2];

    return D;
  }



  vector<vector<double > > LMEtet::ScalarProduct(const double a, const vector<vector<double > > A)
  {
    vector<vector<double > > D(3, vector<double > (3,0.0));

    D[0][0] = A[0][0]*a;   D[0][1] = A[0][1]*a;   D[0][2] = A[0][2]*a;
    D[1][0] = A[1][0]*a;   D[1][1] = A[1][1]*a;   D[1][2] = A[1][2]*a;
    D[2][0] = A[2][0]*a;   D[2][1] = A[2][1]*a;   D[2][2] = A[2][2]*a;

    return D;

  }



  void LMEtet::MatrixSum (vector<vector<double > > & A, const vector<vector<double > > & B)
  {
    A[0][0] += B[0][0];   A[0][1] += B[0][1];   A[0][2] += B[0][2];
    A[1][0] += B[1][0];   A[1][1] += B[1][1];   A[1][2] += B[1][2];
    A[2][0] += B[2][0];   A[2][1] += B[2][1];   A[2][2] += B[2][2];
  }



  void LMEtet::MatrixSub (vector<vector<double > > & A, const vector<vector<double > > & B)
  {
    A[0][0] -= B[0][0];   A[0][1] -= B[0][1];   A[0][2] -= B[0][2];
    A[1][0] -= B[1][0];   A[1][1] -= B[1][1];   A[1][2] -= B[1][2];
    A[2][0] -= B[2][0];   A[2][1] -= B[2][1];   A[2][2] -= B[2][2];
  }


}
