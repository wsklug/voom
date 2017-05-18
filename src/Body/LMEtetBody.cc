#include <string>
#include <fstream>
#include <blitz/array-impl.h>
#include "VoomMath.h"



namespace voom
{
  LMEtetBody::LMEtetBody(vector< DeformationNode<3>* > & defNodes, vector<vector<unsigned int > > & connT, 
			double beta, double search_radius, unsigned int quad_order, double tol, unsigned int nItMax, 
			double rho, double E, double nu): 
    _defNodes(defNodes), _connT(connT), _beta(beta), _search_radius(search_radius), _tol(tol), _nItMax(nItMax), _rho(rho)
  { 
    // Determine evaluation points
    vector<double > nqp(4, 0.25);
    switch (quad_order)
    {
    case 1:
      _NatQP.push_back(nqp);
      _weightQP.push_back(1.0);
      break;
    case 2:
      nqp[0] = 0.585410196624969;
      nqp[1] = 0.138196601125011;
      nqp[2] = 0.138196601125011;
      nqp[3] = 0.138196601125011;
        _NatQP.push_back(nqp);
	_weightQP.push_back(0.25);
      nqp[0] = 0.138196601125011;
      nqp[1] = 0.585410196624969;
        _NatQP.push_back(nqp);
	_weightQP.push_back(0.25);
      nqp[1] = 0.138196601125011;
      nqp[2] = 0.585410196624969;
        _NatQP.push_back(nqp);
	_weightQP.push_back(0.25);
      nqp[2] = 0.138196601125011;
      nqp[3] = 0.585410196624969;
        _NatQP.push_back(nqp);
	_weightQP.push_back(0.25);
      break;
    case 3:
      nqp[0] = 0.25;
      nqp[1] = 0.25;
      nqp[2] = 0.25;
      nqp[3] = 0.25;
        _NatQP.push_back(nqp);
	_weightQP.push_back(-0.8);
      
      nqp[0] = 0.5;
      nqp[1] = 0.166666666666667;
      nqp[2] = 0.166666666666667;
      nqp[3] = 0.166666666666667;
        _NatQP.push_back(nqp);
	_weightQP.push_back(0.45);
      nqp[0] = 0.166666666666667;
      nqp[1] = 0.5;
        _NatQP.push_back(nqp);
	_weightQP.push_back(0.45);
      nqp[1] = 0.166666666666667;
      nqp[2] = 0.5;
        _NatQP.push_back(nqp);
	_weightQP.push_back(0.45);
      nqp[2] = 0.166666666666667;
      nqp[3] = 0.5;
        _NatQP.push_back(nqp);
	_weightQP.push_back(0.45);
      break;
    case 4:
      nqp[0] = 0.25;
      nqp[1] = 0.25;
      nqp[2] = 0.25;
      nqp[3] = 0.25;
        _NatQP.push_back(nqp);
	_weightQP.push_back(-0.013155555555556);

      nqp[0] = 0.785714285714286;
      nqp[1] = 0.071428571428571;
      nqp[2] = 0.071428571428571;
      nqp[3] = 0.071428571428571;
        _NatQP.push_back(nqp);
	_weightQP.push_back(0.007622222222222);
      nqp[0] = 0.071428571428571;
      nqp[1] = 0.785714285714286;
        _NatQP.push_back(nqp);
	_weightQP.push_back(0.007622222222222);
      nqp[1] = 0.071428571428571;
      nqp[2] = 0.785714285714286;
        _NatQP.push_back(nqp);
	_weightQP.push_back(0.007622222222222);
      nqp[2] = 0.071428571428571;
      nqp[3] = 0.785714285714286;
        _NatQP.push_back(nqp);
	_weightQP.push_back(0.007622222222222);

      nqp[0] = 0.399403576166799;
      nqp[1] = 0.399403576166799;
      nqp[2] = 0.100596423833201;
      nqp[3] = 0.100596423833201;
        _NatQP.push_back(nqp);
	_weightQP.push_back(0.024888888888889);
      nqp[0] = 0.399403576166799;
      nqp[1] = 0.100596423833201;
      nqp[2] = 0.399403576166799;
      nqp[3] = 0.100596423833201;
        _NatQP.push_back(nqp);
	_weightQP.push_back(0.024888888888889);
      nqp[0] = 0.399403576166799;
      nqp[1] = 0.100596423833201;
      nqp[2] = 0.100596423833201;
      nqp[3] = 0.399403576166799;
        _NatQP.push_back(nqp);
	_weightQP.push_back(0.024888888888889);
      nqp[0] = 0.100596423833201;
      nqp[1] = 0.399403576166799;
      nqp[2] = 0.399403576166799;
      nqp[3] = 0.100596423833201;
        _NatQP.push_back(nqp);
	_weightQP.push_back(0.024888888888889);
      nqp[0] = 0.100596423833201;
      nqp[1] = 0.399403576166799;
      nqp[2] = 0.100596423833201;
      nqp[3] = 0.399403576166799;
        _NatQP.push_back(nqp);
	_weightQP.push_back(0.024888888888889);
      nqp[0] = 0.100596423833201;
      nqp[1] = 0.100596423833201;
      nqp[2] = 0.399403576166799; 
      nqp[3] = 0.399403576166799;
        _NatQP.push_back(nqp);
	_weightQP.push_back(0.024888888888889);
      break;

    case 5:
      nqp[0] = 0.25;
      nqp[1] = 0.25;
      nqp[2] = 0.25;
      nqp[3] = 0.25;
        _NatQP.push_back(nqp);
	_weightQP.push_back(0.030283678097089);

      nqp[0] = 0.0;
      nqp[1] = 0.333333333333333;
      nqp[2] = 0.333333333333333;
      nqp[3] = 0.333333333333333;
        _NatQP.push_back(nqp);
	_weightQP.push_back(0.006026785714286);
      nqp[0] = 0.333333333333333;
      nqp[1] = 0.0;
        _NatQP.push_back(nqp);
	_weightQP.push_back(0.006026785714286);
      nqp[1] = 0.333333333333333;
      nqp[2] = 0.0;
        _NatQP.push_back(nqp);
	_weightQP.push_back(0.006026785714286);
      nqp[2] = 0.333333333333333;
      nqp[3] = 0.0;
        _NatQP.push_back(nqp);
	_weightQP.push_back(0.006026785714286);

      nqp[0] = 0.727272727272727;
      nqp[1] = 0.090909090909091;
      nqp[2] = 0.090909090909091;
      nqp[3] = 0.090909090909091;
        _NatQP.push_back(nqp);
	_weightQP.push_back(0.011645249086029);
      nqp[0] = 0.090909090909091;
      nqp[1] = 0.727272727272727;
        _NatQP.push_back(nqp);
	_weightQP.push_back(0.011645249086029);
      nqp[1] = 0.090909090909091;
      nqp[2] = 0.727272727272727;
        _NatQP.push_back(nqp);
	_weightQP.push_back(0.011645249086029);
      nqp[2] = 0.090909090909091;
      nqp[3] = 0.727272727272727;
        _NatQP.push_back(nqp);
	_weightQP.push_back(0.011645249086029);

      nqp[0] = 0.066550153573664;
      nqp[1] = 0.066550153573664;
      nqp[2] = 0.433449846426336;
      nqp[3] = 0.433449846426336;
        _NatQP.push_back(nqp);
	_weightQP.push_back(0.010949141561386);
      nqp[0] = 0.066550153573664;
      nqp[1] = 0.433449846426336;
      nqp[2] = 0.066550153573664;
      nqp[3] = 0.433449846426336;
        _NatQP.push_back(nqp);
	_weightQP.push_back(0.010949141561386);
      nqp[0] = 0.066550153573664;
      nqp[1] = 0.433449846426336;
      nqp[2] = 0.433449846426336;
      nqp[3] = 0.066550153573664;
        _NatQP.push_back(nqp);
	_weightQP.push_back(0.010949141561386);
      nqp[0] = 0.433449846426336;
      nqp[1] = 0.066550153573664;
      nqp[2] = 0.066550153573664;
      nqp[3] = 0.433449846426336;
        _NatQP.push_back(nqp);
	_weightQP.push_back(0.010949141561386);
      nqp[0] = 0.433449846426336;
      nqp[1] = 0.066550153573664;
      nqp[2] = 0.433449846426336;
      nqp[3] = 0.066550153573664;
        _NatQP.push_back(nqp);
	_weightQP.push_back(0.010949141561386);
      nqp[0] = 0.433449846426336;
      nqp[1] = 0.433449846426336;
      nqp[2] = 0.066550153573664;
      nqp[3] = 0.066550153573664;
        _NatQP.push_back(nqp);
	_weightQP.push_back(0.010949141561386);
    break;  
    default:
      cout << "Quad order = " << quad_order << " not implemented." << endl;
      exit(1);
    };

    // Compute Neo-Hookean stiffness matrix for F = I
    Tensor3D ID(0.0);
    ID(0,0) = 1.0;
    ID(1,1) = 1.0;
    ID(2,2) = 1.0;
    unsigned int i = 0, J = 0, k = 0, L = 0;
    double lambda = E*nu/((1.0+nu)*(1.0-2.0*nu));
    double mu = 0.5*E/(1+nu);
    
    for (i = 0; i<3; i++) {
      for (J = 0; J<3; J++) {
	for (k = 0; k<3; k++) {
	  for (L = 0; L<3; L++) {
	    _Kmat[i][J][k][L] = mu*ID(J,k)*ID(L,i) + lambda*ID(J,i)*ID(L,k) + mu*ID(i,k)*ID(J,L); 
	  }
	}
      }
    }

  }; 

  void LMEtetBody::computeNormalModes(double* K, double *eval, unsigned int ModesNumber, bool EgvFlag)
  {
    // Prepare space for global stiffness matrix
    int N = _defNodes.size()*3, info = 0;
    int LWORK = 3*N;
    double WORK[LWORK]; 

    char jobz = 'N';
    char uplo = 'U';
    if (EgvFlag)
    {
      jobz = 'V';
    }

    // Loop over all elements
    unsigned int el = 0, i = 0, j = 0, NumEvP = 0;
    Vector3D A(0.0), B(0.0), C(0.0), D(0.0), AB(0.0), BC(0.0), CD(0.0), TetBar(0.0), EvP(0.0);
    double TetVol =0.0;
    
    for (el=0; el < _connT.size(); el++)
    {
      // Compute volume of tetrahedral element
      A(0) = _defNodes[_connT[el][0]]->getPosition(0);
      A(1) = _defNodes[_connT[el][0]]->getPosition(1);
      A(2) = _defNodes[_connT[el][0]]->getPosition(2);

      B(0) = _defNodes[_connT[el][1]]->getPosition(0);
      B(1) = _defNodes[_connT[el][1]]->getPosition(1);
      B(2) = _defNodes[_connT[el][1]]->getPosition(2);

      C(0) = _defNodes[_connT[el][2]]->getPosition(0);
      C(1) = _defNodes[_connT[el][2]]->getPosition(1);
      C(2) = _defNodes[_connT[el][2]]->getPosition(2);

      D(0) = _defNodes[_connT[el][3]]->getPosition(0);
      D(1) = _defNodes[_connT[el][3]]->getPosition(1);
      D(2) = _defNodes[_connT[el][3]]->getPosition(2);

      AB = A-B;
      BC = B-C;
      CD = C-D;

      TetVol = fabs(AB(0)*BC(1)*CD(2) + AB(1)*BC(2)*CD(0) + AB(2)*BC(0)*CD(1) - AB(2)*BC(1)*CD(0) - AB(0)*BC(2)*CD(1) - AB(1)*BC(0)*CD(2))/6.0;

      

      // Compute evaluation points positions and respective weights
      NumEvP = _weightQP.size();
      vector<double > Weights(NumEvP, 0.0);
      vector<Vector3D > EvaluationPoints;
      for (i = 0; i < NumEvP; i++)
      {
	Weights[i] = _weightQP[i]*TetVol;
	EvP = A*_NatQP[i][0] + B*_NatQP[i][1] + C*_NatQP[i][2] + D*_NatQP[i][3];
	EvaluationPoints.push_back(EvP);
      }
      
      // Compute nodes at a distance less than search_R from the tetrahedra element barycenter
      vector<int > NodesInd;
      TetBar = A*0.25 + B*0.25 + C*0.25 + D*0.25;

      for (i = 0; i < _defNodes.size(); i++)
      {
	if (tvmet::norm2(TetBar - _defNodes[i]->position()) < _search_radius)
	{
	  NodesInd.push_back(i);
	}
      }
  
      unsigned int NIS = NodesInd.size();
      LMEtet LMEshape(_beta, _search_radius, _tol, _nItMax, NodesInd, EvaluationPoints, Weights);
      double* LMEK = new double[3*NIS*3*NIS];
      LMEshape.LocalStiffness(_defNodes, _Kmat, LMEK);
      
      // Assemble global matrix from local ones
      for (i = 0; i < NIS; i++) {
	for (j = 0; j < NIS; j++) {
	  K[NodesInd[i]*3   + NodesInd[j]*3*N] += LMEK[i*3   + j*9*NIS];
	  K[NodesInd[i]*3+1 + NodesInd[j]*3*N] += LMEK[i*3+1 + j*9*NIS];
	  K[NodesInd[i]*3+2 + NodesInd[j]*3*N] += LMEK[i*3+2 + j*9*NIS];

	  K[NodesInd[i]*3   + (NodesInd[j]*3+1)*N] += LMEK[i*3   + (j*9+3)*NIS];
	  K[NodesInd[i]*3+1 + (NodesInd[j]*3+1)*N] += LMEK[i*3+1 + (j*9+3)*NIS];
	  K[NodesInd[i]*3+2 + (NodesInd[j]*3+1)*N] += LMEK[i*3+2 + (j*9+3)*NIS];

	  K[NodesInd[i]*3   + (NodesInd[j]*3+2)*N] += LMEK[i*3   + (j*9+6)*NIS];
	  K[NodesInd[i]*3+1 + (NodesInd[j]*3+2)*N] += LMEK[i*3+1 + (j*9+6)*NIS];
	  K[NodesInd[i]*3+2 + (NodesInd[j]*3+2)*N] += LMEK[i*3+2 + (j*9+6)*NIS];
	}
      }
      delete LMEK;

    } // End loop over elements

    /*for (i=0; i<12; i++)
    {
      cout << endl;
      for (j=0; j<12; j++)
      {
	cout << K[i*12+j] << " ";
      }
    }
    */
    // Compute eigenvalues and eigevectors
    /* N = 3;
    LWORK = 3*N;
    double WORKt[LWORK]; 
    double Kt[9];
    
    //Kt[0] = 6.0; Kt[1] =-1.0;  Kt[2] = 2.0;
    //Kt[3] =-1.0; Kt[4] = 22.0; Kt[5] = 15.0;
    //Kt[6] = 2.0; Kt[7] = 15.0; Kt[8] = 17.0;
    
    // Upper triangle using columnwise storage
    Kt[0] = 6.0; Kt[1] = 0.0;  Kt[2] = 0.0;
    Kt[3] =-1.0; Kt[4] = 22.0; Kt[5] = 0.0;
    Kt[6] = 2.0; Kt[7] = 15.0; Kt[8] = 17.0;
    */
    dsyev_(&jobz, &uplo, &N, Kt, &N, eval, WORKt, &LWORK, &info);

  };
  


} // namespace voom
