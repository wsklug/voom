#include "MFShape.h"
namespace voom {
  // commonly used cubic B-spline kernel
  void MFShape::compute(const CoordinateArray & s, const std::vector<DeformationNode<3>* >& nodes, const std::vector<double>& SuppSize, const std::vector<int> & list, double radius){
    //given array of all nodes and it's number and support size etc., it calculates the shape functions and its derivatives
    //at point s given the radius of sphere around s (required for calculating the derivatives using SNNI)
    tvmet::Matrix<double,4,4> M(0);
    int nNodes = nodes.size();
    bool flag;
    
    //number of neighbours
    int nNeigh = 0;
    for(int i=0;i<list.size();i++){
      CoordinateArray x(nodes[list[i]]->getPosition(0),nodes[list[i]]->getPosition(1),nodes[list[i]]->getPosition(2));
      CoordinateArray temp(x-s);
      double dist = tvmet::norm2(temp);
      if(dist<=SuppSize[list[i]]+radius){
	nNeigh++;
	_nodes.push_back(list[i]);
      }
    }
    
    //calculate the shape functions
    flag=compute_function(s,nodes,SuppSize,nNeigh,_functions);
    if(flag==false) std::cout<<"Matrix singular calculating shape function at the node"<<std::endl;
    
    //calculating the derivative of shape functions
    for(int i=0;i<nNeigh;i++){
      _xderivative.push_back(0.);
      _yderivative.push_back(0.);
      _zderivative.push_back(0.);}
    //take the vertices of dodecahedron as the integration points (20 in number)
    std::vector<CoordinateArray> x(20);
    double phi=0.5+sqrt(5./4.);
    x[0]=1.,1.,1.;
    x[1]=1.,1.,-1.;
    x[2]=1.,-1.,1.;
    x[3]=1.,-1.,-1.;
    x[4]=-1.,1.,1.;
    x[5]=-1.,1.,-1.;
    x[6]=-1.,-1.,1.;
    x[7]=-1.,-1.,-1.;

    x[8]=0.,1./phi,phi;
    x[9]=0.,1./phi,-phi;
    x[10]=0.,-1./phi,phi;
    x[11]=0.,-1./phi,-phi;

    x[12]=1./phi,phi,0.;
    x[13]=1./phi,-phi,0.;
    x[14]=-1./phi,phi,0.;
    x[15]=-1./phi,-phi,0.;

    x[16]=phi,0.,1./phi;
    x[17]=-phi,0.,1./phi;
    x[18]=phi,0.,-1./phi;
    x[19]=-phi,0.,-1./phi;
    /*std::vector<CoordinateArray> x(6);
    x[0]=1.,0.,0.;
    x[1]=-1.,0.,0.;
    x[2]=0.,1.,0.;
    x[3]=0.,-1.,0.;
    x[4]=0.,0.,1.;
    x[5]=0.,0.,-1.;*/
    //double smooth=1.;
    for(int i=0;i<x.size();i++){
      //normalize x
      x[i]=x[i]/tvmet::norm2(x[i]);
      //integration point
      CoordinateArray p(s+x[i]*_smooth*radius);
      FunctionContainer temp;
      temp.reserve(nNeigh);
      flag=compute_function(p,nodes,SuppSize,nNeigh,temp);
      if(flag==false) std::cout<<"Matrix singular calculating derivative of shape function and the radius is "<<radius<<std::endl;
      for(int j=0;j<nNeigh;j++){
        _xderivative[j] += temp[j]*x[i](0);
        _yderivative[j] += temp[j]*x[i](1);
        _zderivative[j] += temp[j]*x[i](2);}
      }
   
    for(int i=0;i<nNeigh;i++){
      _xderivative[i] *= (3./_smooth/radius/x.size());
      _yderivative[i] *= (3./_smooth/radius/x.size());
      _zderivative[i] *= (3./_smooth/radius/x.size());}
  }

  //! uses singular kernel
  void MFShape::compute_sing(const CoordinateArray & s, const std::vector<DeformationNode<3>* >& nodes, const std::vector<double>& SuppSize, double radius, double sing_order){
    //given array of all nodes and it's number and support size etc., it calculates the shape functions and its derivatives
    //at point s given the radius of sphere around s (required for calculating the derivatives using SNNI)
    int nNodes = nodes.size();
    bool flag;
    
    //number of neighbours
    int nNeigh = 0;
    for(int i=0;i<nNodes;i++){
      CoordinateArray x(nodes[i]->getPosition(0),nodes[i]->getPosition(1),nodes[i]->getPosition(2));
      CoordinateArray temp(x-s);
      double dist = tvmet::norm2(temp);
      if(dist<=SuppSize[i]+radius){
	nNeigh++;
	_nodes.push_back(i);
      }
    }
    //can't calculate the shape functions numerically due to the singularity, but they are trivialy kronecker delta
    //flag=compute_function_sing(s,nodes,SuppSize,sing_order,nNeigh,_functions);
    //if(flag==false) std::cout<<"Matrix singular calculating shape function at the node"<<std::endl;
    for(int i=0;i<nNeigh;i++){
      CoordinateArray x(nodes[_nodes[i]]->getPosition(0),nodes[_nodes[i]]->getPosition(1),nodes[_nodes[i]]->getPosition(2));
      CoordinateArray temp(x-s);
      double dist = tvmet::norm2(temp);
      if(dist<1e-10)
        _functions.push_back(1.);
      else
        _functions.push_back(0.);
    }

    
    //calculating the derivative of shape functions
    for(int i=0;i<nNeigh;i++){
      _xderivative.push_back(0.);
      _yderivative.push_back(0.);
      _zderivative.push_back(0.);}

    //take the vertices of dodecahedron as the integration points (20 in number)
    std::vector<CoordinateArray> x(20);
    double phi=0.5+sqrt(5./4.);
    x[0]=1.,1.,1.;
    x[1]=1.,1.,-1.;
    x[2]=1.,-1.,1.;
    x[3]=1.,-1.,-1.;
    x[4]=-1.,1.,1.;
    x[5]=-1.,1.,-1.;
    x[6]=-1.,-1.,1.;
    x[7]=-1.,-1.,-1.;

    x[8]=0.,1./phi,phi;
    x[9]=0.,1./phi,-phi;
    x[10]=0.,-1./phi,phi;
    x[11]=0.,-1./phi,-phi;

    x[12]=1./phi,phi,0.;
    x[13]=1./phi,-phi,0.;
    x[14]=-1./phi,phi,0.;
    x[15]=-1./phi,-phi,0.;

    x[16]=phi,0.,1./phi;
    x[17]=-phi,0.,1./phi;
    x[18]=phi,0.,-1./phi;
    x[19]=-phi,0.,-1./phi;
    /*std::vector<CoordinateArray> x(6);
    x[0]=1.,0.,0.;
    x[1]=-1.,0.,0.;
    x[2]=0.,1.,0.;
    x[3]=0.,-1.,0.;
    x[4]=0.,0.,1.;
    x[5]=0.,0.,-1.;*/
    //normalize x
    for(int i=0;i<x.size();i++)
      x[i]=x[i]/tvmet::norm2(x[i]);
    int numQuadrants=8;
    const double rCentroid=radius*3.*sqrt(3.)/8.;
    std::vector<CoordinateArray> QuadCenter(numQuadrants);
  //  QuadCenter[0]=0.,0.,0.;
    QuadCenter[0]=1.,1.,1.;
    QuadCenter[1]=1.,1.,-1.;
    QuadCenter[2]=1.,-1.,1.;
    QuadCenter[3]=1.,-1.,-1.;
    QuadCenter[4]=-1.,1.,1.;
    QuadCenter[5]=-1.,1.,-1.;
    QuadCenter[6]=-1.,-1.,1.;
    QuadCenter[7]=-1.,-1.,-1.;
    double radiusQuad=radius/2.; //radius/pow(numQuadrants,1./3.);
    //normalize QuadCenter 
    for(int i=0;i<QuadCenter.size();i++)
      QuadCenter[i]=QuadCenter[i]/tvmet::norm2(QuadCenter[i]);

    //double smooth=1.0;
    for(int numQ=0; numQ<numQuadrants;numQ++){
      // Get center of the quadrant
      CoordinateArray C(rCentroid*QuadCenter[numQ]);
      for(int i=0;i<x.size();i++){
        //integration point
        CoordinateArray p(s+C+x[i]*_smooth*radiusQuad);
        FunctionContainer temp;
        temp.reserve(nNeigh);
        flag=compute_function_sing(p,nodes,SuppSize,sing_order,nNeigh,temp);
        if(flag==false) std::cout<<"Matrix singular calculating derivative of shape function and the radius is "<<radius<<std::endl;
        for(int j=0;j<nNeigh;j++){
          _xderivative[j] += temp[j]*x[i](0);
          _yderivative[j] += temp[j]*x[i](1);
          _zderivative[j] += temp[j]*x[i](2);}
        }
      }
   
    for(int i=0;i<nNeigh;i++){
      _xderivative[i] *= (3./_smooth/radiusQuad/x.size()/numQuadrants);
      _yderivative[i] *= (3./_smooth/radiusQuad/x.size()/numQuadrants);
      _zderivative[i] *= (3./_smooth/radiusQuad/x.size()/numQuadrants);}
  }

  //! uses modified kernel
  void MFShape::compute_mod(const CoordinateArray & s, const std::vector<DeformationNode<3>* > & nodes, const std::vector<double> & SuppSize, const std::vector<double> & SuppSizeHat, const std::vector<int> & list, double radius){
    //given array of all nodes and it's number and support size etc., it calculates the shape functions and its derivatives with modified kernel having kronecker delta property (therefore needs two support sizes, one for primitive and one for enrichment function
    //at point s given the radius of sphere around s (required for calculating the derivatives using SNNI)

    int nNodes = nodes.size();
    bool flag;
    //std::cout<<"Size of the list is "<<list.size()<<std::endl; 
    //number of neighbours
    int nNeigh = 0;
    for(int i=0;i<list.size();i++){
      CoordinateArray x(nodes[list[i]]->getPosition(0),nodes[list[i]]->getPosition(1),nodes[list[i]]->getPosition(2));
      CoordinateArray temp(x-s);
      double dist = tvmet::norm2(temp);
      if(dist<=SuppSize[list[i]]+radius){
	nNeigh++;
	_nodes.push_back(list[i]);
      }
    }
    /*for(int i=0;i<nNodes;i++){
      CoordinateArray x(nodes[i]->getPosition(0),nodes[i]->getPosition(1),nodes[i]->getPosition(2));
      CoordinateArray temp(x-s);
      double dist = tvmet::norm2(temp);
      if(dist<=SuppSize[i]+radius){
	nNeigh++;
	_nodes.push_back(i);
      }
    }*/
    //calculate the shape functions
    flag=compute_function_mod(s,nodes,SuppSize,SuppSizeHat,nNeigh,_functions);
    if(flag==false) std::cout<<"Matrix singular calculating shape function at the node"<<std::endl;
    
    //calculating the derivative of shape functions
    for(int i=0;i<nNeigh;i++){
      _xderivative.push_back(0.);
      _yderivative.push_back(0.);
      _zderivative.push_back(0.);}

    //take the vertices of dodecahedron as the integration points (20 in number)
    std::vector<CoordinateArray> x(20);
    double phi=0.5+sqrt(5./4.);
    x[0]=1.,1.,1.;
    x[1]=1.,1.,-1.;
    x[2]=1.,-1.,1.;
    x[3]=1.,-1.,-1.;
    x[4]=-1.,1.,1.;
    x[5]=-1.,1.,-1.;
    x[6]=-1.,-1.,1.;
    x[7]=-1.,-1.,-1.;

    x[8]=0.,1./phi,phi;
    x[9]=0.,1./phi,-phi;
    x[10]=0.,-1./phi,phi;
    x[11]=0.,-1./phi,-phi;

    x[12]=1./phi,phi,0.;
    x[13]=1./phi,-phi,0.;
    x[14]=-1./phi,phi,0.;
    x[15]=-1./phi,-phi,0.;

    x[16]=phi,0.,1./phi;
    x[17]=-phi,0.,1./phi;
    x[18]=phi,0.,-1./phi;
    x[19]=-phi,0.,-1./phi;
    /*std::vector<CoordinateArray> x(6);
    x[0]=1.,0.,0.;
    x[1]=-1.,0.,0.;
    x[2]=0.,1.,0.;
    x[3]=0.,-1.,0.;
    x[4]=0.,0.,1.;
    x[5]=0.,0.,-1.;*/
    //normalize x
    for(int i=0;i<x.size();i++)
      x[i]=x[i]/tvmet::norm2(x[i]);
    int numQuadrants=1;
    const double rCentroid=radius*3.*sqrt(3.)/8.;
    std::vector<CoordinateArray> QuadCenter(numQuadrants);
    QuadCenter[0]=0.,0.,0.;
    /*QuadCenter[0]=1.,1.,1.;
    QuadCenter[1]=1.,1.,-1.;
    QuadCenter[2]=1.,-1.,1.;
    QuadCenter[3]=1.,-1.,-1.;
    QuadCenter[4]=-1.,1.,1.;
    QuadCenter[5]=-1.,1.,-1.;
    QuadCenter[6]=-1.,-1.,1.;
    QuadCenter[7]=-1.,-1.,-1.;*/
    double radiusQuad=radius;///2.; //radius/pow(numQuadrants,1./3.);
    //normalize QuadCenter 
  //  for(int i=0;i<QuadCenter.size();i++)
  //    QuadCenter[i]=QuadCenter[i]/tvmet::norm2(QuadCenter[i]);

    //double smooth=1.0;
    for(int numQ=0; numQ<numQuadrants;numQ++){
      // Get center of the quadrant
      CoordinateArray C(rCentroid*QuadCenter[numQ]);
      for(int i=0;i<x.size();i++){
        //integration point
        CoordinateArray p(s+C+x[i]*_smooth*radiusQuad);
        FunctionContainer temp;
        temp.reserve(nNeigh);
        flag=compute_function_mod(p,nodes,SuppSize,SuppSizeHat,nNeigh,temp);
        if(flag==false) std::cout<<"Matrix singular calculating derivative of shape function and the radius is "<<radius<<std::endl;
        for(int j=0;j<nNeigh;j++){
          _xderivative[j] += temp[j]*x[i](0);
          _yderivative[j] += temp[j]*x[i](1);
          _zderivative[j] += temp[j]*x[i](2);}
        }
      }
   
    for(int i=0;i<nNeigh;i++){
      _xderivative[i] *= (3./_smooth/radiusQuad/x.size()/numQuadrants);
      _yderivative[i] *= (3./_smooth/radiusQuad/x.size()/numQuadrants);
      _zderivative[i] *= (3./_smooth/radiusQuad/x.size()/numQuadrants);}
  }

//modified to use 2nd order basis
  bool MFShape::compute_function_sing(const CoordinateArray & s, const std::vector<DeformationNode<3>* >& nodes, const std::vector<double>& SuppSize, double sing_order, int nNeigh, FunctionContainer & function){
    tvmet::Matrix<double,4,4> M(0);
    int tnNeigh=0;
    for(int i=0;i<nNeigh;i++){
      CoordinateArray x(nodes[_nodes[i]]->getPosition(0),nodes[_nodes[i]]->getPosition(1),nodes[_nodes[i]]->getPosition(2));
      CoordinateArray temp(x-s);
      double dist = tvmet::norm2(temp);  
      double phi = cal_phi_sing(dist/SuppSize[_nodes[i]],sing_order);  if(phi>0) tnNeigh++;
      tvmet::Matrix<double,4,1> Ht;
      Ht = 1,temp(0),temp(1),temp(2);//,temp(0)*temp(0),temp(1)*temp(1),temp(2)*temp(2),temp(0)*temp(1),temp(1)*temp(2),temp(2)*temp(0);
      M = M + Ht*tvmet::trans(Ht)*phi;
    }
   //invert the moment matrix with LAPACK (10*10 matrix) 
    tvmet::Matrix<double,4,4> invM(0);

    //LAPACK related stuff
    double A[4][4];
    int n = 4, info;
    int lda = n;
    int ipiv[n];
    int lwork = n;
    double work[lwork];
    for(int i = 0; i < n; i++)
      for(int j = 0; j < n; j++) A[i][j] = M(i,j);

    // LU decomposition first
    dgetrf_(&n, &n, &A[0][0],&lda, ipiv, &info);
    // Inversion then
    dgetri_(&n, &A[0][0], &lda, ipiv, work, &lwork, &info);
    bool flag=true;

    if (info != 0) {
      flag=false;
      std::cerr << "DGETRI: Matrix Inversion Failed\n";
    }
    for(int i = 0; i < n; i++)
      for(int j = 0; j < n; j++) invM(i,j) = A[i][j];

    //bool flag=MFShape::inverse(M,invM);
    //if(info != 0) std::cout<<"Matrix is singular and number of neighbours: "<<tnNeigh<<std::endl;
    
    for(int i=0;i<nNeigh;i++){
      CoordinateArray X(nodes[_nodes[i]]->getPosition(0),nodes[_nodes[i]]->getPosition(1),nodes[_nodes[i]]->getPosition(2));
      CoordinateArray temp(X-s);
      double dist = tvmet::norm2(temp);
      double phi=cal_phi_sing(dist/SuppSize[_nodes[i]],sing_order);
      tvmet::Matrix<double,4,1> H;
      H = 1,temp(0),temp(1),temp(2);//,temp(0)*temp(0),temp(1)*temp(1),temp(2)*temp(2),temp(0)*temp(1),temp(1)*temp(2),temp(2)*temp(0);
      tvmet::Matrix<double,1,4> H0(0);
      H0=1.0,0.0,0.0,0.0;//,0.,0.,0.,0.,0.,0.;
      tvmet::Matrix<double,1,1> temp2;
      temp2 = H0*invM*H*phi;
      function.push_back(temp2(0,0));
    }
    return flag;
  }
  
  bool MFShape::compute_function(const CoordinateArray & s, const std::vector<DeformationNode<3>* >& nodes, const std::vector<double>& SuppSize, int nNeigh, FunctionContainer & function){
    tvmet::Matrix<double,4,4> M(0);
    int tnNeigh=0;
    for(int i=0;i<nNeigh;i++){
      CoordinateArray x(nodes[_nodes[i]]->getPosition(0),nodes[_nodes[i]]->getPosition(1),nodes[_nodes[i]]->getPosition(2));
      CoordinateArray temp(x-s);
      double dist = tvmet::norm2(temp);  
      double phi = cal_phi(dist/SuppSize[_nodes[i]]);  if(phi>0) tnNeigh++;
      tvmet::Matrix<double,4,1> Ht;
      Ht = 1,temp(0),temp(1),temp(2);
      M = M + Ht*tvmet::trans(Ht)*phi;
    }
    
    tvmet::Matrix<double,4,4> invM(0);
    bool flag=MFShape::inverse(M,invM);
    if(flag==false) std::cout<<"Matrix is singular and number of neighbours: "<<tnNeigh<<std::endl;
    
    for(int i=0;i<nNeigh;i++){
      CoordinateArray X(nodes[_nodes[i]]->getPosition(0),nodes[_nodes[i]]->getPosition(1),nodes[_nodes[i]]->getPosition(2));
      CoordinateArray temp(X-s);
      double dist = tvmet::norm2(temp);
      double phi=cal_phi(dist/SuppSize[_nodes[i]]);
      tvmet::Matrix<double,4,1> H;
      H=1.0,temp(0),temp(1),temp(2);
      tvmet::Matrix<double,1,4> H0;
      H0=1.0,0.0,0.0,0.0;
      tvmet::Matrix<double,1,1> temp2;
      temp2 = H0*invM*H*phi;
      function.push_back(temp2(0,0));
    }
    return flag;
  }
  bool MFShape::compute_function_mod(const CoordinateArray & s, const std::vector<DeformationNode<3>* >& nodes, const std::vector<double>& SuppSize, const std::vector<double>& SuppSizeHat, int nNeigh, FunctionContainer & function){
    tvmet::Matrix<double,4,4> M(0);
    tvmet::Matrix<double,4,1> Fhat(0);
    int tnNeigh=0;
    for(int i=0;i<nNeigh;i++){
      CoordinateArray x(nodes[_nodes[i]]->getPosition(0),nodes[_nodes[i]]->getPosition(1),nodes[_nodes[i]]->getPosition(2));
      CoordinateArray temp(x-s);
      double dist = tvmet::norm2(temp);  
      double phi = cal_phi(dist/SuppSize[_nodes[i]]);  if(phi>0) tnNeigh++;
      double phihat = cal_phi(dist/SuppSizeHat[_nodes[i]]);
      tvmet::Matrix<double,4,1> Ht;
      Ht = 1,temp(0),temp(1),temp(2);
      M = M + Ht*tvmet::trans(Ht)*phi;
      Fhat = Fhat + Ht*phihat;
    }
    
    tvmet::Matrix<double,4,4> invM(0);
    bool flag=MFShape::inverse(M,invM);
    if(flag==false) std::cout<<"Matrix is singular and number of neighbours: "<<tnNeigh<<std::endl;
    
    for(int i=0;i<nNeigh;i++){
      CoordinateArray X(nodes[_nodes[i]]->getPosition(0),nodes[_nodes[i]]->getPosition(1),nodes[_nodes[i]]->getPosition(2));
      CoordinateArray temp(X-s);
      double dist = tvmet::norm2(temp);
      double phi=cal_phi(dist/SuppSize[_nodes[i]]);
      double phihat = cal_phi(dist/SuppSizeHat[_nodes[i]]);
      tvmet::Matrix<double,1,4> H;
      H=1.0,temp(0),temp(1),temp(2);
      tvmet::Matrix<double,4,1> H0;
      H0=1.0,0.0,0.0,0.0;
      tvmet::Matrix<double,1,1> temp2;
      temp2 = H*invM*(H0-Fhat)*phi;
      function.push_back(temp2(0,0)+phihat);
    }
    return flag;
  }
  
  float MFShape::partition_of_unity(const CoordinateArray & s, const std::vector<DeformationNode<3>* >& nodes, const std::vector<double>& SuppSize){
    //tvmet::Matrix<double,4,4> M(0);
    double M=0.; //the moment matrix is just a scalar because only a constant basis is used
    
    int nNodes = nodes.size();
    //number of neighbours
    int nNeigh = 0;
    for(int i=0;i<nNodes;i++){
      CoordinateArray x(nodes[i]->getPosition(0),nodes[i]->getPosition(1),nodes[i]->getPosition(2));
      CoordinateArray temp(x-s);
      double dist = tvmet::norm2(temp);
      double phi = cal_phi(dist/SuppSize[i]);  
      if(phi>0){
	nNeigh++;
	_nodes.push_back(i);
      }
      M += phi;
    }
    
    if(M==0) return 0;
    double sum_shapefn=0.;
    
    for(int i=0;i<nNeigh;i++){
      CoordinateArray X(nodes[_nodes[i]]->getPosition(0),nodes[_nodes[i]]->getPosition(1),nodes[_nodes[i]]->getPosition(2));
      CoordinateArray temp(X-s);
      double dist = tvmet::norm2(temp);
      double phi=cal_phi(dist/SuppSize[_nodes[i]]);
      sum_shapefn += phi/M;
    }
    return sum_shapefn;
  }
  
  // Parition of Unity Check
  bool MFShape::checkPartitionOfUnity(const CoordinateArray& s, const std::vector<DeformationNode<3>* >& nodes) {
    double sumShape = 0;
    const double TOLER = 1E-6;
    for(int i = 0; i < _functions.size(); i++)
      sumShape += _functions[i];

    double x = 0., y = 0., z=0.;
    for(int i = 0; i < _functions.size(); i++) {
	x += _functions[i] * nodes[_nodes[i]]->getPosition(0);
	y += _functions[i] * nodes[_nodes[i]]->getPosition(1);
	z += _functions[i] * nodes[_nodes[i]]->getPosition(2);
    }
    CoordinateArray Pos1(x,y,z);
    CoordinateArray error1(Pos1 - s);
   
    double xsumDer = 0., ysumDer = 0., zsumDer=0.;
    for(int i = 0; i < _xderivative.size(); i++) {
	xsumDer += _xderivative[i];
	ysumDer += _yderivative[i];
        zsumDer += _zderivative[i];
      //  std::cout<<_xderivative[i]<<"\t"<<_yderivative[i]<<"\t"<<_zderivative[i]<<std::endl;
    }
    CoordinateArray Pos2(xsumDer, ysumDer,zsumDer);
    CoordinateArray error2(Pos2);
    x = y = z = 0.;
    for(int i = 0; i < _xderivative.size(); i++) {
	x += _xderivative[i] * nodes[_nodes[i]]->getPosition(0);
	y += _yderivative[i] * nodes[_nodes[i]]->getPosition(1);
	z += _zderivative[i] * nodes[_nodes[i]]->getPosition(2);
    }
    CoordinateArray Pos3(x, y, z);
    CoordinateArray One(1.,1.,1.);
    CoordinateArray error3(Pos3 - One);
    
    if ( fabs(sumShape - 1.) < TOLER && tvmet::norm2(error1) < TOLER 
	&& tvmet::norm2(error2) < TOLER && tvmet::norm2(error3) < TOLER ) 
	return true;
    else {
	std::cout << "Check 1: " << fabs(sumShape - 1.) << "\n";
	std::cout << "Check 2: " << tvmet::norm2(error1) << "\n";
	std::cout << "Check 3: " << tvmet::norm2(error2) << "\n";
	std::cout << "Check 4: " << tvmet::norm2(error3) << "\n";
	return false;
    }
  }
  
  bool MFShape::inverse(tvmet::Matrix<double,4,4> M, tvmet::Matrix<double,4,4> & invM){
    //determinant of M
    bool flag=true;
    double det 
      = M(0,0)*M(1,1)*M(2,2)*M(3,3) 
      - M(0,0)*M(1,1)*M(2,3)*M(3,2) 
      - M(0,0)*M(1,2)*M(2,1)*M(3,3) 
      + M(0,0)*M(1,2)*M(2,3)*M(3,1) 
      + M(0,0)*M(1,3)*M(2,1)*M(3,2) 
      - M(0,0)*M(1,3)*M(2,2)*M(3,1) 
      - M(0,1)*M(1,0)*M(2,2)*M(3,3) 
      + M(0,1)*M(1,0)*M(2,3)*M(3,2) 
      + M(0,1)*M(1,2)*M(2,0)*M(3,3) 
      - M(0,1)*M(1,2)*M(2,3)*M(3,0) 
      - M(0,1)*M(1,3)*M(2,0)*M(3,2) 
      + M(0,1)*M(1,3)*M(2,2)*M(3,0) 
      + M(0,2)*M(1,0)*M(2,1)*M(3,3) 
      - M(0,2)*M(1,0)*M(2,3)*M(3,1) 
      - M(0,2)*M(1,1)*M(2,0)*M(3,3) 
      + M(0,2)*M(1,1)*M(2,3)*M(3,0) 
      + M(0,2)*M(1,3)*M(2,0)*M(3,1) 
      - M(0,2)*M(1,3)*M(2,1)*M(3,0) 
      - M(0,3)*M(1,0)*M(2,1)*M(3,2) 
      + M(0,3)*M(1,0)*M(2,2)*M(3,1) 
      + M(0,3)*M(1,1)*M(2,0)*M(3,2) 
      - M(0,3)*M(1,1)*M(2,2)*M(3,0) 
      - M(0,3)*M(1,2)*M(2,0)*M(3,1) 
      + M(0,3)*M(1,2)*M(2,1)*M(3,0);
    if(det<1e-9) {//std::cout<<"The matrix is singular"<<std::endl; 
      flag=false;}
    //closed form inverse of M
    invM = M(1,1)* M(2,2)* M(3,3) - M(1,1) *M(2,3)* M(3,2) - M(1,2)* M(2,1)* M(3,3) + M(1,2)* M(2,3)* M(3,1) + M(1,3)* M(2,1)* M(3,2) - M(1,3)* M(2,2)* M(3,1), M(0,1)* M(2,3)* M(3,2) - M(0,1) *M(2,2)* M(3,3) + M(0,2)* M(2,1)* M(3,3) - M(0,2)* M(2,3) *M(3,1) - M(0,3)* M(2,1)* M(3,2) + M(0,3) *M(2,2)* M(3,1), M(0,1)* M(1,2)* M(3,3) - M(0,1)* M(1,3) *M(3,2) - M(0,2)* M(1,1)* M(3,3) + M(0,2) *M(1,3) *M(3,1) + M(0,3)* M(1,1) *M(3,2) - M(0,3)* M(1,2)* M(3,1), M(0,1)* M(1,3)* M(2,2) - M(0,1)* M(1,2) *M(2,3) + M(0,2) *M(1,1)* M(2,3) - M(0,2)* M(1,3) *M(2,1) - M(0,3)* M(1,1) *M(2,2) + M(0,3)* M(1,2) *M(2,1),         
      M(1,0) *M(2,3)* M(3,2) - M(1,0) *M(2,2)* M(3,3) + M(1,2)* M(2,0)* M(3,3) - M(1,2)* M(2,3)*M(3,0) - M(1,3)* M(2,0)* M(3,2) + M(1,3)* M(2,2)* M(3,0), M(0,0) *M(2,2)* M(3,3) - M(0,0)* M(2,3)* M(3,2) - M(0,2)* M(2,0) *M(3,3) + M(0,2)* M(2,3)* M(3,0) + M(0,3)* M(2,0)* M(3,2) - M(0,3)* M(2,2)* M(3,0), M(0,0)* M(1,3) *M(3,2) - M(0,0)* M(1,2)* M(3,3) + M(0,2)* M(1,0)* M(3,3) - M(0,2) *M(1,3) *M(3,0) - M(0,3)* M(1,0) *M(3,2) + M(0,3)* M(1,2)* M(3,0), M(0,0)* M(1,2)* M(2,3) - M(0,0) *M(1,3)* M(2,2) - M(0,2) *M(1,0) *M(2,3) + M(0,2)* M(1,3)* M(2,0) + M(0,3)* M(1,0)* M(2,2) - M(0,3) *M(1,2) *M(2,0),
      M(1,0) *M(2,1)* M(3,3) - M(1,0)* M(2,3)* M(3,1) - M(1,1)* M(2,0)* M(3,3) + M(1,1)* M(2,3)* M(3,0) + M(1,3)* M(2,0) *M(3,1) - M(1,3)* M(2,1)* M(3,0), M(0,0) *M(2,3)* M(3,1) - M(0,0)* M(2,1) *M(3,3) + M(0,1)* M(2,0)* M(3,3) - M(0,1) *M(2,3) *M(3,0) - M(0,3)* M(2,0)* M(3,1) + M(0,3) *M(2,1) *M(3,0), M(0,0)* M(1,1) *M(3,3) - M(0,0) *M(1,3) *M(3,1) - M(0,1) *M(1,0) *M(3,3) + M(0,1) *M(1,3)* M(3,0) + M(0,3) *M(1,0) *M(3,1) - M(0,3) *M(1,1) *M(3,0), M(0,0) *M(1,3)* M(2,1) - M(0,0) *M(1,1) *M(2,3) + M(0,1) *M(1,0) *M(2,3) - M(0,1) *M(1,3) *M(2,0) - M(0,3) *M(1,0) *M(2,1) + M(0,3) *M(1,1)* M(2,0),
      M(1,0) *M(2,2) *M(3,1) - M(1,0) *M(2,1)* M(3,2) + M(1,1)* M(2,0) *M(3,2) - M(1,1) *M(2,2)* M(3,0) - M(1,2)* M(2,0) *M(3,1) + M(1,2)* M(2,1) *M(3,0), M(0,0)* M(2,1)* M(3,2) - M(0,0) *M(2,2)* M(3,1) - M(0,1)* M(2,0) *M(3,2) + M(0,1)* M(2,2)* M(3,0) + M(0,2)* M(2,0)* M(3,1) - M(0,2) *M(2,1)* M(3,0), M(0,0) *M(1,2) *M(3,1) - M(0,0)* M(1,1)* M(3,2) + M(0,1)* M(1,0)* M(3,2) - M(0,1)* M(1,2)* M(3,0) - M(0,2)* M(1,0)* M(3,1) + M(0,2)* M(1,1) *M(3,0), M(0,0) *M(1,1) *M(2,2) - M(0,0)* M(1,2)* M(2,1) - M(0,1)* M(1,0)* M(2,2) + M(0,1) *M(1,2)* M(2,0) + M(0,2) *M(1,0) *M(2,1) - M(0,2)* M(1,1)* M(2,0) ;
    invM = invM/det;
    return flag;
  }
  
}
