namespace voom
{
  template<class Material_t>
  void InteractionElement<Material_t>::compute(bool f0, bool f1, bool f2)

  {
    std::vector<Vector3D> _x1; //quad points of e1
    std::vector<Vector3D> _x2; //quad points of e2

    std::vector<double> _W1; //weight of each quad point of e1
    std::vector<double> _W2; //weight of each quad point of e2

    //shape function at each quad point for each node, e1
    std::vector< blitz::Array< double, 1> > _N1; 
    //shape function at each quad point for each node, e2
    std::vector< blitz::Array< double, 1> > _N2;
    
    int i=0, j=0;
    Vector3D rhat; //unit distance vector 
    double r; //distance
	
    //store all the quad points, weights and shape functions for e1
    for(QuadPointIterator p=_e1->_quadPoints.begin();p!=_e1->_quadPoints.end();p++){

      Vector3D x(0.0), X(0.0);
      tvmet::Vector< Vector3D, 2 > a;
      tvmet::Matrix< Vector3D, 2, 2 > aPartials;
      Vector3D zero(0);
      a = zero, zero;
      aPartials = zero, zero, zero, zero;

      const LoopShellShape::FunctionArray & N
	= p->shape.functions();
      const LoopShellShape::DerivativeArray & DN
	= p->shape.derivatives();
      const LoopShellShape::SecondDerivativeArray & DDN
	= p->shape.secondDerivatives();

      for (int b = 0; b < _e1->_nodes.size(); b++){
	const DeformationNode<3>::Point & xb = _e1->_nodes[b]->point();
	const DeformationNode<3>::Point & Xb = _e1->_nodes[b]->position();
	x 	       +=   N(b)     * xb;
	
	X 	       +=   N(b)     * Xb;
	a(0)           +=  DN(b,0)   * Xb;
	a(1)           +=  DN(b,1)   * Xb;
	aPartials(0,0) += DDN(b,0,0) * Xb;
	aPartials(0,1) += DDN(b,0,1) * Xb;
	aPartials(1,1) += DDN(b,1,1) * Xb;

      }
           
      aPartials(1,0) = aPartials(0,1);//by symmetry
      
      Material_t& material = p->material;
	// compute shell geometry
      ShellGeometry refgeometry( a, aPartials );
      
	// store the ref geometry in shell geometry class
      material.setRefGeometry(refgeometry);
      
      const double refMetric = ( material.refShellGeometry()).metric();
      const double weight =  refMetric * p->weight;
      
      
      _x1.push_back(x);     //ith quad point position
      _N1.push_back(N); //shape functions of nodes at ith quad point
      _W1.push_back(weight);//weight of ith quad point
      
      i++;
    }

    //repeat for e2
    for(QuadPointIterator p=_e2->_quadPoints.begin();p!=_e2->_quadPoints.end();p++){

      Vector3D x(0.0), X(0.0);
      tvmet::Vector< Vector3D, 2 > a;
      tvmet::Matrix< Vector3D, 2, 2 > aPartials;
      Vector3D zero(0);
      a = zero, zero;
      aPartials = zero, zero, zero, zero;

      const LoopShellShape::FunctionArray & N
	= p->shape.functions();
      const LoopShellShape::DerivativeArray & DN
	= p->shape.derivatives();
      const LoopShellShape::SecondDerivativeArray & DDN
	= p->shape.secondDerivatives();
      
      for (int b = 0; b < _e2->_nodes.size(); b++){
	const DeformationNode<3>::Point & xb = _e2->_nodes[b]->point();
	const DeformationNode<3>::Point & Xb = _e2->_nodes[b]->position();
	x              +=   N(b)     * xb;
	
	X 	       +=   N(b)     * Xb;
	a(0)           +=  DN(b,0)   * Xb;
	a(1)           +=  DN(b,1)   * Xb;
	aPartials(0,0) += DDN(b,0,0) * Xb;
	aPartials(0,1) += DDN(b,0,1) * Xb;
	aPartials(1,1) += DDN(b,1,1) * Xb;

      }
      aPartials(1,0) = aPartials(0,1);//by symmetry
      
      Material_t& material = p->material;
      // compute shell geometry
      ShellGeometry refgeometry( a, aPartials );
	  			
      // store the ref geometry in shell geometry class
      material.setRefGeometry(refgeometry);
      
      const double refMetric = ( material.refShellGeometry()).metric();
      const double weight =  refMetric * p->weight;
      
      _x2.push_back(x);     //jth quad point position
      _N2.push_back(N); //shape functions of nodes at jth quad point
      _W2.push_back(weight);//weight of jth quad point

      j++;
    }
    
    
    if( f0 ){
      _energy = 0.0;
      _adhesionEnergy = 0.0;
    }
    
    if( f1 ){
      Vector3D zero(0);
      _internalForce1 = zero;
      _internalForce2 = zero;
    }

    if( f0 ){
      for(int i=0; i<_x1.size(); i++)
	for(int j=0; j<_x2.size(); j++){
	  r=tvmet::norm2(_x1[i]-_x2[j]);
	  if (r<=_rc)
	    _adhesionEnergy += (_k*( r-_rZero )*( r-_rZero ) - _epsilon)*_W1[i]*_W2[j];
	}

      _energy = _adhesionEnergy;      
    }
    

    if( f1 ){
      for(int i=0; i<_x1.size(); i++)
	for(int j=0; j<_x2.size(); j++){
	  
	  r=tvmet::norm2(_x2[j]-_x1[i]);
	  rhat=(_x2[j]-_x1[i])/r;
	  
	  if (r<=_rc){
	    for(int a = 0; a < _e1->_nodes.size(); a++){
	      _internalForce1(a) += - _N1[i](a) *2.0*_k*(r-_rZero) * rhat * _W1[i]*_W2[j]; 	    
	    }
	  
	    for(int b = 0; b < _e2->_nodes.size(); b++){
	      _internalForce2(b) += _N2[j](b) *2.0*_k*(r-_rZero) * rhat * _W1[i]*_W2[j]; 
	    }    
	  }

	}
	

	int a=0;
	for(NodeIterator na=_e1->_nodes.begin();  na!=_e1->_nodes.end(); na++, a++)
	  for(int i=0; i<3; i++) {
	    double f_ia = _internalForce1(a)(i); 
	    (*na)->addForce( i, f_ia );
	  }
      
	a=0;
	for(NodeIterator na=_e2->_nodes.begin();  na!=_e2->_nodes.end(); na++, a++)
	  for(int j=0; j<3; j++) {
	    double f_ja = _internalForce2(a)(j);
 
	    (*na)->addForce( j, f_ja );
	  }
      
    }
  }
    


  
}// namespace voom
