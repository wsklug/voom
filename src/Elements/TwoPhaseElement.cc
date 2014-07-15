// -*- C++ -*-
//----------------------------------------------------------------------
//
//                   William S. Klug, Feng Feng
//                University of California Los Angeles
//                   (C) 2004 All Rights Reserved
//
//----------------------------------------------------------------------


////////////////////////////////////////////////////////////////////
// LAPACK subroutine for computing eigenvalues and eigenvectors
//
// define prototype of LAPACK functions
extern "C" void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda,
		       double *w, double *work, int *lwork, int *info);

//#define _DEBUG_


namespace voom
{
  using std::endl;
  using std::cout;

  template<class Material_t>
  void TwoPhaseElement<Material_t>::compute(bool f0, bool f1, bool f2)
  {
    //
    // initialize things
    //
    // if applied volume constraint, _volume is needed for computing
    // if applied area constraint, _area is needed for computing
    // both energy and forces
    _volume = 0.0;
    _area = 0.0;
    _areaOne =0.0;
    // _trC = 0.0;

    

    //          3
    //        a/ \c
    //        /	  \
    //       1-----2
    //          b
    // Compute lengths of sides from current positions
    //
    const XCNode<3>::Point & x1 = _nodes[0]->point();
    const XCNode<3>::Point & x2 = _nodes[1]->point();
    const XCNode<3>::Point & x3 = _nodes[2]->point(); 
    double aa 
      = (x1[0]-x3[0])*(x1[0]-x3[0])
      + (x1[1]-x3[1])*(x1[1]-x3[1])
      + (x1[2]-x3[2])*(x1[2]-x3[2]);
    aa = pow(aa,0.5);

    double bb 
      = (x1[0]-x2[0])*(x1[0]-x2[0])
      + (x1[1]-x2[1])*(x1[1]-x2[1])
      + (x1[2]-x2[2])*(x1[2]-x2[2]);
    bb = pow(bb,0.5);

    double cc 
      = (x2[0]-x3[0])*(x2[0]-x3[0])
      + (x2[1]-x3[1])*(x2[1]-x3[1])
      + (x2[2]-x3[2])*(x2[2]-x3[2]);
    cc = pow(cc,0.5);


    //
    // Compute s, A, h from ref positions
    //
    const XCNode<3>::Point & X1 = _nodes[0]->position();
    const XCNode<3>::Point & X2 = _nodes[1]->position();
    const XCNode<3>::Point & X3 = _nodes[2]->position(); 
    double AA 
      = (X1[0]-X3[0])*(X1[0]-X3[0])
      + (X1[1]-X3[1])*(X1[1]-X3[1])
      + (X1[2]-X3[2])*(X1[2]-X3[2]);
    AA = pow(AA,0.5);

    double BB 
      = (X1[0]-X2[0])*(X1[0]-X2[0])
      + (X1[1]-X2[1])*(X1[1]-X2[1])
      + (X1[2]-X2[2])*(X1[2]-X2[2]);
    BB = pow(BB,0.5);

    double CC 
      = (X2[0]-X3[0])*(X2[0]-X3[0])
      + (X2[1]-X3[1])*(X2[1]-X3[1])
      + (X2[2]-X3[2])*(X2[2]-X3[2]);
    CC = pow(CC,0.5);

    double s = (AA+BB+CC)/2.0;

    double A = s*(s-AA)*(s-BB)*(s-CC);
    A = pow(A,0.5);

    double fourOverRoot3 = 4.0/sqrt(3.0);
    double h = pow(fourOverRoot3*A, 0.5);
		
    double shapeCri=1.0/2.0*((AA-h)*(AA-h)+(BB-h)*(BB-h)+(CC-h)*(CC-h))/h/h;

    if( f0 ) {
      _energy = 0.0;
      _localConstraintEnergy = 0.0;
      _strainEnergy = 0.0;
      _constraintEnergy = 0.0;
      _chemicalConstraintEnergy = 0.0;
      _work = 0.0;
      _springEnergy = 0.0;
    }

    if( f1 ) {
      Vector3D zero(0.0); Vector4D zero4(0.0);
      _internalForce = zero4;
      _chemicalTensionForce = zero4;
      _pressureForce = zero;
      _tensionForce = zero;
    }
		
    // loop for every quadrature point
    for(QuadPointIterator p=_quadPoints.begin(); 
	p!=_quadPoints.end(); p++){
      //
      // compute Shell geometry
      //
			
      // compute position, basis vector, and derivatives of basis vector
      Vector3D x(0.0);
      tvmet::Vector< Vector3D, 2 > a;
      tvmet::Matrix< Vector3D, 2, 2 > aPartials;
      Vector3D zero(0);
      a = zero, zero;
      //ra = zero, zero;
      aPartials = zero, zero, zero, zero;
      //raPartials = zero, zero, zero, zero;			

      //concentration and gradient
      double C  = 0.0;
      tvmet::Vector< double, 2 > dC(0.0);

      const LoopShellShape::FunctionArray & N 
	= p->shape.functions();
      const LoopShellShape::DerivativeArray & DN 
	= p->shape.derivatives();
      const LoopShellShape::SecondDerivativeArray & DDN 
	= p->shape.secondDerivatives();

      for (int b = 0; b < _nodes.size(); b++){
	const XCNode<3>::Point & xb = _nodes[b]->point();
	const double & cb = _nodes[b]->getPoint(3);

	x 	       +=   N(b)     * xb;
	a(0)           +=  DN(b,0)   * xb;
	a(1)           +=  DN(b,1)   * xb;
	aPartials(0,0) += DDN(b,0,0) * xb;
	aPartials(0,1) += DDN(b,0,1) * xb;
	aPartials(1,1) += DDN(b,1,1) * xb;

	C     += N(b)    * cb;
	dC(0) += DN(b,0) * cb;
	dC(1) += DN(b,1) * cb; 

      }

      aPartials(1,0)  = aPartials(0,1); // by symmetry
			
      Material_t& material = p->material;
      // compute shell geometry
      ShellGeometry geometry( a, aPartials );
      const Vector3D& d = geometry.d();
      const tvmet::Vector< Vector3D, 2 >& aDual = geometry.aDual();
      const tvmet::Vector<Vector3D, 2>& dPartials = geometry.dPartials();
      const tvmet::Matrix<Vector3D, 2, 2>& dualPartials = geometry.aDualPartials();
			
      //double material.getConcentration(c);
      //double material.getConDerivative(cd); 

      // store the deformed geometry in shell geometry class
      material.setGeometry(geometry);

      //put concentration information to material
      material.setConcentration(C, dC); 

      // compute strain energy, stress and moment resultants
      material.updateState(f0, f1, f2); 
			
      const double metric = /*0.5 */ geometry.metric();
      // const double refMetric = 0.5 *refgeometry.metric();
      const double refMetric = /*0.5 */( material.refShellGeometry()).metric();
      const double jacobian = metric/refMetric;
      const double weight =  metric * p->weight;

      // compute area for the area constaint energy
      _area += weight;

      // compute arar if the first component
      _areaOne += C*weight;

      // compute volume for the volume constraint energy
      _volume +=  dot(d,x) * weight / 3.0;

      
      // compute energy
      if ( f0 ){

	// compute strain energy 
	_strainEnergy += material.energyDensity() * weight;
       
        _localConstraintEnergy += material.constraintEnergyDensity()*weight;
       
      }


      // compute forces
      if ( f1 ) {
	//
	// get stress and moment resultants
	const tvmet::Vector< Vector3D, 3 >& sr = material.stressResultants();
	const tvmet::Vector< Vector3D, 2 >& mr = material.momentResultants();
	const tvmet::Vector< Vector3D, 3 >& cr = material.concentrationStressResultants();
	const double&                      chr = material.chemicalResultants();
	const tvmet::Vector< double, 2 >& chgr = material.chemicalGradientResultants();

	double pressure = _pressureNode->point();
	
	double tension = _tensionNode->point();

	double chemicalTension = _chemicalTensionNode->point();

	// loop for all nodes to compute forces 
	for (int a=0; a<_nodes.size(); a++) {
	  Vector3D _tempInternalForce(0.0);
	  Vector3D _tempChemicalTensionForce(0.0);

	  // compute internal forces

	  // calculate the gradient of the derivatives of the director
	  // w.r.t curvilinear coords
	  for ( int alpha = 0; alpha < 2; alpha++){
	    // Stress & Concentration Resultant part
	    Vector3D ftmp;
	    ftmp = (sr(alpha) + cr(alpha)) *  DN(a,alpha) * weight;
 	    _tempInternalForce  += ftmp;
	    // Moment Resultant part
	    for(int beta=0; beta<2; beta++) {
	      ftmp = -  dot(mr(alpha),aDual[beta])*
			( DDN(a,alpha,beta)*d + DN(a,beta)*dPartials(alpha) )
		-  dot(mr(alpha),dualPartials(beta,alpha))*DN(a,beta)*d;
	      ftmp *= weight;
	      _tempInternalForce += ftmp;
	    }
	  } 

	  for(int i=0;i<3;i++){
	    _internalForce(a)(i) += _tempInternalForce(i);
	  }

	  _internalForce(a)(3) += chr * N(a)*weight;
	  for(int beta=0; beta<2; beta++){
	    _internalForce(a)(3) += chgr(beta)*DN(a,beta)*weight;
	  } 

          // global area constraint
	  _tensionForce(a) += 
	    tension * (DN(a,0) * aDual[0] + DN(a,1) * aDual[1]) * weight;

	  // compute pressure/volume constraint forces 
// 	  _pressureForce(a) -= pressure * d * N(a) * weight;
	  _pressureForce(a) -= 
	    pressure*( d * N(a) 
		       + dot(x,d)*( aDual[0]*DN(a,0)+
				    aDual[1]*DN(a,1) )
		       - ( dot(x,aDual[0])*DN(a,0) + 
			   dot(x,aDual[1])*DN(a,1)  )*d
		       )*weight/3.0;


	    _tempChemicalTensionForce +=
	      chemicalTension * (DN(a,0) * aDual[0] + DN(a,1) * aDual[1]) * C * weight;

	  for(int i=0;i<3;i++){
	    _chemicalTensionForce(a)(i) += _tempChemicalTensionForce(i);
	  }

	  _chemicalTensionForce(a)(3) += chemicalTension * N(a) * weight;


	  
	} // end nodes loop
	
      } // end force calcs

        
    } // end quadrature loop

    QuadPointIterator p=_quadPoints.begin();
    Material_t& material = p->material;
    //get spring constant from material
    _kSpring = material.getKSpring();

    if(f0) {
      if (shapeCri>1.0)
        _springEnergy = 0.5*_kSpring*((aa-h)*(aa-h)+(bb-h)*(bb-h)+(cc-h)*(cc-h));
      else
	_springEnergy = 0.5*_kSpring*((aa-AA)*(aa-AA)+(bb-BB)*(bb-BB)+(cc-CC)*(cc-CC));

      _strainEnergy += _springEnergy;
      _energy = _strainEnergy;
      _localConstraintEnergy += _springEnergy;

    }

    if(f1) {
      int a=0, ia=0;
      for(NodeIterator na=_nodes.begin();  na!=_nodes.end(); na++, a++){
	for(int i=0; i<3; i++, ia++) {
	  double f_ia = _internalForce(a)(i) 
	    	      + _pressureForce(a)(i) 
	    	      + _tensionForce(a)(i)
	              + _chemicalTensionForce(a)(i);
	  (*na)->addForce( i, f_ia );
	}

	double f_ia = _internalForce(a)(3) + _chemicalTensionForce(a)(3);
	(*na)->addForce( 3, f_ia );

      }

      typename XCNode<3>::Point f;
      NodeIterator na= _nodes.begin();

      if (shapeCri > 1.0){
	f = _kSpring*(aa-h)*(x1-x3)/aa;
	(*na)->updateForce(f);
	f = _kSpring*(bb-h)*(x1-x2)/bb;
	(*na)->updateForce(f);

	na++;
	f = _kSpring*(bb-h)*(x2-x1)/bb;
	(*na)->updateForce(f);
	f = _kSpring*(cc-h)*(x2-x3)/cc;
	(*na)->updateForce(f);

	na++;
	f = _kSpring*(aa-h)*(x3-x1)/aa;
	(*na)->updateForce(f);
	f = _kSpring*(cc-h)*(x3-x2)/cc;
	(*na)->updateForce(f);
	}

      else{
	f = _kSpring*(aa-AA)*(x1-x3)/aa;
	(*na)->updateForce(f);
	f = _kSpring*(bb-BB)*(x1-x2)/bb;
	(*na)->updateForce(f);

	na++;
	f = _kSpring*(bb-BB)*(x2-x1)/bb;
	(*na)->updateForce(f);
	f = _kSpring*(cc-CC)*(x2-x3)/cc;
	(*na)->updateForce(f);

	na++;
	f = _kSpring*(aa-AA)*(x3-x1)/aa;
	(*na)->updateForce(f);
	f = _kSpring*(cc-CC)*(x3-x2)/cc;
	(*na)->updateForce(f);
	}
    }


  }


  template<class Material_t>
  void TwoPhaseElement<Material_t>::updateRefConfiguration() {
    for(QuadPointIterator p=_quadPoints.begin(); 
	p!=_quadPoints.end(); p++){
      //
      // compute Shell geometry
      //
			
      // compute position, basis vector, and derivatives of basis vector
      Vector3D X(0.0);
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

      for (int b = 0; b < _nodes.size(); b++){
	const DeformationNode<3>::Point & Xb = _nodes[b]->position();
	X 	       +=   N(b)     * Xb;
	a(0)           +=  DN(b,0)   * Xb;
	a(1)           +=  DN(b,1)   * Xb;
	aPartials(0,0) += DDN(b,0,0) * Xb;
	aPartials(0,1) += DDN(b,0,1) * Xb;
	aPartials(1,1) += DDN(b,1,1) * Xb;
      }
      aPartials(1,0)  = aPartials(0,1); // by symmetry
			
      Material_t& material = p->material;
      // compute shell geometry
      ShellGeometry refgeometry( a, aPartials );
			
      // store the reference geometry in shell geometry class
      material.setRefGeometry(refgeometry);

      //material.updateState(true, true, true);

    }    
    compute(true,true,true);
    
    return;
  }
   

} // namespace voom
