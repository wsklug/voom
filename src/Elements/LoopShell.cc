// -*- C++ -*-
//----------------------------------------------------------------------
//
//                 William S. Klug, Feng Feng, Lin Ma
//                University of California Los Angeles
//                 (C) 2004-2008 All Rights Reserved
//
//----------------------------------------------------------------------

////////////////////////////////////////////////////////////////////
// LAPACK subroutine for computing eigenvalues and eigenvectors
//
// define prototype of LAPACK functions
//extern "C" void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda,
//	double *w, double *work, int *lwork, int *info);

//#define _DEBUG_

namespace voom
{

template<class Material_t>
void LoopShell<Material_t>::compute(bool f0, bool f1, bool f2)
{
	//
	// initialize things
	//
	// if applied volume constraint, _volume is needed for computing
	// if applied area constraint, _area is needed for computing
	// both energy and forces
	_volume = 0.0;
	_area = 0.0;
	_totalCurvature = 0.0;
	// _trC = 0.0;

	// **************************************************************
	// Lin: I have removed all of the *spring // stuff from the
	// element.  Really the appropriate way *to // implement the
	// springs is in a new element class.  Then when you want to add
	// your spring energy you can just add spring elements to the body
	// (using the *Body::pushBack(Element *) method. - WSK
	// **************************************************************

	if( f0 ) {
		_energy = 0.0;
		_strainEnergy = 0.0;
		_stretchingEnergy = 0.0;
		_bendingEnergy = 0.0;
		_work = 0.0;
	}

	if( f1 ) {
		Vector3D zero(0);
		_internalForce = zero;
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

		const LoopShellShape::FunctionArray & N
		= p->shape.functions();
		const LoopShellShape::DerivativeArray & DN
		= p->shape.derivatives();
		const LoopShellShape::SecondDerivativeArray & DDN
		= p->shape.secondDerivatives();

		for (int b = 0; b < _nodes.size(); b++){
			const DeformationNode<3>::Point & xb = _nodes[b]->point();
			x 	       +=   N(b)     * xb;
			a(0)           +=  DN(b,0)   * xb;
			a(1)           +=  DN(b,1)   * xb;
			aPartials(0,0) += DDN(b,0,0) * xb;
			aPartials(0,1) += DDN(b,0,1) * xb;
			aPartials(1,1) += DDN(b,1,1) * xb;
		}
		aPartials(1,0)  = aPartials(0,1); // by symmetry


		Material_t& material = p->material;
		// compute shell geometry
		ShellGeometry geometry( a, aPartials );
		const Vector3D& d = geometry.d();
		const tvmet::Vector< Vector3D, 2 >& aDual = geometry.aDual();
		const tvmet::Vector<Vector3D, 2>& dPartials = geometry.dPartials();
		const tvmet::Matrix<Vector3D, 2, 2>& dualPartials = geometry.aDualPartials();

		// store the deformed geometry in shell geometry class
		material.setGeometry(geometry);

		// compute strain energy, stress and moment resultants
		material.updateState(f0, f1, f2);

		const double metric = geometry.metric();
		// const double refMetric = refgeometry.metric();
		const double refMetric = ( material.refShellGeometry()).metric();
		const double jacobian = metric/refMetric;
		const double weight =  metric * p->weight;

		// compute area for the area constraint energy
		_area += weight;

		// compute volume for the volume constraint energy
		_volume +=  dot(d,x) * weight / 3.0;

		// compute total curvature for constraint energy
		_totalCurvature += material.meanCurvature()*weight;

		// compute energy
		if ( f0 ){
			// compute strain energy
			_strainEnergy += material.energyDensity() * weight;
			_stretchingEnergy += material.stretchingEnergy() * weight;
			_bendingEnergy += material.bendingEnergy() * weight;
		}

		// compute forces
		if ( f1 ) {
			//
			// get stress and moment resultants
			const tvmet::Vector< Vector3D, 3 >& sr = material.stressResultants();
			const tvmet::Vector< Vector3D, 2 >& mr = material.momentResultants();
			// ***********************************************************
			// Lin: this is a poor design.  The materials shouldn't need
			// to compute total curvature stuff (otherwise ALL shell
			// materials would HAVE to do that calculation, and that would
			// be silly).  You can do this right here in element
			// code. -WSK
			// ***********************************************************

			const tvmet::Vector< Vector3D, 3 >& stcr = material.totalCurvatureStressResultants();

			double pressure = _pressureNode->point();
			double tension = _tensionNode->point();
			double totalCurvatureForce = _totalCurvatureNode->point();

			// loop for all nodes to compute forces
			for (int a=0; a<_nodes.size(); a++) {

				// compute internal forces

				// calculate the gradient of the derivatives of the director
				// w.r.t curvilinear coords
				for ( int alpha = 0; alpha < 2; alpha++){
					// Stress Resultant part
					Vector3D ftmp;
					ftmp = sr(alpha) *  DN(a,alpha) * weight;
					_internalForce(a)  += ftmp;
					// Moment Resultant part
					for(int beta=0; beta<2; beta++) {
						ftmp = -  dot(mr(alpha),aDual[beta])*
								( DDN(a,alpha,beta)*d + DN(a,beta)*dPartials(alpha) )
								-  dot(mr(alpha),dualPartials(beta,alpha))*DN(a,beta)*d;
						ftmp *= weight;
						_internalForce(a) += ftmp;
					}
				}

				// *********************************************************
				// Lin: you should only do these calculations if there is a
				// total curvature constraint.  Check before wasting
				// time. -WSK
				// *********************************************************

				for ( int alpha = 0; alpha < 2; alpha++){
					Vector3D ftmp;
					ftmp = stcr(alpha) *  DN(a,alpha) * weight * totalCurvatureForce;
					_internalForce(a)  += ftmp;

					for(int beta=0; beta<2; beta++) {
						ftmp = 1.0/2.0*  dot(aDual[alpha],aDual[beta])*
								( DDN(a,alpha,beta)*d + DN(a,beta)*dPartials(alpha) )*totalCurvatureForce
								+1.0/2.0*  dot(aDual[alpha],dualPartials(beta,alpha))* DN(a,beta)*d *totalCurvatureForce;
						ftmp *= weight;
						_internalForce(a) += ftmp;
					}
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


			} // end nodes loop

		} // end force calcs

		// compute stiffness matrix
		if( f2 ) {
			// 	std::cerr << std::endl << std::endl << "\t"
			// 		  << "Aaaarrrrrrggggggh!  No stiffness matrix yet in voom::LoopShell!"
			// 		  << std::endl << std::endl;
		}

	} // end quadrature loop

	if(f0) {
		_energy = _strainEnergy;
	}



	if(f1) {
		int a=0, ia=0;
		for(NodeIterator na=_nodes.begin();  na!=_nodes.end(); na++, a++)
			for(int i=0; i<3; i++, ia++) {
				double f_ia = _internalForce(a)(i)
	    	    		  + _pressureForce(a)(i)
						  + _tensionForce(a)(i);
				(*na)->addForce( i, f_ia );
			}
	}


}


template<class Material_t>
void LoopShell<Material_t>::updateRefConfiguration() {
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

		bool flat=false;
		if( flat ) {
			// this is a hack to treat the reference geometry as if the
			// element were its current shape but embedded locally within
			// a flat sheet.
			//
			//     ^ s^2
			//    /
			//   2
			//  / \
			// 0 - 1 --> s^1
			//
			// X = (1-s^1-s^2)*X_0 + s^1*X_1 + s^2*X_2
			// G_1 = X_{,1} = X_1-X_0
			// G_2 = X_{,2} = X_2-X_0
			//
			a(0) = _nodes[1]->position() - _nodes[0]->position();
			a(1) = _nodes[2]->position() - _nodes[0]->position();

		} else {
			// compute reference geometry correctly by interpolation from
			// the reference coordinates of the nodes
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

		}

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

//added to set the reference configuration explicitly
template<class Material_t>
void LoopShell<Material_t>::SetRefConfiguration(double edgelen) {
	for(QuadPointIterator p=_quadPoints.begin();
			p!=_quadPoints.end(); p++){
		tvmet::Vector< Vector3D, 2 > a;
		tvmet::Matrix< Vector3D, 2, 2 > aPartials;
		a(0) = edgelen, 0.0, 0.0;
		// WSK: an example of what not to do... Seriously, when you can
		// evaluate a complicated expression in closed form, you should
		// do so.
		//  a(1) = edgelen*cos(M_PI/3.),edgelen*sin(M_PI/3.),0.;
		a(1) = edgelen*0.5, edgelen*0.5*sqrt(3.0), 0.0;

		Material_t& material = p->material;

		// compute shell geometry
		ShellGeometry refgeometry( a, aPartials );

		// store the reference geometry in shell geometry class
		material.setRefGeometry(refgeometry);
	}
	compute(true,true,true);
	return;
}

template<class Material_t>
Vector3D LoopShell<Material_t>::computePosition(const double s1, const double s2)
{
	//CornerValences v(6,6,6);
	tvmet::Vector<unsigned int, 3> v(6,6,6);
	Array1D paraCoords(2);
	paraCoords = s1, s2;
	const int nodes = v(0) + v(1) + v(2) - 6;
	//
	// create a loop shell shape function object
	LoopShellShape ls( nodes, v, paraCoords );
	//
	// new position
	Vector3D pos(0.0);
	NodeIterator p = _nodes.begin();
	for(; p != _nodes.end(); p ++){
		const int i = std::distance(_nodes.begin(), p);
		pos += ls.functions()(i) * (*p)->point();
	}

	return pos;
}


template < class Material_t >
void LoopShell<Material_t>::checkPositions()
{
	blitz::Array< tvmet::Vector<double, 2>, 1 > para(7);
	para(0) = 1.0/3.0;
	para(1)(0) = 1.0; para(1)(1) = 0.0;
	para(2)(0) = 0.0; para(2)(1) = 0.0;
	para(3)(0) = 0.0; para(3)(1) = 1.0;
	para(4)(0) = 0.0; para(4)(1) = 0.5;
	para(5)(0) = 0.5; para(5)(1) = 0.5;
	para(6)(0) = 0.5; para(6)(1) = 0.0;

	Vector3D pos;
	for (int i = 0; i < 7; i ++){
		pos = computePosition(para(i)(0), para(i)(1));
		std::cout << "new position at ( " << para(i)(0)
				<< ", " << para(i)(1)
				<< " ) are "
				<< "( " << pos(0)
				<< ", " << pos(1)
				<< ", " << pos(2)
				<< ") " << std::endl;
	}

	std::cout << "input parametric coords (two doubles):" << std::endl;
	double s1, s2;
	std::cin >> s1 >> s2;

	pos = computePosition(s1, s2);

	std::cout << "new position at ( " << s1
			<< ", " << s2
			<< " ) are "
			<< "( " << pos(0)
			<< ", " << pos(1)
			<< ", " << pos(2)
			<< ") " << std::endl;

}

template<class Material_t>
double LoopShell<Material_t>::meancurvature()
{
	// Assumption: Only one quadrature point
	if(_quadPoints.size()!=1) std::cout<<"Expected only one quadrature point for returning mean curvature"<<std::endl;
	for(QuadPointIterator p=_quadPoints.begin();
			p!=_quadPoints.end(); p++){
		// compute position, basis vector, and derivatives of basis vector
		//Vector3D x(0.0);
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
			const DeformationNode<3>::Point & xb = _nodes[b]->point();
			//x 	       +=   N(b)     * xb;
			a(0)           +=  DN(b,0)   * xb;
			a(1)           +=  DN(b,1)   * xb;
			aPartials(0,0) += DDN(b,0,0) * xb;
			aPartials(0,1) += DDN(b,0,1) * xb;
			aPartials(1,1) += DDN(b,1,1) * xb;
		}
		aPartials(1,0)  = aPartials(0,1); // by symmetry
		ShellGeometry defgeometry( a, aPartials );

		const tvmet::Vector< Vector3D, 2 >& dual = defgeometry.aDual();
		const tvmet::Vector<Vector3D, 2>& dPartials = defgeometry.dPartials();
		//mean curvature
		double H = - 0.5 * ( dot(dual(0), dPartials(0)) + dot(dual(1), dPartials(1)) );
		return H;
	}
}

template<class Material_t>
Vector3D LoopShell<Material_t>::PushForwardOperator(Vector3D & Nbar)
{
	int b = 0, i = 0, j = 0, alpha = 0;
	Vector3D zero(0), nbar(0);

	// Loop over quadrature point
	for(QuadPointIterator p = _quadPoints.begin(); p != _quadPoints.end(); p++)
	{
		// compute position, basis vector, and derivatives of basis vector
		tvmet::Vector< Vector3D, 2 > a;
		a = zero, zero;

		const LoopShellShape::DerivativeArray & DN
		= p->shape.derivatives();

		for (b = 0; b < _nodes.size(); b++){
			const DeformationNode<3>::Point & xb = _nodes[b]->point();
			a(0)           +=  DN(b,0)   * xb;
			a(1)           +=  DN(b,1)   * xb;
		}

		const tvmet::Vector<Vector3D, 2> & refDual = (p->material.refShellGeometry()).aDual();

		Tensor3D F(0.0);
		for(alpha = 0; alpha<2; alpha++){
			for(i=0; i<3; i++) {
				for(j=0; j<3; j++) {
					F(i,j) += a(alpha)(i)*refDual(alpha)(j);
				}
			}
		}

		// nbar += F*Nbar*(p->weight);
		nbar = F*Nbar;
	}

	nbar /= norm2(nbar);

	return nbar;
}

} // namespace voom
