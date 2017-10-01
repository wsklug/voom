// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                 (C) 2004-2005 All Rights Reserved
//
//----------------------------------------------------------------------
//
// $Log$
// Created 2015/11/25 03:06:00  amit
// 
//
//
//----------------------------------------------------------------------


#include "VoomMath.h"
#include "OPSBrownianKick.h"

namespace voom
{
	//! Constructor 1
    OPSBrownianKick::OPSBrownianKick(const NodeContainer &opsNodes,
                                     double v):_nodes(opsNodes) {

		// seed random number generator
        _rng.seed((unsigned int)time(0));

        _coefficient = 1.41421356237*v; //sqrt(2)

		// set the number of nodes
		_nodeCount = _nodes.size();
        _randomKick.resize(_nodeCount);

		//seed the nodal random number generator
		ranlib::DiscreteUniform<int> dis(_nodeCount);
		_dis = &dis;
		_dis->seed((unsigned int)time(0));

		_cutOff = 1.0;
	}

    //! Constructor 2
    OPSBrownianKick::OPSBrownianKick(const NodeContainer &opsNodes, double v,
                                     double cut): _nodes(opsNodes), _cutOff(cut) {
        
        _coefficient = 1.41421356237*v; //sqrt(2)

        // seed random number generator
        _rng.seed((unsigned int)time(0));
        
        // set the number of nodes
        _nodeCount = _nodes.size();
        _randomKick.resize(_nodeCount);
        
        //seed the nodal random number generator
        ranlib::DiscreteUniform<int> dis( _nodeCount);
        _dis = &dis;
        _dis->seed((unsigned int)time(0));
        
    }
    
    //Destructor
    OPSBrownianKick::~OPSBrownianKick(){}
    
    void OPSBrownianKick::updateSerialKick(){
        unsigned int randomNode = _dis->random();
        for(int i=0; i < _nodeCount; i++){
            if(i == randomNode){
                Vector3D xi(_rng.random(),
                            _rng.random(),_rng.random());
                _randomKick[i] = xi;
            }
            else{
                _randomKick[i] = 0.0,0.0,0.0;
            }
        }
    }

    std::vector<double> OPSBrownianKick::getKickStats()
	{
		double avg = 0.0;//Mean and std-dev
		double max = 0.0, min = 1e10;
		for (int i = 0; i < _nodeCount; i++) {
            double currNorm = tvmet::norm2(_randomKick[i]);
			avg += currNorm;
			max = (currNorm > max)? currNorm : max;
			min = (currNorm < min) ? currNorm : min;
		}
		avg /= _nodeCount;
		std::vector<double> out = {max, min, avg};
		return out;
	}

    /*
     * updateParallelKick() - kicks all particles simultaneously
     */

    void OPSBrownianKick::updateParallelKick(){
        for(int i=0; i < _nodeCount; i++){
            Vector3D xi(_rng.random(),
                        _rng.random(),_rng.random());
            _randomKick[i] = xi;
        }
    }
    
    void OPSBrownianKick::updateProjectedKick(){
        
        for(int i=0; i < _nodeCount; i++){
            Vector3D xi(_rng.random(),
                        _rng.random(),_rng.random());
        
            Vector3D currPoint(0.0),tempVec(0.0);
            currPoint = _nodes[i]->deformedPosition();

            tempVec = currPoint + xi;
            xi = tempVec;
            tempVec = ( xi / tvmet::norm2( xi ) )*tvmet::norm2( currPoint );
            xi = tempVec - currPoint;
            
            _randomKick[i] = xi;
        }
    }
    
    void OPSBrownianKick::truncatedProjectedKick() {

		for (int i = 0; i < _nodeCount; i++) {
            Vector3D xi(_rng.random(),
                _rng.random(), _rng.random());

			Vector3D currPoint(0.0),tempVec(0.0);
			for(int z=0; z < 3; z++) currPoint[z] = _nodes[i]->getPoint(z);

            tempVec = currPoint + xi;
			xi = tempVec;

			tempVec = (xi / tvmet::norm2(xi))*tvmet::norm2(currPoint);

			xi = tempVec - currPoint;

			if (tvmet::norm2(xi) > _cutOff) {
				Vector3D temp(0,0,0);
				temp = xi / tvmet::norm2(xi)*_cutOff;
				xi = temp;
			}
            _randomKick[i] = xi;
		}
	}

    //Rigid rotations as kicks
    void OPSBrownianKick::updateRotationKick(){
        //Vector3D tempAxis(_rng.random(), _rng.random(),_rng.random());
    	Vector3D tempAxis(0.0, 0.0,1.0);
        Vector3D axis;
        axis = tvmet::normalize( tempAxis );
        double angle = M_PI_2;
        double cos_t = cos( angle );
        double sin_t = sin( angle );
        
        for( int i=0; i < _nodeCount; i++ ){
        	Vector3D v(0.0);
        	for(int z=0; z < 3; z++) v[z] = _nodes[i]->getPoint(z);
            Vector3D v_rot, kick;
            v_rot = v*cos_t + tvmet::cross(axis,v)*sin_t +
                tvmet::dot(axis,v)*(1-cos_t)*axis;            
            kick = v_rot - v;            
            _randomKick[i] = kick;//*sqrt( _D*_dt );
        }
    }
    
    void OPSBrownianKick::update2DKick() {
		for (int i = 0; i < _nodeCount; i++) {
            Vector3D xi(_rng.random(),
                _rng.random(), 0.0);
            _randomKick[i] = xi;
		}
	}


    void OPSBrownianKick::update1DKick() {
		for (int i = 0; i < _nodeCount; i++) {
            Vector3D xi(_rng.random(), 0.0, 0.0);
            _randomKick[i] = xi;
		}
	}

    // Do mechanics on element; compute energy, forces, and/or stiffness.
    void OPSBrownianKick::compute(bool f0, bool f1, bool f2) {
        
        if( f0 ) {
            _energy = 0;
            for(int i=0; i < _nodeCount; i++){
                Vector3D x;
                x = _nodes[i]->deformedPosition();
                    _energy += -_coefficient*tvmet::dot(_randomKick[i],x);
            }
        }
        if( f1 ) {
            Vector3D f;
            for(int i=0; i < _nodeCount; i++){
                f = -_coefficient * _randomKick[i];
				for(int j=0; j<3; j++){
					_nodes[i]->addForce(j,f[j]);
				}
            }
        }
        
        return;
        
    }
    
} // namespace voom
