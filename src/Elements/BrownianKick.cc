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
#include "BrownianKick.h"

namespace voom
{
    
    // Constructor
    BrownianKick::BrownianKick(const NodeContainer &defNodes, 
                               double Cd, double D, double dt )
    : _nodes(defNodes), _Cd(Cd), _D(D), _dt(dt) {
        
        // seed random number generator
        _rng.seed((unsigned int)time(0));
        
        // set the number of nodes
        _nodeCount = _nodes.size();
        _delta_xB.resize(_nodeCount);
        
        //seed the nodal random number generator
        ranlib::DiscreteUniform<int> dis( _nodeCount);
        _dis = &dis;
        _dis->seed((unsigned int)time(0));
        
    }
    
    //Destructor
    BrownianKick::~BrownianKick(){}
    
    void BrownianKick::updateSerialKick(){
        unsigned int randomNode = _dis->random();
        for(int i=0; i < _nodeCount; i++){
            if(i == randomNode){
                Vector3D xi(_rng.random(),
                            _rng.random(),_rng.random());
                _delta_xB[i] = xi*sqrt(_D*_dt);
            }
            else{
                _delta_xB[i] = 0.0,0.0,0.0;
            }
        }
    }
    
    void BrownianKick::updateParallelKick(){
        
        for(int i=0; i < _nodeCount; i++){
            Vector3D xi(_rng.random(),
                        _rng.random(),_rng.random());
            _delta_xB[i] = xi*sqrt(_D*_dt);
        }
    }
    
    void BrownianKick::updateProjectedKick(){
        
        for(int i=0; i < _nodeCount; i++){
            Vector3D xi(_rng.random(),
                        _rng.random(),_rng.random());
        
            Vector3D currPoint = _nodes[i]->point();
            
            Vector3D tempVec( 0, 0, 0);
            
            tempVec = currPoint + xi;
            xi = tempVec;
            
            tempVec = ( xi / tvmet::norm2( xi ) )*tvmet::norm2( currPoint );
            
            xi = tempVec - currPoint;
            
            _delta_xB[i] = xi*sqrt(_D*_dt);
        }
    }
    
    //Rigid rotations as kicks
    void BrownianKick::updateRotationKick(){
        Vector3D tempAxis(_rng.random(), _rng.random(),_rng.random());
        Vector3D axis;
        axis = tvmet::normalize( tempAxis );
        double angle = M_PI_2;        
        double cos_t = cos( angle );
        double sin_t = sin( angle );
        
        for( int i=0; i < _nodeCount; i++ ){
            Vector3D v = _nodes[i]->point();
            Vector3D v_rot, kick;
            v_rot = v*cos_t + tvmet::cross(axis,v)*sin_t +
                tvmet::dot(axis,v)*(1-cos_t)*axis;            
            kick = v_rot - v;            
            _delta_xB[i] = kick;//*sqrt( _D*_dt );
        }
    }
    
	void BrownianKick::update2DKick() {
		for (int i = 0; i < _nodeCount; i++) {
			Vector3D xi(_rng.random(),
				_rng.random(), 0.0);
			_delta_xB[i] = xi*sqrt(_D*_dt);
		}
	}

    // Do mechanics on element; compute energy, forces, and/or stiffness.
    void BrownianKick::compute(bool f0, bool f1, bool f2) {
        
        if( f0 ) {
            _energy = 0;
            for(int i=0; i < _nodeCount; i++){
                for(int j=0; j < 3; j++){
                    _energy += -_D*_delta_xB[i][j]*_nodes[i]->getPoint(j);
                }
            }
        }
        
        if( f1 ) {
            Vector3D f;
            for(int i=0; i < _nodeCount; i++){
                //f = -_Cd * _delta_xB[i];
				f = -_D * _delta_xB[i];
                _nodes[i]->updateForce(f);
            }
        }
        
        return;
        
    }
    
} // namespace voom
