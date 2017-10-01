#ifndef _OPSBODY_H_
#define _OPSBODY_H_

#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <regex>

#include <tvmet/Matrix.h>
#include "Body.h"
#include "voom.h"
#include "Node.h"
#include "VoomMath.h"

#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkKdTree.h>
#include <vtkIdFilter.h>
#include <vtkDelaunay3D.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkLoopSubdivisionFilter.h>
#include <vtkMassProperties.h>

#if VTK_MAJOR_VERSION < 6
#define SetInputData SetInput
#endif

using namespace std;

namespace voom {
/*!
 OPSBody is a concrete class derived from Body, implementing
 a oriented-particle system with no co-circularity or co-normality terms.
 Reference:
         Szeliski, Richard and Tonnesen, D. (1992). Surface Modeling with
         Oriented Particle Systems.	Siggraph ’92, 26(2), 160.
         https://doi.org/10.1017/CBO9781107415324.004
 */

struct OPSParams{
        //Default values
        OPSParams():D_e(1.0), r_e(1.0), a(6.9314718056),//Fracture strain = 10%
                    b(1.0), alpha(1.0), beta(1.0), gamma(1.0){}

        //Copy from another instance
        OPSParams(const OPSParams &p){
                D_e = p.D_e; r_e = p.r_e;
                a = p.a; b = p.b;
                alpha = p.alpha; beta = p.beta; gamma = p.gamma;
        }
        // Initialize from a vector of doubles
        OPSParams( double E, double r, double av, double bv,
                   double p, double q, double s): D_e(E), r_e(r), a(av), b(bv),
            alpha(p), beta(q), gamma(s){}
        /* Morse potential parameters:
         * epsilon = equilibrium energy, r_e = equilibrium separation
         * a = Morse potential width-controlling parameter
         * b = standard deviation of the Gaussian kernel in phi_N and phi_C
         */
        double D_e, r_e, a, b;
        double alpha, beta, gamma;
};

class OPSBody: public Body {
public:

        typedef vector<OPSNode*>::const_iterator opsNodeIterator;
        typedef tvmet::Vector<double,6> Vector6D;
        typedef tvmet::Matrix<double,3,3> Matrix3X3;

        enum Property {E, r, av, bv, A, B, G};

        OPSBody(){}

        //! Construct body from
        OPSBody(const vector<OPSNode*> & nodes, OPSParams &p, double r);

        //! Destructor
        ~OPSBody() {}

        //!Construct vtkPolyData from nodes
        void updatePolyDataAndKdTree();

        //! Do mechanics on Body
        void compute(bool f0, bool f1, bool f2);

        void updateNeighbors();

        double getAverageEdgeLength();

        double getAverageRadius();

        double getVolume(){ return _volume;}

        void calcVolumeAndDerivative();

        double calcVolume();

        double calcAvgVolume(){
            //Calculate (4/3)*M_PI*Ravg^3
            double R = getAverageRadius();
            double avgVol = 4.1887902047863905*R*R*R;
            return avgVol;
        }

        double getAvgNumNeighbors(){
            double num = 0.0;
            for(int i=0; i < _numNodes; i++){
                num += _neighbors[i]->GetNumberOfIds();
            }
            num /= _numNodes;
            return num;
        }

        std::vector<Vector3D> calcVolumeDerivative();

        std::vector<Vector3D> calcAvgVolDerivative();

        double getAsphericity();

        double getLoopAsphericity();

        //! Return the energy of the body
        double energy() const { return _energy;}

        //!Return the current OPS properties
        OPSParams getProperties(){return _prop;}

        //!Update the current OPS properties
        void updateProperties(OPSParams p){_prop = p;}

        //!Return the current search radius
        double getSearchRadius(){return _searchR;}

        //!Update the search radius
        void updateSearchRadius( double r ){ _searchR = r;}

        //!Get the polydata associated with the OPS
        vtkSmartPointer<vtkPolyData> getPolyData(){ return _polyData;}

        //! General printing of a Paraview file
        void printParaview(const string name) const;

        //! Printing of a Paraview file after applying Kabsch algorithm
        void printParaview(const string name, Eigen::Matrix3Xd) const;

        void pushBack(Element* e) { _elements.push_back(e);}

        //! Return the mean squared displacement
        double msd();

        //!Get nearest neighbor for rmsd calculation
        std::vector<int> getInitialNearestNeighbor(){
                return _initialNearestNeighbor;
        }

        void setVolumeConstraint(MultiplierNode *n, double vol){
            _PV = n;
            _volConstraintOn = true;
            _volConstraint = vol;
        }

        void updateVolumeConstraint(double vol){
            _volConstraint = vol;
        }

        void updateProperty(Property p, double val);

        double getMorseEnergy(){return _morseEn;}

        double getNormalityEnergy(){return _normalEn;}

        double getCircularityEnergy(){return _circEn;}

        double getViscoEnergy(){return _viscoEn;}

        double getBrownEnergy(){return _brownEn;}

        double getPVEnergy(){return _PVen;}

        static Vector3D rotVecToNormal( Vector3D u);

        static Matrix3X3 diffNormalRotVec( Vector3D vi);

        //!OPS kernel and potential functions
        double morse(Vector3D xi, Vector3D xj);
        double psi(Vector3D xi, Vector3D xj);        
        double phi_n(Vector3D vi, Vector3D vj);
        double phi_c(Vector3D vi, Vector3D vj, Vector3D xi, Vector3D xj);

        //!OPS kernel and potential derivatives
        Vector3D DmorseDxi(Vector3D xi, Vector3D xj);
        Vector3D DmorseDxj(Vector3D xi, Vector3D xj){
            Vector3D temp(0.0), ans(0.0) ;
            temp = DmorseDxi(xi,xj);
            ans = temp*(-1.0);
            return ans;
        }

        Vector3D DpsiDxi(Vector3D xi, Vector3D xj);
        Vector3D DpsiDxj(Vector3D xi, Vector3D xj){
            Vector3D temp(0.0), ans(0.0) ;
            temp = DpsiDxi(xi,xj);
            ans = temp*(-1.0);
            return ans;
        }

        Vector3D Dphi_nDvi(Vector3D vi, Vector3D vj);
        Vector3D Dphi_nDvj(Vector3D vi, Vector3D vj);

        Vector3D Dphi_cDvi(Vector3D vi, Vector3D vj, Vector3D xi, Vector3D xj);
        Vector3D Dphi_cDvj(Vector3D vi, Vector3D vj, Vector3D xi, Vector3D xj);
        Vector3D Dphi_cDxi(Vector3D vi, Vector3D vj, Vector3D xi, Vector3D xj);
        Vector3D Dphi_cDxj(Vector3D vi, Vector3D vj, Vector3D xi, Vector3D xj);

protected:
        const vector<OPSNode*> _opsNodes; // Nodes
        int _numNodes;
        double _searchR; // Search radius
        double _radius;
        double _volume;
        OPSParams _prop; // Parameters for the OPS
        double _morseEn;        
        double _normalEn;
        double _circEn;
        double _brownEn;
        double _viscoEn;
        bool _volConstraintOn;
        double _volConstraint;
        double _PVen;
        MultiplierNode *_PV;
        std::vector<Vector3D> _volDiff;
        vtkSmartPointer<vtkKdTree> _kdTree;
        vtkSmartPointer<vtkPolyData> _polyData;
        vector<vtkSmartPointer<vtkIdList> > _neighbors;
        std::vector<int> _initialNearestNeighbor; // Needed to find rmsd
};//  OPS body

}
#endif /* _OPSBODY_H_ */

