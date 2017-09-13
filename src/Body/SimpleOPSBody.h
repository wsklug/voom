#ifndef _SIMPLEOPSBODY_H_
#define _SIMPLEOPSBODY_H_

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
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>

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

struct SimpleOPSParams{
        //Default values
        SimpleOPSParams():alphaM(1.0), alphaP(1.0),
                        epsilon(1.0), r_e(1.0), s(6.9314718056),//Fracture strain = 10%
                        K(1.0), a(1.0), b(1.0){}

        //Copy from another instance
        SimpleOPSParams(const SimpleOPSParams &p){
                alphaM = 1.0; alphaP = p.alphaP;
                epsilon = p.epsilon; r_e = p.r_e;
                s = p.s; K = p.K; a = p.a; b = p.b;
        }
        // Initialize from a vector of doubles
        SimpleOPSParams( double aM, double aP, double E, double r,
                        double sv, double Kv, double av, double bv):alphaM(aM),alphaP(aP),
                                        epsilon(E),r_e(r),s(sv),K(Kv),a(av),b(bv){}
        // Weights for the potentials involved in total energy
        double alphaM, alphaP;
        /* Morse potential parameters:
         * epsilon = equilibrium energy, r_e = equilibrium separation
         * s = Morse potential width-controlling parameter
         */
        double epsilon, r_e, s;
        // Kernel parameters
        double K, a, b;
};

class SimpleOPSBody: public Body {
public:

        typedef vector<OPSNode*>::const_iterator opsNodeIterator;
        typedef tvmet::Vector<double,6> Vector6D;

        enum Property {aM, aP, E, r, sv, Kv, av, bv};

        SimpleOPSBody(){}

        //! Construct body from
        SimpleOPSBody(const vector<OPSNode*> & nodes, SimpleOPSParams &p, double r);

        //! Destructor
        ~SimpleOPSBody() {}

        //!Construct vtkPolyData from nodes
        void updatePolyDataAndKdTree();

        //! Do mechanics on Body
        void compute(bool f0, bool f1, bool f2);

        void updateNeighbors();

        double getAverageEdgeLength();

        double getAverageRadius();

        double getAsphericity();

        double getLoopAsphericity();

        //! Return the energy of the body
        double energy() const { return _energy;}

        //!Return the current OPS properties
        SimpleOPSParams getProperties(){return _prop;}

        //!Update the current OPS properties
        void updateProperties(SimpleOPSParams p){_prop = p;}

        //!Return the current search radius
        double getSearchRadius(){return _searchR;}

        //!Update the search radius
        void updateSearchRadius( double r ){ _searchR = r;}

        //!Get the polydata associated with the OPS
        vtkSmartPointer<vtkPolyData> getPolyData(){ return _polyData;}

        //!Get the K-d tree associated with the OPS
        vtkSmartPointer<vtkKdTree> getKdTree(){ return _kdTree;}

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

        void updateProperty(Property p, double val);

        double getMorseEnergy(){return _morseEn;}
        double getPlanarityEnergy(){return _planarEn;}

        //!OPS kernel and potential functions
        double morse(Vector3D xi, Vector3D xj);
        double psi(Vector3D vi, Vector3D xi, Vector3D xj);
        double phi_p(Vector3D vi, Vector3D xi, Vector3D xj);

        //!OPS kernel and potential derivatives
        Vector3D DmorseDxi(Vector3D xi, Vector3D xj);
        Vector3D DmorseDxj(Vector3D xi, Vector3D xj);

        Vector3D DpsiDvi(Vector3D vi, Vector3D xi, Vector3D xj);
        Vector3D DpsiDxi(Vector3D vi, Vector3D xi, Vector3D xj);
        Vector3D DpsiDxj(Vector3D vi, Vector3D xi, Vector3D xj);

        Vector3D Dphi_pDvi(Vector3D vi, Vector3D xi, Vector3D xj);
        Vector3D Dphi_pDxi(Vector3D vi, Vector3D xi, Vector3D xj);
        Vector3D Dphi_pDxj(Vector3D vi, Vector3D xi, Vector3D xj);

protected:
        const vector<OPSNode*> _opsNodes; // Nodes
        int _numNodes;
        int _numBadTri;
        double _searchR; // Search radius
        SimpleOPSParams _prop; // Parameters for the OPS
        double _morseEn;
        double _planarEn;
        vtkSmartPointer<vtkKdTree> _kdTree;
        vtkSmartPointer<vtkPolyData> _polyData;
        vector<vtkSmartPointer<vtkIdList> > _neighbors;
        std::vector<int> _initialNearestNeighbor; // Needed to find rmsd
};// Simple OPS body

}
#endif /* _SIMPLEOPSBODY_H_ */

