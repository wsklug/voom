/*
 * OPSAsphericity.cc
 *
 *  Created on: Aug 7, 2017
 *      Author: amit
 */
#include <string>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <vector>
#include <fstream>

#include <tvmet/Vector.h>
#include <limits>
#include "Node.h"
#include "Model.h"
#include "Lbfgsb.h"
#include "OPSBody.h"

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkSmartPointer.h>
#include "HelperFunctions.h"

using namespace tvmet;
using namespace std;
using namespace voom;

int main(int argc, char* argv[])
{
    clock_t t1, t2;
    t1 = clock();
    if (argc != 2) {
        cout << "usage: " << argv[0] << " <filename>\n";
        return -1;
    }

    string inputFileName = argv[1];
    string fname = inputFileName.substr(0, inputFileName.find("."));
    string rName;
    std::stringstream sstm;

    //For Morse material
    double D_e, re, s, b, gamma;
    double percentStrain = 10;
    bool projectOnSphere = false, volumeConstraintOn = false;
    double initialSearchRad = 1.0, finalSearchRad = 1.2;       

    //Default Values
    s = log(2.0)*(100/14);
    D_e = 1.0;
    re = 1.0;
    gamma = 1.0;

    //Read epsilon and percentStrain from input file. percentStrain is
    //calculated so as to set the inflection point of Morse potential
    //at a fixed distance relative to the equilibrium separation
    //e.g. 1.1*R_eq, 1.5*R_eq etc.
    std::ifstream miscInpFile("miscInp.dat");
    assert(miscInpFile);
    string temp;
    miscInpFile
            >> temp >> D_e
            >> temp >> re
            >> temp >> percentStrain            
            >> temp >> b           
            >> temp >> initialSearchRad
            >> temp >> finalSearchRad
            >> temp >> projectOnSphere
            >> temp >> volumeConstraintOn;
    miscInpFile.close();

    s = (100 / (re*percentStrain))*log(2.0);   
    struct OPSParams props = {D_e, re, s, b, 0.0, 0.0, gamma};

    vtkSmartPointer<vtkPolyDataReader> reader =
            vtkSmartPointer<vtkPolyDataReader>::New();
    vtkSmartPointer<vtkPolyData> mesh;

    reader->SetFileName(inputFileName.c_str());
    reader->Update();
    mesh = reader->GetOutput();

    // create vector of nodes
    int dof = 0;
    std::vector<OPSNode* > nodes;
    std::vector<NodeBase*> baseNodes;
    double Ravg = 0.0;

    // read in points
    for (int i = 0; i < mesh->GetNumberOfPoints(); i++) {
        int id = i;
        NodeBase::DofIndexMap idx(6);
        Vector3D X;
        OPSNode* n;
        for (int j = 0; j < 6; j++) idx[j] = dof++;
        mesh->GetPoint(i, &(X[0]));
        n = new OPSNode(id, idx, X);
        Ravg += tvmet::norm2(X);
        nodes.push_back(n);
        baseNodes.push_back(n);
    }
    assert(nodes.size() != 0);
    int numNodes = nodes.size();
    Ravg /= numNodes;
    cout << "Number of nodes: " << numNodes << endl
         << "Initial radius: " << Ravg << endl;

    OPSBody* bd = new OPSBody( nodes, props, initialSearchRad );

    // Calculate side lengths average of the imaginary equilateral triangles
    double EdgeLength = bd->getAverageEdgeLength();

    // Rescale size of the capsid by the average equilateral edge length
    for (int i = 0; i < numNodes; i++) {
        Vector3D X;
        X = nodes[i]->referencePosition();
        X *= 1.0 / EdgeLength;
        nodes[i]->setReferencePosAndRotVec(X);
        nodes[i]->setDeformedPosAndRotVec(X);
    }
    //Recalculate edge lengths and capsid radius
    bd->updateSearchRadius(finalSearchRad);
    EdgeLength = bd->getAverageEdgeLength();
    bd->updatePolyDataAndKdTree();
    bd->updateNeighbors();

    Ravg = bd->getAverageRadius();

    int numSolverDOFs = 6*numNodes;
    //********************************* ENABLE VOLUME CONSTRAINT ***********************//
    double volConstraint = 0.0;
    if(volumeConstraintOn){
        volConstraint = bd->calcAvgVolume();
        cout << "Prescribed Volume = "<< volConstraint << std::endl;
        NodeBase::DofIndexMap idx(1);
        idx[0] = dof++;
       MultiplierNode *pressNode = new MultiplierNode(numNodes,idx,1.0);
       baseNodes.push_back( pressNode );
       bd->setVolumeConstraint( pressNode, volConstraint);
       numSolverDOFs++;
    }
    //----------------------------------------------------------------------------------//


    std::cout << "Radius of capsid after rescaling = " << Ravg << endl;

    if (projectOnSphere) {
        //Project points to a sphere of radius Ravg
        for (int i = 0; i < numNodes; i++) {
            Vector3D X;
            X = nodes[i]->referencePosition();
            X *= Ravg / (tvmet::norm2(X));
            nodes[i]->setReferencePosAndRotVec(X);
            nodes[i]->setDeformedPosAndRotVec(X);
        }

        //Recalculate edge lengths and capsid radius
        bd->updatePolyDataAndKdTree();
        bd->updateSearchRadius(finalSearchRad);
        bd->updateNeighbors();
        EdgeLength = bd->getAverageEdgeLength();
    }
    //***************** READ ENERGY COEFFICIENTS FROM FILE ********************
    std::ifstream coeffsFile("coeffs.dat");
    assert(coeffsFile);
    std::vector<double> coeffVec;
    double currGamma;
    while (coeffsFile >> currGamma) {
        coeffVec.push_back(currGamma);
    }
    coeffsFile.close();

    std::string dataOutputFile;
    sstm << fname << "-output.dat";
    dataOutputFile = sstm.str();
    sstm.str("");
    sstm.clear();

    ofstream myfile;
    myfile.open(dataOutputFile.c_str());
    myfile << "#Step" << "\t"
           << "Gamma" << "\t"
           << "Radius" << "\t"
           << "Asphericity" << "\t"
           << "Volume" << "\t"
           << "MorseEnergy" << "\t"           
           << "NormalityEn" << "\t"
           << "CircularityEn" << "\t"
           << "PressureEnergy" << "\t"
           << "TotalFunctional"
           << std::endl;

    // Update the Morse parameters
    s = (100 / (EdgeLength*percentStrain))*log(2.0);    
    bd->updateProperty( OPSBody::r, EdgeLength );
    bd->updateProperty( OPSBody::av, s );

    //Create Model
    Model::BodyContainer bdc;
    bdc.push_back(bd);
    Model model(bdc, baseNodes);

    //Create a l-BFGS-b solver
    int m = 7;
    int maxIter = 1e5;
    double factr = 10.0;
    double pgtol = 1e-7;
    int iprint = 2000;
    Lbfgsb solver(numSolverDOFs, m, factr, pgtol, iprint, maxIter, false);


    bool checkConsistency = true;
    if (checkConsistency) {
        std::cout << "Checking consistency......" << std::endl;
        //bd->checkConsistency(true);
        model.checkConsistency(true,false);
    }


    //***************************  SOLUTION LOOP ***************************
    for(int z=0; z < coeffVec.size(); z++){

        cout<<"ITERATION: Z = "<< z <<std::endl;
        cout<<std::endl;

        currGamma = coeffVec[z];
        bd->updateProperty(OPSBody::G, currGamma);

        solver.solve( &model );
        bd->updatePolyDataAndKdTree();
        bd->updateNeighbors();

        myfile << z << "\t"
               << currGamma << "\t"
               << bd->getAverageRadius() << "\t"
               << bd->getAsphericity() << "\t"               
               << bd->calcAvgVolume() << "\t"
               << bd->getMorseEnergy() << "\t"               
               << bd->getNormalityEnergy() << "\t"
               << bd->getCircularityEnergy() << "\t"
               << bd->getPVEnergy() << "\t"
               << solver.function()
               << std::endl;

        sstm << fname << "-step-" << z << ".vtk";
        rName = sstm.str();
        bd->printParaview(rName);
        sstm.str("");
        sstm.clear();

    }

    myfile.close();
    t2 = clock();
    float diff((float)t2 - (float)t1);
    std::cout << "Solution loop execution time: " << diff / CLOCKS_PER_SEC
              << " seconds" << std::endl;

    //Release the dynamically allocated memory
    delete bd;
    for(vector<OPSNode*>::const_iterator i=nodes.begin(); i!=nodes.end(); ++i){
        delete *i;
    }
    nodes.clear();
}
