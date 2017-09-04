#include "HelperFunctions.h"
#include <string>
#include <iostream>
#include <vector>
#include <tvmet/Vector.h>
#include <vtkCell.h>
#include <vtkDataSetReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <omp.h>


using namespace tvmet;
using namespace std;
using namespace voom;

typedef tvmet::Vector<double,3> Vector3D;

int main(int argc, char* argv[])
{
    if( argc != 5 ) {
        cout << "argc = " << argc << endl
             << "Usage: dual2 baseName numStart numEnd numOPENMPThreads" << endl;
        return(0);
    }

    ////////////////////////////////////////////////////////////////////
    // Input section
    ////////////////////////////////////////////////////////////////////

    // read in vtk file
    string baseFileName = argv[1];
    int numStart = atoi(argv[2]);
    int numEnd = atoi(argv[3]);
    int numThreads = atoi(argv[4]);

    //Set the number of threads
    omp_set_num_threads(numThreads);

#pragma omp parallel for
    for(int bigI=numStart; bigI <= numEnd; bigI++){
        stringstream sstm;
        string inputFileName, outFileName;
        sstm << baseFileName << "-" << bigI <<".vtk";
        inputFileName = sstm.str();
        sstm.str("");
        sstm.clear();

        vtkSmartPointer<vtkPolyDataReader> reader =
                vtkSmartPointer<vtkPolyDataReader>::New();
        reader->SetFileName( inputFileName.c_str() );
        reader->Update();
        vtkSmartPointer<vtkPolyData> inputMesh = reader->GetOutput();
        int npts = inputMesh->GetNumberOfPoints();

        // get vertex positions
        std::vector< Vector3D > points( npts );
        for(int a=0; a < npts; a++) {
            inputMesh->GetPoint(a, &(points[a](0)));
        }

        //Can we add displacements to vertex positions before generatinng
        //the dual mesh - Amit

        //get displacements, if they exist
        //std::vector< tvmet::Vector<double,3> > displacements( npts );
        string vectorName="displacements";
        vtkSmartPointer<vtkDoubleArray> displacements =
                vtkDoubleArray::SafeDownCast(inputMesh->GetPointData()->
                                             GetVectors(vectorName.c_str()));
        if( displacements.GetPointer() != NULL){
            double currDisp[3];
            for(int a=0; a < npts; a++){
                displacements->GetTuple(a,currDisp);
                points[a][0] += currDisp[0];
                points[a][1] += currDisp[1];
                points[a][2] += currDisp[2];
            }
        }

        vtkSmartPointer<vtkKdTree> kdt = vtkSmartPointer<vtkKdTree>::New();
        vtkSmartPointer<vtkPoints> newPts =
                vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkIntArray> valence =
                vtkSmartPointer<vtkIntArray>::New();
        valence->SetName("Valence");
        valence->SetNumberOfComponents(1);

        vtkSmartPointer<vtkCellArray> cells = inputMesh->GetPolys();
        vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();
        double avgSearchRad = 0.0;

        while( cells->GetNextCell( cellPointIds ) ){
            int numCellPoints = cellPointIds->GetNumberOfIds();
            Vector3D centroid(0.0);
            for(int i=0; i < numCellPoints; i++){
                int currCellPoint = cellPointIds->GetId(i);
                centroid += points[currCellPoint];
            }
            centroid /= numCellPoints;
            newPts->InsertNextPoint( &centroid[0] );
            avgSearchRad += norm2( centroid - points[ cellPointIds->GetId(0) ] );
        }
        avgSearchRad /= cells->GetNumberOfCells();
        kdt->BuildLocatorFromPoints( newPts );

        //factor of safety in search radius;
        avgSearchRad *= 1.5;
        //std::cout<<"Average search radius = " << avgSearchRad << std::endl;

        //Connect the centroids into polygons
        vtkSmartPointer<vtkCellArray> newPolys =
                vtkSmartPointer<vtkCellArray>::New();

        for(int i=0; i < points.size(); i++){

            Vector3D center;
            center = points[i];

            vtkSmartPointer<vtkIdList> currPolyPtIds =
                    vtkSmartPointer<vtkIdList>::New();
            kdt->FindPointsWithinRadius( avgSearchRad, &center[0], currPolyPtIds );

            std::list< neighbors > currPoly;
            // For the 0th centroid
            Vector3D centroid(0.0), vec0(0.0), vecj(0.0), currCross(0.0), axis(0.0);
            double vec0_norm, vecj_norm, sign, currSin, currCos, currAngle;
            int numCellPoints = currPolyPtIds->GetNumberOfIds();
            int currId = currPolyPtIds->GetId(0);
            newPts->GetPoint( currId, &centroid[0] );
            vec0 = centroid - center;
            vec0_norm = norm2( vec0 );
            vec0 /= vec0_norm;
            neighbors pt0(currId, 0.0);
            currPoly.push_back( pt0 );

            // For remaining centroids
            for(int j=1; j < numCellPoints; j++){
                currId = currPolyPtIds->GetId(j);
                newPts->GetPoint( currId, &centroid[0] );
                vecj = centroid - center;
                vecj_norm = norm2( vecj );
                vecj /= vecj_norm;
                currCross = tvmet::cross(vec0, vecj);
                currSin = tvmet::norm2(currCross);
                axis = currCross / currSin;
                sign = tvmet::dot(axis, center);
                currSin = (sign > 0.0) ? currSin : -1.0 * currSin;
                currCos = tvmet::dot(vec0, vecj);
                currAngle = (180 / M_PI) * atan2(currSin, currCos);
                currAngle = (currAngle < 0) ? (360 + currAngle) : currAngle;
                neighbors ptj(currId, currAngle);
                currPoly.push_back(ptj);
            }

            //Sort the list of neigbors and make a polygon
            currPoly.sort();
            newPolys->InsertNextCell( numCellPoints );
            for(std::list<neighbors>::iterator t = currPoly.begin();
                t != currPoly.end(); ++t){
                neighbors n = *t;
                newPolys->InsertCellPoint( n._id );
            }
            valence->InsertNextTuple1( numCellPoints );
        }
        //Assign points and polygons to a new polydata and write it out
        vtkSmartPointer<vtkPolyData> newPolyData =
                vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkPolyDataWriter> writer =
                vtkSmartPointer<vtkPolyDataWriter>::New();

        newPolyData->SetPoints( newPts );
        newPolyData->SetPolys( newPolys );
        newPolyData->GetCellData()->AddArray( valence );
        writer->SetInputData( newPolyData );
        sstm << baseFileName << "-dual-" << bigI <<".vtk";
        outFileName = sstm.str();
        sstm.str("");
        sstm.clear();
        writer->SetFileName(outFileName.c_str());
        writer->Write();
    }
    return 0;
}

