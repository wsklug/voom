#if defined(_MSC_VER) || defined(__INTEL_COMPILER)
#define _USE_MATH_DEFINES
#include <cmath>
#endif
#include <string>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <unistd.h>
#include <tvmet/Vector.h>
#include <limits>
#include <algorithm>

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkIdList.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkSphereSource.h>

#if VTK_MAJOR_VERSION < 6
# define SetInputData SetInput
#endif

using namespace std;

void densityPlotter( std::vector<std::string> fileNames,
                     std::string outputFile,
                     int lat_res, int long_res );

int main(int argc, char* argv[]){
    
    if( argc < 7 ) {
        cout << "Usage: densityPlot <base-name> <num of files> <interval>" 
        <<" <outputFile> <lat_res> <long_res>" << endl;
        return(0);
    }
    
    clock_t t1,t2;
    t1=clock();
    
    string baseFileName = argv[1];
    int lastFileNum = std::atoi(argv[2]);
    int interval = std::atoi(argv[3]);
    string outputFile = argv[4];
    int lat_res = std::atoi(argv[5]);
    int long_res = std::atoi(argv[6]);
    
    std::stringstream sstm;
    
    std::vector<std::string> allVTKFiles;
    allVTKFiles.reserve(lastFileNum/interval);
    
    for(int fileNum=0 ; fileNum < lastFileNum; fileNum++){    
        sstm << baseFileName <<"-"<< (fileNum*interval) <<".vtk";
        std::string tempString = sstm.str();
        allVTKFiles.push_back(tempString);
        sstm.str("");
        sstm.clear();
    }
    
    densityPlotter(allVTKFiles, outputFile , lat_res, long_res);
    t2=clock();
    double diff = ((float)t2-(float)t1);
    std::cout<<"Post-processing execution time: "<<diff/CLOCKS_PER_SEC
    <<" seconds"<<std::endl;
    
    return 0;
}

void densityPlotter( std::vector<std::string> fileNames,
                     std::string outputFile,
                     int lat_res, int long_res ){
    
    bool debug = false;
    
    /*
     * Create a sphere that is meshed along latitude and longitude
     */
    vtkSmartPointer<vtkSphereSource> sp = 
    vtkSmartPointer<vtkSphereSource>::New();
    sp->SetRadius( 1.0 );
    sp->SetThetaResolution( lat_res );
    sp->SetPhiResolution( long_res );
    sp->LatLongTessellationOn();
    
    vtkSmartPointer<vtkPolyData> pd = 
    sp->GetOutput();
    sp->Update();
    
    
    /* 
     * In the next few lines we will identify spherical co-ordinate
     * limits for each cell in the sphere
     */  
    std::vector<vector<double> > cellLimits;
    
    vtkSmartPointer<vtkCellArray> bins = pd->GetPolys();
    
    vtkSmartPointer<vtkIdList> points 
    = vtkSmartPointer<vtkIdList>::New();
    
    bins->InitTraversal();
    int cellId = 0;
    
    while( bins->GetNextCell( points ) ){
        
        if ( debug ){
            std::cout<<"Cell Id: " << cellId++
            << std::endl;
        }        
        
        double theta_max = 0, theta_min = 180,
        phi_max = 0, phi_min = 360;
        
        for( int i=0; i < points->GetNumberOfIds(); i++ ){
            
            if ( debug ){
                std::cout<< "\tPoint Id : " << points->GetId(i)
                << std::endl << "\t\t";
            }
            
            //Get Cartesian coordinates for each point
            double *xyz = pd->GetPoint( points->GetId( i ) );
            
            if( debug ){
                std::cout<< xyz[0] << "," << xyz[1] << ","
                << xyz[2] << std::endl << "\t\t";
            }
            
            //Skip the "poles" of the sphere
            if( (std::abs(xyz[0]) < 1e-8) && 
                (std::abs( xyz[1]) < 1e-8 ) ){
                if( std::abs( xyz[2] - 1) < 1e-8 )
                    theta_min = 0;
                if( std::abs( xyz[2] + 1) < 1e-8 )
                    theta_max = 180;
                if ( debug ){
                    std::cout<< std::endl;
                }
                continue;
                }
                
                //Convert to spherical coordinates (phi,theta)
                double phi = ( 180/M_PI )*atan2( xyz[1], xyz[0] );
                double theta = ( 180/M_PI )*acos( xyz[2] );
                
                phi = ( phi < 0 )? (360 + phi) : phi; 
                
                if ( debug ){
                    std::cout<< "Phi = "<< phi << " Theta = "<< theta 
                    << std::endl;
                }
                
                //Compare to update max and min values
                theta_max = std::max( theta, theta_max );
                theta_min = std::min( theta, theta_min );
                phi_max = std::max( phi, phi_max );
                phi_min = std::min( phi, phi_min );
                
        }
        vector<double> temp;
        temp.push_back( phi_min );
        temp.push_back( phi_max );
        temp.push_back( theta_min );
        temp.push_back( theta_max );
        
        if ( debug ){
            std::cout<< "\tPhi_min_max: "<< phi_min <<","<< phi_max
            << std::endl << "\t"
            << "Theta_min_max: "<< theta_min <<","<< theta_max
            << std::endl;
        }
        
        cellLimits.push_back( temp );
        
    }
    
    if ( debug ){
        
        std::cout<< "Printing the bins : "<< std::endl;
        std::cout<< "\tBinId\tPhi_min\tPhi_max\tTheta_min\tTheta_max"<<std::endl;
        
        for (int binIter = 0; binIter < cellLimits.size(); binIter++){
            std::cout<< "\t" << binIter 
            << "\t" << cellLimits[ binIter ][0]
            << "\t" << cellLimits[ binIter ][1] 
            << "\t" << cellLimits[ binIter ][2] 
            << "\t" << cellLimits[ binIter ][3]
            << std::endl;
        }
    }
    /* Read a VTK file, project points to a unit sphere, identify which
     *        cell of the unit sphere created above does each point belong to
     */
    
    vtkSmartPointer<vtkDoubleArray> binDensity =
    vtkSmartPointer<vtkDoubleArray>::New();
    binDensity->SetNumberOfComponents( 1 );
    binDensity->SetNumberOfTuples( pd->GetNumberOfCells() );
    binDensity->SetName( "Density" );
    binDensity->FillComponent( 0, 0.0 );
    
    #ifdef _OPENMP
    # pragma omp parallel for
    #endif
    for(int i=0; i < fileNames.size() ; i++){
        
        std::string fileName = fileNames[ i ];
        assert( ifstream( fileName.c_str() ) );
        
        vtkSmartPointer<vtkPolyDataReader> reader = 
        vtkSmartPointer<vtkPolyDataReader>::New();
        reader->SetFileName( fileName.c_str() );
        
        reader->ReadAllScalarsOn();
        reader->ReadAllVectorsOn();
        
        vtkSmartPointer<vtkPolyData> mesh = reader->GetOutput();  
        reader->Update();
        
        vtkSmartPointer<vtkDataArray> displacements 
        = mesh->GetPointData()->GetVectors("displacements");
        
        bool displacementsExist = false;
        
        if ( displacements.GetPointer() ){
            displacementsExist = true;
        }

        
        if ( debug ){
            std::cout<< "Reading from file : "<< fileName << std::endl;
        }
        
        for( int z=0; z <  mesh->GetNumberOfPoints(); z++ ){
            
            if ( debug ){
                std::cout<< "\tPoint Id = " << z << std::endl;
            }
            
            //Get Cartesian coordinates for each point
            double * xyz = mesh->GetPoint( z );
            tvmet::Vector<double, 3> pos( xyz[ 0 ], xyz[ 1 ], xyz[ 2 ]);
            
            //Get displacements
            double disp[3] = {0};
            
            if ( displacementsExist ){
                displacements->GetTuple( z, disp );
                
                for (int k=0 ; k < 3; k++){
                    xyz[ k ] += disp[ k ];
                }
            }
            
            if ( debug ){
                std::cout<< "\t\tOriginal : " <<xyz[0] << "," << xyz[1] << ","
                << xyz[2] << std::endl;
            }
            
            //Project x,y,z to surface of a unit sphere
            for(int p=0; p < 3; p++){
                xyz[ p ] /= tvmet::norm2( pos ); 
            }
            
            if ( debug ){
                std::cout<< "\t\tNormalized : " <<xyz[0] << "," << xyz[1] << ","
                << xyz[2] << std::endl;
            }
            
            //Convert to spherical coordinates (phi,theta)
            double phi = ( 180/M_PI )*atan2( xyz[1], xyz[0] );
            double theta = ( 180/M_PI )*acos( xyz[2] );
            
            phi = ( phi < 0 )? (360 + phi) : phi; 
            
            if ( debug ){
                std::cout<< "\t\tPhi = "<< phi << " Theta = "<< theta 
                << std::endl;
            }
            
            for(int binId=0; binId < cellLimits.size(); binId++  ){
                double p_min, p_max, t_min, t_max;
                p_min = cellLimits[ binId ][ 0 ];
                p_max = cellLimits[ binId ][ 1 ];
                t_min = cellLimits[ binId ][ 2 ];
                t_max = cellLimits[ binId ][ 3 ];
                
                if ( (p_min <= phi && phi < p_max) &&
                    (t_min <= theta && theta < t_max) ){
                    
                    if ( debug ){
                        std::cout<< "\t\tBin found : "<< binId <<std::endl;
                    }
                    
                    double tempCount = binDensity->GetTuple1( binId );
                //Size of fileNames vector corresponds to number of time steps
                binDensity->SetTuple1( binId, tempCount + (1.0/fileNames.size() ));
                break;
                
                    }
            }
            
        }
        
    }
    
    pd->GetCellData()->AddArray( binDensity );
    
    vtkSmartPointer<vtkPolyDataWriter> wr =
    vtkSmartPointer<vtkPolyDataWriter>::New();
    wr->SetFileName( outputFile.c_str() );
    //wr->SetFileTypeToBinary();
    wr->SetInputData( pd );
    wr->Update();
    wr->Write();
    
                     }
                     
                     
                     
                     
