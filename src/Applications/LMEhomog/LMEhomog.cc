// This file reads in the configurations from MD
// either in .dat or .pdf formats
// Then calculates the deformation gradient at each node
// Then calculates the three invariants at each node
// and writes them to the file invariants1.dat *2.dat and *3.dat
// in the format : at each configuration a row of invariants at all the nodes
#include <iostream>
#include <vector>
#include <fstream>
#include <getopt.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace std;
#include "voom.h"
#include "LMEbody.h"
#include "LMEshape.h"
#include "NodeBase.h"
#include "Node.h"
#include "Model.h"
#include "CompNeoHookean.h"
#include<tvmet/Vector.h>
#include<tvmet/Matrix.h>

using namespace voom;

void readconfig_pdb(char* filename, float* x, float* y, float* z, int npts);
void readconfig_dat(char* filename, float* x, float* y, float* z, int npts);
void readconfig_dat(char* filename, float* x, float* y, float* z, int npts, int average);

int main(int argc, char* argv[])
{
  if( argc < 2 ) {
    cout << "Usage: homog inputFileName [-n start_config -N end_config -s searchR -b beta -a average_number]." << endl;
    return(0);
  }
  
  bool verbose=true;
  
  for(int i=0; i<argc; i++) {
    std::cout << std::setw(8) << i << "\t"
	      << argv[i] << std::endl;
  }
  
  //Input file name and other parameters from run command
  string objFileName = argv[1];
  ifstream ifs;
  string modelName = argv[1];
  string inputFileName = modelName;
  bool datfile = true;
  bool pdbfile = false;
  bool badCommandLine = false;
  double searchR = 1.0;
  double beta = 0.8;
  int nItMax = 100;
  double tol = 1.0e-10;
  int StartConfig = 0, EndConfig = 0, avg = 1;
 
  
  for(char option; (option=getopt(argc,argv,"n:N:s:b:a:t:i:")) != EOF; ) 
    switch (option) {
    default :
      badCommandLine=true;
      break;
    case 'n' :
      StartConfig = std::atoi(optarg);
      std::cout << "Starting No. of configurations: " << StartConfig << std::endl;
      break;
    case 'N' :
      EndConfig = std::atoi(optarg);
      std::cout << "Ending No. of configurations: " << EndConfig << std::endl;
      break;
    case 's' :
      searchR = std::atof(optarg);
      std::cout << "Global support size: " << searchR << std::endl;
      break;
    case 'b' :
      beta = std::atof(optarg);
      std::cout << "Beta: " << beta << std::endl;
      break;
    case 'a' :
      avg = std::atoi(optarg);
      std::cout << "Average of "<< avg << " trajectory points taken" << std::endl;
      break;
    case 't' :
      tol = std::atof(optarg);
      std::cout << "LME tolerance " << tol << std::endl;
      break;
    case 'i' :
      nItMax = std::atoi(optarg);
      std::cout << "LME max iterations " << nItMax << std::endl;
      break;
    }

  if( badCommandLine || argc < 2 ) {
    cout << "Usage: homog inputFileName [-n start_configs -N end_config -s support_size - b LME beta -a average_number -t LME tol -i LME max iter]." << endl;
    return(0);
  }
  


  // Create input stream
  ifs.open( inputFileName.c_str(), ios::in);
  if (!ifs) {
    cout << "Cannot open input file: " << inputFileName << endl;
    exit(0);
  }

  // Create vector of nodes
    int dof = 0;
    std::vector< DeformationNode<3>* > Nodes;
    vector<double> node_volume, supp_size;
    double Gsupp = 0.0;

    // Input .dat file containing nodes and weights
    std::string token;
    ifs >> token;
    while( token != "POINTS" ) ifs >> token;
    int npts=0;
    ifs >> npts; 
    Nodes.reserve(npts);
  
    // Read in points
    for(int i = 0; i < npts; i++) {
      int id = i;
      int temp;
      double weight, supp;

      DeformationNode<3>::Point x;
      ifs >> temp >> x(0) >> x(1) >> x(2) >> weight >> supp;
      // cout << x(0) << " " << x(1) << " " << x(2) << " " << weight << " " << supp << endl;
      NodeBase::DofIndexMap idx(3);
      for(int j=0; j<3; j++) idx[j] = dof++;
      DeformationNode<3>* n = new DeformationNode<3>(id, idx, x);

      Nodes.push_back( n );

      node_volume.push_back(weight);
      supp_size.push_back(supp);
      Gsupp += supp;
    }

    Gsupp /= npts;
    beta /= pow(Gsupp,2.0);
  cout << "Number of nodes: " << Nodes.size() << endl;
  cout << "Avg supp size: " << Gsupp << endl;
  cout << "LME const beta used: " << beta << endl;
  ifs.close();



  // Create material, values don't matter here
  double rho = 1.0;
  double E = 200.0;
  double nu = 0.4;
  double viscosity = 1.0e-2;
  
  // Create and initialize body
  typedef CompNeoHookean MaterialType;
  MaterialType protein( rho, E, nu );
  
  // create Body
  typedef LMEbody<MaterialType, LMEshape> body;
  body LMEbd(protein, Nodes, node_volume, supp_size, beta, searchR, tol, nItMax);
 
  // LMEbd.setOutput(paraview);
  // LMEbd.compute(true,true,false);
  
  std::cout << std::endl;
  std::cout << "Body energy = " << LMEbd.energy() << std::endl;
  std::cout << "Body volume = " << LMEbd.volume() << std::endl;
  std::cout << std::endl;

  // Open the output files
  ofstream file1, file2, file3;
  file1.open("inv1.dat",ios::binary);
  file2.open("inv2.dat",ios::binary);
  file3.open("inv3.dat",ios::binary);
  char tab = '\t', newline = '\n';
  int col = npts+1, row = EndConfig-StartConfig+1;

  // Write header file
  file1.write((char*)&row,sizeof(int)); file1.write((char*)&col,sizeof(int));
  file2.write((char*)&row,sizeof(int)); file2.write((char*)&col,sizeof(int));
  file3.write((char*)&row,sizeof(int)); file3.write((char*)&col,sizeof(int));

  // Read in the configuration files in a loop
  for(int i=StartConfig;i<=EndConfig;i++)
  {
    char inpfile[20];
    std::cout<< "Reading " << i << " th configuration." << std::endl;
    float x[npts], y[npts], z[npts];
    if(datfile){
      // sprintf(inpfile,"../../com-traj/%d.dat",i);
      sprintf(inpfile,"/u/home/cardio/ankush/research-new/Homog/semv/com-traj/%d.dat", i);
      // sprintf(inpfile,"./TestTraj/%d.dat", i);
      readconfig_dat(inpfile,x,y,z,npts,avg);
    }
    
    if(pdbfile){
      sprintf(inpfile,"../calpha-traj/%d.pdb",i);
      readconfig_pdb(inpfile,x,y,z,npts);
    }
 
    // std::cout<< "File has been read" << std::endl;

    for(int inode = 0; inode < npts; inode++)
    {
      DeformationNode<3>::Point newx;
      newx(0) = x[inode]; newx(1) = y[inode]; newx(2) = z[inode];
      Nodes[inode]->setPoint(newx);
    }

    // Compute the invariants
    std::vector<double> I1(npts,0.), I2(npts,0.), I3(npts,0.);
    // std::cout<< "Calculating the invariants" << std::endl;
    LMEbd.cal_invariants(I1,I2,I3);
 
    // Print the results 
    double di = double(i);
    // std::cout<<"double of i is "<<di<<std::endl;
    file1.write((char*)&di,sizeof(double));
    file2.write((char*)&di,sizeof(double));
    file3.write((char*)&di,sizeof(double));
    
    for(int inode=0;inode<npts;inode++)
    {
      file1.write((char*)&I1[inode],sizeof(double)); 
      file2.write((char*)&I2[inode],sizeof(double));
      file3.write((char*)&I3[inode],sizeof(double));
      
      // cout << Nodes[inode]->position() << endl;
      // cout << I1[inode] << " " <<  I2[inode] << " " <<  I3[inode] << " " << endl;
    }
    
  }
  
  std::cout<<"Everything done!"<<std::endl;
  return 0;

}



void readconfig_dat(char* filename, float* x, float* y, float* z, int npts){
  ifstream inpfile;
  inpfile.open(filename);
  if (!inpfile) {
    cout << "Cannot open input file: " << filename << endl;
    exit(0);
  }
  int id;
  
  for(int i=0; i<npts; i++){
    inpfile >> id >> x[i] >> y[i] >> z[i];
  }
  inpfile.close();
}



void readconfig_pdb(char* filename, float* x, float* y, float* z, int npts){
  ifstream inpfile;
  inpfile.open(filename);
  if (!inpfile) {
    cout << "Cannot open input file: " << filename << endl;
    exit(0);
  }
  int i=0;
  string line;
  
  while (! inpfile.eof() )
    {
      getline (inpfile,line);
      if (line.length()>0 && line.substr(0,4)=="ATOM"){
	x[i] = atof(line.substr(30,8).c_str()); 
	y[i] = atof(line.substr(38,8).c_str()); 
	z[i] = atof(line.substr(46,8).c_str()); 
	i++;
      }
    }
	inpfile.close();
}



void readconfig_dat(char* filename, float* x, float* y, float* z, int npts, int average)
{
  // Average the N-closest points //
  ifstream inpfile;
  inpfile.open(filename);
  if (!inpfile) {
    cout << "Cannot open input file: " << filename << endl;
    exit(0);
  }
  int id;
  
  for(int i = 0; i < npts; i++)
  {
    // Average nodal positions
    double xm = 0.0, ym = 0.0, zm = 0.0;
    for(int j = 0; j < average; j++)
    {
      double xj, yj, zj;
      inpfile >> id >> xj >> yj >> zj;
      xm += xj; ym += yj; zm += zj;
    }
    xm /= average; ym /= average; zm /= average;
    x[i] = xm; y[i] = ym; z[i] = zm;
  }
  inpfile.close();
}
