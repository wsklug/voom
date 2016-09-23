// This file reads in the configurations from MD
// either in .dat or .pdf formats and the triangulation of the domain
// Then calculates the deformation gradient at the barycenter of each element
// Then calculates the three invariants in each element and writes them to the
//  file invariants1.dat *2.dat and *3.dat
// in the format : at each configuration a row of invariants per each element
#include <iostream>
#include <vector>
#include <fstream>
#include <getopt.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace std;
#include "voom.h"
#include "Body3D.h"
#include "NodeBase.h"
#include "Node.h"
#include "Model.h"
#include "Node.h"
#include "CompNeoHookean.h"
#include "StVenant.h"
#include "TetQuadrature.h"
#include "ShapeTet4CP.h"
#include "Lbfgsb.h"
#include <tvmet/Vector.h>
#include <fstream>
// #include "CompNeoHookean.h"
#include "HomogMP.h"
#include<tvmet/Vector.h>
#include<tvmet/Matrix.h>

using namespace voom;
using namespace std;

void readconfig_pdb(char* filename, float* x, float* y, float* z, int npts);
void readconfig_dat(char* filename, float* x, float* y, float* z, int npts);
void readconfig_dat(char* filename, float* x, float* y, float* z, int npts, int average);
double most_prob(blitz::Array<double,1> x);
double most_probJ(blitz::Array<double,1> x);

int main(int argc, char* argv[])
{
  if( argc < 2 ) {
    cout << "Usage: homog [-M mesh_file -n start_config -N end_config -a average_number]." << endl;
    return(0);
  }
  
  bool verbose=true;
  
  for(int i=0; i<argc; i++) {
    cout << setw(8) << i << "\t"
	 << argv[i] << endl;
  }
  
  //Input file name and other parameters from run command
  ifstream ifsV;
  string meshFileName;
  bool datfile = true;
  bool pdbfile = false;
  bool badCommandLine = false;
  int StartConfig = 0, EndConfig = 0, avg = 1;
  double k = -1.0;
 
  
  for(char option; (option=getopt(argc,argv,"M:n:N:a:k:")) != EOF; ) 
    switch (option) {
    default :
      badCommandLine=true;
      break;
    case 'M' :
      meshFileName = string(optarg);
      cout << "Mesh file name: " << meshFileName << endl;
      break;
    case 'n' :
      StartConfig = atoi(optarg);
      cout << "Starting No. of configurations: " << StartConfig << endl;
      break;
    case 'N' :
      EndConfig = atoi(optarg);
      cout << "Ending No. of configurations: " << EndConfig << endl;
      break;
    case 'a' :
      avg = atoi(optarg);
      cout << "Average of "<< avg << " trajectory points taken" << endl;
      break;
    case 'k' :
      k = atof(optarg);
      cout << "Minimization parameter k = "<< k << endl;
      break;
    }

  if( badCommandLine || argc < 2 ) {
    cout << "Usage: homog [-M mesh_file -n start_configs -N end_config -a average_number -k minimization parameter]." << endl;
    return(0);
  }




  
  // Create connectivity table
  ifsV.open( meshFileName.c_str(), ios::in);
  if (!ifsV) {
    cout << "Cannot open mesh file: " << meshFileName << endl;
    exit(0);
  }

  // Create vector of nodes and connectivities
  int dof = 0;
  vector< NodeBase* > nodes;
  vector< DeformationNode<3>* > DefNodes;
  vector<vector< int> > Connectivities;
  
  // Input .vtk file containing nodes and connectivities
  string token;
  ifsV >> token; 
  while( token != "POINTS" ) ifsV >> token;

  int npts = 0;
  ifsV >> npts; 
  DefNodes.reserve(npts);
  nodes.reserve(npts);

  ifsV >> token;   // skip number type
  cout<< "npts = " << npts << endl;

  // Read in points - Position will be adjusted according to inputFile
  for(int i = 0; i < npts; i++) {
    int id = i;
    DeformationNode<3>::Point x;
    ifsV >> x(0) >> x(1) >> x(2);
    NodeBase::DofIndexMap idx(3);
    for(int j = 0; j < 3; j++) idx[j]=dof++;
    DeformationNode<3>* n = new DeformationNode<3>(id,idx,x);
    nodes.push_back( n );
    DefNodes.push_back( n );
  }
  assert(DefNodes.size()!=0);


  vector<Material * > Mat;
  // Read in tet connectivities
  while( token != "CELLS" ) ifsV >> token;
  int ntet = 0; 
  ifsV >> ntet;
  Connectivities.reserve(ntet);
  vector<int > c(4,0);

  cout << "Number of tetrahedra: " << ntet << endl;
  Mat.reserve(ntet);
  int tmp; 
  ifsV >> tmp;
  // cout << " tmp = " << tmp << endl;
  for (int i = 0; i < ntet; i++){
    int tmp = 0;
    ifsV >> tmp;
    if(tmp != 4) { 
      std::cout<<"Some mistake reading the elements connectivity from file. Check again."<<std::endl;
      exit(0);
    }
    for(int a = 0; a < 4; a++) ifsV >> c[a];
    Connectivities.push_back(c);
    Material * mat = new HomogMP();
    Mat.push_back(mat);
  }
  cout << "Connectivity table has been read " << endl;

  ifsV.close();



  // Create quadrature rule and shape function
  TetQuadrature Quad(1);
  ShapeTet4 Sh;
 
  // Create body
  Body3D BodyHomog(Mat, Connectivities, DefNodes, Quad, Sh, k);
  cout << "Body has been initialized once and for all :) " << endl;


  /*
  // Creating the model to solve for the minima
  Model::BodyContainer bdc;
  bdc.push_back(&BodyHomog);
  Model model(bdc, nodes);

  //model.checkConsistency(true,false);
  //model.checkRank(model.dof()-6,true);
  //return 0;

  //Input lbfgsb solver parameters
  int m = 5;
  double factr = 1.0e+1;
  double pgtol = 1.0e-5;
  int iprint = 0;
  ifstream lbfgsbinp("lbfgsb.inp");
  lbfgsbinp >> iprint >> factr >> pgtol >> m ;
  if(verbose) 
    cout << "Input iprint: " << iprint << endl
	 << "Input factr:  " << factr  << endl
	 << "Input pgtol:  " << pgtol  << endl
	 << "Input m:      " << m      << endl;
  Lbfgsb solver(model.dof(), m, factr, pgtol, iprint ); //(true);

  // add a viscous regularizer
  // double viscosity = 10.;
  // ViscousRegularizer vr(bd.nodes(), viscosity);
  // BodyHomog.addElement( &vr ); 



  int f = 0;
  // Initial invariants and QP positions //
  vector<pair<Vector3D, vector<double > > > Initial;
  Initial = BodyHomog.invariants(f);
  cout << "Elements with detF < 0.0 = " << f << endl;
  f = 0;
  ofstream fileRef;
  fileRef.open("Reference.dat");

  // fileRef << "x y z I1 I2 I3 V" << endl;
  Body::ElementContainer TetEl = BodyHomog.elements();
  double TotVol = 0.0, TetVol = 0.0;
  for(int el = 0; el < ntet; el++)
  {
    Vector3D location = Initial[el].first;
    vector<double > inv = Initial[el].second;
    TetVol = dynamic_cast<Element3D *>(TetEl[el])->volume();

    TotVol += TetVol;
    
    fileRef << location[0] << " " << location[1] << " " << location[2] << " " 
	    << inv[0]      << " " << inv[1]      << " " << inv [2]     << " "
	    << TetVol << endl;
    
  }
  fileRef.close();

  // cout << "Initial invariants printed " << endl;
  cout << "Body volume from file = " << TotVol << endl;



// Start minimization loop
 int Iteration = 0;
 double EnergyError = 1.0, VolIt = 0.0, EnIt = 0.0;
 int col = ntet+1, row = EndConfig-StartConfig+1;

// Needed for minimization
  blitz::Array<double,2> I1bar(row,col-1);
  blitz::Array<double,2> I2bar(row,col-1);
  blitz::Array<double,2> J(row,col-1);
  I1bar=0.0; I2bar=0.0; J=0.0;

 vector<double > EnergyIt(3, 0.0);
 ofstream fileMinHist;
 fileMinHist.open("MinimizationHistory.dat");

while (Iteration < 100 && EnergyError > 0.001)
{ 
  BodyHomog.compute(true,false,false);
  EnIt = BodyHomog.totalStrainEnergy();
  VolIt = BodyHomog.volume();
  EnergyIt[0] = EnIt;
  std::cout << "Body strain energy = " << EnIt << std::endl;
  std::cout << "Body volume = " <<  VolIt << std::endl;
  std::cout << std::endl;


  // Read in the configuration files in a loop
  DeformationNode<3>::Point newx;
  int irow = 0;
  for(int i = StartConfig; i <= EndConfig; i++)
  {
    char inpfile[20];
    if (i%1000 == 0) cout<< "Reading " << i << " th configuration." << std::endl;
    float x[npts], y[npts], z[npts];
    if(datfile){
      // sprintf(inpfile,"../../com-traj/%d.dat",i);
      sprintf(inpfile,"/u/home/cardio/ankush/research-new/Homog/semv/com-traj/%d.dat", i);
      // printf(inpfile,"/u/home/cardio/ankush/research-new/Homog/semv/com-traj/%d.dat", i);
      // sprintf(inpfile,"./TestTraj/%d.dat", i);
      readconfig_dat(inpfile,x,y,z,npts,avg);
    }
    
    if(pdbfile){
      sprintf(inpfile,"../calpha-traj/%d.pdb",i);
      readconfig_pdb(inpfile,x,y,z,npts);
    }
    
    for(int inode = 0; inode < npts; inode++)
    {
      newx(0) = x[inode]; newx(1) = y[inode]; newx(2) = z[inode];
      DefNodes[inode]->setPoint(newx);
    }

    // Compute the invariants
    std::vector<double> I1(ntet, 0.0), I2(ntet, 0.0), I3(ntet, 0.0);
    f = BodyHomog.invariants(I1,I2,I3);
    // cout << "Elements with detF < 0.0 " << f << " over " << ntet << " .In percentage = " << double(f)/double(ntet) <<  endl;
    
    for(int el = 0; el < ntet; el++)
    {
      I1bar(irow, el) = I1[el];
      I2bar(irow, el) = I2[el];
      J(irow, el) = I3[el];
    }
    irow++;
  }


  //-------------------//
  // calculates a functional = (sum over nodes) (I1bar_prob-3)^2 where I1bar_prob is the most probable I1bar
  vector<double> I1MostProb(ntet, 3.0), I2MostProb(ntet, 3.0), JMostProb(ntet, 1.0);

  for(int el = 0; el < ntet; el++)
  {
    // Call the function to get most probable invariant 
    I1MostProb[el] = most_prob(I1bar(blitz::Range::all(),el)); // cout << I1MostProb[el] << endl;
    I2MostProb[el] = most_prob(I2bar(blitz::Range::all(),el)); // cout << I2MostProb[el] << endl;
    JMostProb[el] = most_probJ(J(blitz::Range::all(),el));     // cout << JMostProb[el] << endl;
  }

  // cout << "Calculated the most probable strain state invariants" << endl;
  BodyHomog.setMPinv(I1MostProb,I2MostProb,JMostProb);

  // reset the reference state
  for(int inode = 0; inode < npts; inode++)
  {
    DefNodes[inode]->setPoint(DefNodes[inode]->position());
  }

  BodyHomog.compute(true,false,false);
  EnIt = BodyHomog.energy();
  EnergyIt[1] = EnIt;
  cout << endl 
       << "Body energy = "
       << EnIt << endl;





  // solve to minimize the functional
  solver.solve(&model);


  BodyHomog.compute(true,false,false);
  EnIt = BodyHomog.energy();
  EnergyIt[2] = EnIt;
  cout << "Body energy = " << EnIt << endl
       << " -------------------------------- " 
       << endl << endl;


  // Set the new reference configuration
  for(int inode = 0; inode < npts; inode++)
  {
    DefNodes[inode]->resetPosition();
  }
  BodyHomog.reset();



  Iteration++;
  EnergyError = (EnergyIt[1]-EnergyIt[2])/EnergyIt[2];
   
  fileMinHist << VolIt   << " " 
	      << EnergyIt[0] << " " 
	      << EnergyIt[1] << " " 
	      << EnergyIt[2] << endl;
 }

 fileMinHist.close();
*/
   int col = ntet+1, row = EndConfig-StartConfig+1, f = 0;

cout << "Read in files for the last time " << endl;
// Open the output files
  ofstream file1, file2, file3;
  file1.open("inv1.dat",ios::binary);
  file2.open("inv2.dat",ios::binary);
  file3.open("inv3.dat",ios::binary);
  char tab = '\t', newline = '\n';

  // Write header file
  file1.write((char*)&row,sizeof(int)); file1.write((char*)&col,sizeof(int));
  file2.write((char*)&row,sizeof(int)); file2.write((char*)&col,sizeof(int));
  file3.write((char*)&row,sizeof(int)); file3.write((char*)&col,sizeof(int));

 // Read in the configuration files in a loop
  DeformationNode<3>::Point newx;
  int irow = 0;
  for(int i = StartConfig; i <= EndConfig; i++)
  {
    char inpfile[20];
    if (i%1000 == 0) cout<< "Reading " << i << " th configuration." << std::endl;
    float x[npts], y[npts], z[npts];
    if(datfile){
      // sprintf(inpfile,"../../com-traj/%d.dat",i);
      // sprintf(inpfile,"/u/home/cardio/ankush/research-new/Homog/semv/com-traj/%d.dat", i);
      sprintf(inpfile,"/u/home/cardio/ankush/research-new/namd-runs/t4-enm/allmode-trajectory/%d.dat", i);
      // sprintf(inpfile,"./TestTraj/%d.dat", i);
      readconfig_dat(inpfile,x,y,z,npts,avg);
    }
    
    if(pdbfile){
      sprintf(inpfile,"../calpha-traj/%d.pdb",i);
      readconfig_pdb(inpfile,x,y,z,npts);
    }
    
    for(int inode = 0; inode < npts; inode++)
    {
      newx(0) = x[inode]; newx(1) = y[inode]; newx(2) = z[inode];
      DefNodes[inode]->setPoint(newx);
    }

    // Compute the invariants
    std::vector<double> I1(ntet, 0.0), I2(ntet, 0.0), I3(ntet, 0.0);
    f = BodyHomog.invariants(I1,I2,I3);
    // cout << "Elements with detF < 0.0 " << f << " over " << ntet << " .In percentage = " << double(f)/double(ntet) <<  endl;
    
    // Print the results 
    double di = double(i);
    file1.write((char*)&di,sizeof(double));
    file2.write((char*)&di,sizeof(double));
    file3.write((char*)&di,sizeof(double));
    for(int el = 0; el < ntet; el++)
    {
      file1.write((char*)&I1[el],sizeof(double)); 
      file2.write((char*)&I2[el],sizeof(double));
      file3.write((char*)&I3[el],sizeof(double));
    }
  }

  file1.close();
  file2.close();
  file3.close();



  std::cout<<"Everything done!"<<std::endl;
  


  for (int i = 0; i < npts; i++) 
    delete nodes[i];
  for (int i = 0; i < ntet; i++) 
    delete Mat[i];

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


double most_prob(blitz::Array<double,1> x){
  //find max, min and size of vector x (cut the limits at 3 and 6)
  double max,min;
  max=blitz::max(x(blitz::Range::all()));
  min=blitz::min(x(blitz::Range::all()));
  if(max>min+3.) {max=min+3.;}//this is only for I1bar, to prevent bad binning due to one large value of I1bar
  int size=x.rows();
  //create the bins
  int nbins=30;
  double interval=(max-min)/nbins;
  std::vector<double> bin_centers(nbins,0.);
  std::vector<int> counts(nbins,0);
 
  for(int i=0;i<nbins;i++){
    bin_centers[i] = min+(i+0.5)*interval;
  }

  //count the frequency in them
  for(int i=0;i<size;i++){
    int bin=int((x(i)-min)/interval);
    if(bin>nbins-1 || bin < 0) continue; //bin=nbins-1;
    counts[bin]++;
  }

  //find the bin with maximum counts
  int maxcount=0; int maxid=0;
  for(int i=0;i<nbins;i++){
    if(counts[i]>maxcount){
      maxcount=counts[i];
      maxid=i;
    }
  }
  //return the bin of maxid i.e. the most probable x
  return bin_centers[maxid];
}

double most_probJ(blitz::Array<double,1> x){
  //find max, min and size of vector x (cut the limits at 0.5 and 2)
  double max,min;
  max=blitz::max(x(blitz::Range::all()));
  min=blitz::min(x(blitz::Range::all()));
  if(max>2.) max=2.; 
  if(min<0.5) min=0.5; //this is only for J, to prevent bad binning due to one large or one small value of J
  int size=x.rows();
  //create the bins
  int nbins=30;
  double interval=(max-min)/nbins;
  std::vector<double> bin_centers(nbins,0.);
  std::vector<int> counts(nbins,0);
 
  for(int i=0;i<nbins;i++){
    bin_centers[i] = min+(i+0.5)*interval;
  }

  //count the frequency in them
  for(int i=0;i<size;i++){
    int bin=int((x(i)-min)/interval);
    if(bin>nbins-1 || bin < 0) continue; //bin=nbins-1;
    counts[bin]++;
  }

  //find the bin with maximum counts
  int maxcount=0; int maxid=0;
  for(int i=0;i<nbins;i++){
    if(counts[i]>maxcount){
      maxcount=counts[i];
      maxid=i;
    }
  }
  //return the bin of maxid i.e. the most probable x
  return bin_centers[maxid];
} 
