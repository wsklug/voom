#include <string>
#include <iostream>
#include <vector>
#include "Node.h"
#include "SCElastic.h"
#include "EvansElastic.h"
#include <tvmet/Vector.h>
#include <fstream>
#include "LoopShellBody.h"
#include "Model.h"
#include "Solver.h"
#include "ConjugateGradientWSK.h"
#include "SimulatedAnnealing.h"
#include "Lbfgsb.h"

//#define NODENUMBER 602
// #define ELEMNUMBER 80

using namespace tvmet;
using namespace std;
using namespace voom;

void ioSetting(int argc, char* argv[], ifstream&, string&);



class PenaltyBC : public Element
{
public:

  // typedefs
  typedef std::vector<NodeBase*> BaseNodeContainer;
  typedef std::vector<NodeBase*>::iterator BaseNodeIterator;
  typedef std::vector<NodeBase*>::const_iterator ConstBaseNodeIterator;
  typedef tvmet::Vector<double,3> Vector3D;
  typedef tvmet::Vector<bool,3> VectorBC;
  typedef DeformationNode<3> DefNode;

  PenaltyBC( DefNode * node, Vector3D x0, VectorBC bc, double k ) {
    _defNode = node;
    _baseNodes.push_back(node);
    _k = k;
    _x0 = x0;
    _bc = bc;
  }

  // E = 0.5 * k * ( x3-Z )^2
  // f3 = dE/dx3 = k * ( x3-Z )
  //! Do mechanics on Body
  virtual void compute( bool f0, bool f1, bool f2 ) {
    //std::cout << "RigidBC::compute()" << std::endl;

    if(f0) _energy = 0.0;

    Vector3D dx(0.0);
    for(int i=0; i<3; i++) {
      if(_bc(i)) dx(i) = _defNode->getPoint(i) - _x0(i);
    }
    
    if(f0) {
      _energy += 0.5*_k*tvmet::dot(dx,dx);
    }
    if(f1) {
      for(int i=0; i<3; i++) {
	if(_bc(i)) _defNode->addForce(i, _k*dx(i));
      }
    }
  }
    
private:

  double _k;
  Vector3D _x0;
  VectorBC _bc;
  DefNode * _defNode;
  BaseNodeContainer _baseNodes;
  };


class PointLoad : public Element
{
public:

  // typedefs
  typedef std::vector<NodeBase*> BaseNodeContainer;
  typedef std::vector<NodeBase*>::iterator BaseNodeIterator;
  typedef std::vector<NodeBase*>::const_iterator ConstBaseNodeIterator;
  typedef tvmet::Vector<double,3> Vector3D;
  typedef DeformationNode<3> DefNode;

  PointLoad( DefNode * node, Vector3D f ) {
    _defNode = node;
    _baseNodes.push_back(node);
    _force = f;
  }

  virtual void compute( bool f0, bool f1, bool f2 ) {

    if(f0) _energy = 0.0;
    
    if(f0) {
      for(int i=0; i<3; i++)
	_energy += -_force(i)*_defNode->getPoint(i);
    }
    if(f1) {
      for(int i=0; i<3; i++) {
	_defNode->addForce(i, -_force(i));
      }
    }
  }

  void setForce(Vector3D & nf) {
//     for(int i=0; i<3; i++) {
//       _force(i)=nf(i);
//     }
    _force = nf;
  }
    
private:

  Vector3D _force;
  DefNode * _defNode;
  BaseNodeContainer _baseNodes;
  };


int main(int argc, char* argv[])
{
         
  ifstream ifs;
  string ofn;
 
  ioSetting(argc, argv, ifs, ofn);

	
  
  int NODENUMBER, ELEMNUMBER;
	
  ifs >> NODENUMBER >> ELEMNUMBER;
	
  //
  // create vector of nodes
  int dof=0;
  double alpha =1.2;
  double beta  =0.5;
  double Radius=0.53;
  double x0, x1, x2;
  std::vector< NodeBase* > nodes;
  for ( int i = 0; i < NODENUMBER; i++){
    int id=-1;
    DeformationNode<3>::Point x;
    ifs >> id >> x(0) >> x(1) >> x(2);
    //x(0)*=0.8;

    NodeBase::DofIndexMap idx(3);

    
   for(int j=0; j<3; j++) idx[j]=dof++;
   nodes.push_back(new DeformationNode<3>(id,idx,x));    

}
    

  //
  // create connectivities
  vector< tvmet::Vector<int,3> > connectivities;
  tvmet::Vector<int, 3> c;
  for (int i = 0; i < ELEMNUMBER; i++){
    ifs >> c[0];c[0]-=1;
    ifs >> c[2];c[2]-=1;
    ifs >> c[1];c[1]-=1;

    connectivities.push_back(c);
  }

  ifs.close();
  /*  
  
          
  ifstream ifs; 
  string ofn;
 
  ioSetting(argc, argv, ifs, ofn);
  

  //
  // create vector of nodes
  int dof=0;
  std::vector< NodeBase* > nodes;
  char key;
  ifs>>key;
  double radius=0.0;
  //ofstream outfile("volkmarn.obj");
  for ( int i = 0; key=='v'; i++, ifs>>key){
    int id=i;
    DeformationNode<3>::Point x;
    ifs >> x(0) >> x(1) >> x(2);
    //x(2)*=1.45;

    //x(0)*=1.75;
    //x(1)*=1.75;   
    
    NodeBase::DofIndexMap idx(3);
    for(int j=0; j<3; j++) idx[j]=dof++;
    nodes.push_back(new DeformationNode<3>(id,idx,x));
  }
  
  // constrain top and bottom nodes

  //
  // create connectivities
  vector< tvmet::Vector<int,3> > connectivities;
  tvmet::Vector<int, 3> c;
  for (int i = 0; key=='f'; i++){
    int tmp=0;
    ifs >> tmp; c[1]=tmp-1;
    ifs >> tmp; c[0]=tmp-1;
    ifs >> tmp; c[2]=tmp-1;

    connectivities.push_back(c);
    if(!(ifs>>key)) break;
  }
  cout << connectivities.size() << endl;
  
  ifs.close();
  */
        
  
  double KC = 1.0e0;
  double KG = 0.0e0;
  double C0 = 0.0;
  double mu = 3.5e0;
  double KS = 2.0e0;
  double kSpring = 0.0;
  double viscosity = 0.0e0;
  

  // create Body
  int nBoundaries=0;
  double pressure=0.0;
  double tension=0.0;
  double totalCurvatureForce=0.0;
  int quadOrder=1;
  double penaltyVolume=1.0e4;
  double penaltyArea=1.0e4;
  double penaltyTotalCurvature=1.0e4;
  
  ifstream bdinp("constants.inp");
  bdinp >> mu >> KS >>penaltyVolume >> penaltyArea >> penaltyTotalCurvature >> kSpring >> viscosity;
  std::cout << "mu ="  << mu <<  std::endl
            << "KS ="  << KS <<  std::endl
            << "pV ="  << penaltyVolume << std::endl
            << "pA ="  << penaltyArea   << std::endl
	    << "pTC="  << penaltyTotalCurvature<< std::endl
	    << "kSpring ="<< kSpring    << std::endl
            << "viscosity =" << viscosity  << std::endl;

  EvansElastic bilayer( KC, KG, C0, mu, KS, kSpring );
  //SCElastic bilayer(KC, KG, C0);
  std::cout << "EvansElastic Material has been created." << std::endl;
  //std::cout << "SCElastic Material has been created." << std::endl;	
  //

  LoopShellBody<EvansElastic> bd(bilayer, connectivities, nodes, quadOrder, 
				 nBoundaries, pressure, tension, totalCurvatureForce, penaltyVolume, penaltyArea, penaltyTotalCurvature, viscosity); 

  //LoopShellBody<SCElastic> bd(bilayer, connectivities, nodes, ngp, 
  //		      nBoundaries, pressure, penaltyVolume, penaltyArea);
  cout << "Created a body." << endl
       << "  # nodes: " << bd.nodes().size()<<std::endl
       << "  # shell elements: " << bd.elements().size() << std::endl;
  
  bd.setOutput(Body::paraview);
  bd.compute(false,false,false);
  
//   srand(time(0));	       
//   for(int i=0; i<NODENUMBER; i++) {
//   for(int j=0; j<3; j++) 
//    nodes[i]->addPoint(j,0.01*((double)(rand())/RAND_MAX-0.5) );
//   }	
  
  
  /*double Zmin=std::numeric_limits<double>::max();
  double Zmax=-std::numeric_limits<double>::max();
  cout << nodes.size() << endl;
  for(int i=0; i<nodes.size(); i++) {
    Zmin = std::min(Zmin,nodes[i]->getPoint(2));
    Zmax = std::max(Zmax,nodes[i]->getPoint(2));
  }

  int nn = 24;
  
  tvmet::Vector<int,12> connectT;
  connectT(0)=1;
  connectT(1)=51;
  connectT(2)=101;
  connectT(3)=151;
  connectT(4)=201;
  connectT(5)=251;
  connectT(6)=301;
  connectT(7)=351;
  connectT(8)=401;
  connectT(9)=451;
  connectT(10)=501;
  connectT(11)=551;

  tvmet::Vector<int,12> connectB;  
  connectB(0)=50;
  connectB(1)=100;
  connectB(2)=150;
  connectB(3)=200;
  connectB(4)=250;
  connectB(5)=300;
  connectB(6)=350;
  connectB(7)=400;
  connectB(8)=450;
  connectB(9)=500;
  connectB(10)=550;
  connectB(11)=600;  
  

  tvmet::Vector<PointLoad*,24> neighborTopLoad(0);
  tvmet::Vector<PointLoad*,24> neighborBottomLoad(0);  

  cout<<"size"<<nodes.size()<<endl;

  //double f2 = 1.0e3;
  //bdinp >> f2;
  PointLoad * topLoad = 0;
  PointLoad * bottomLoad = 0;

  tvmet::Vector<double,3> f(0.0);
  
  for(int i=0; i<nodes.size(); i++) {
    double Z = nodes[i]->getPoint(2);
    DeformationNode<3> * point = dynamic_cast< DeformationNode<3>* >(nodes[i]);
    if( false point != 0 ) {

      if( Z == Zmax ) {        
	std::cout << "Loading node at Z = " << Z << std::endl;
	//f(2)=f2;				
	topLoad = 
	  new PointLoad( point, f );
	bd.pushBack( topLoad );        
      }
     
      if( Z == Zmin ) {
	std::cout << "Fixing node at Z = " << Z << std::endl;
	//f(2)=-f2;
	bottomLoad = 
	  new PointLoad( point, f );
	bd.pushBack( bottomLoad );        	
      }
    
  
      
      for(int j=0; j<nn; j++){
	if( i==connectT(j) ){
	  neighborTopLoad(j) =
	    new PointLoad( point, f );
	  bd.pushBack( neighborTopLoad(j) ); 
	}
      }
      for(int j=0; j<nn; j++){
	if( i==connectB(j) ){
	  neighborBottomLoad(j) =
	    new PointLoad( point, f );
	  bd.pushBack( neighborBottomLoad(j) ); 
	}
      }
    }
  }*/

  //topLoad->checkConsistency();
  std::cout << "  # all elements: " << bd.elements().size() << std::endl;

  //
  // create Body Container
  Model::BodyContainer bdc;
  bdc.push_back(&bd);

  //
  // create Model
  Model model(bdc);

  //model.checkConsistency(true,false);
  //return 0;

  int m=5;
  double factr=1.0e+7;
  double pgtol=1.0e-5;
  int iprint = -1;
  double pentol=1.0e-4;
  //bool boundaryConditionFlag = true;

  //ifstream lbfgsbinp("lbfgsb.inp");
  //lbfgsbinp >> iprint >> factr >> pgtol >> m >> pentol;
  //if(verbose) 
    std::cout << "Input iprint: " << iprint << std::endl
	      << "Input factr: " << factr << std::endl
	      << "Input pgtol: " << pgtol << std::endl
	      << "Input m: " << m << std::endl
	      << "Input pentol: " << pentol << std::endl;
  Lbfgsb solver(model.dof(), m, factr, pgtol, iprint );//(true);

  
  //model.checkConsistency(true,false);
//   model.checkRank(model.dof()-6,true);
  //return 0;
/*   
  ConjugateGradientWSK CGsolver;//(true);
  int maxIter = 1000*model.dof();
  int restartStride = 2*model.dof();
  int printStride = 200;//10*restartStride;
  double tol = 1.0e-6;
  double absTol = 1.0e-5;
  double tolLS = 1.0e-6;
  int maxIterLS = 20;
    
  CGsolver.setParameters(voom::ConjugateGradientWSK::Secant,
			 voom::ConjugateGradientWSK::PR,
			 maxIter, restartStride, printStride, 
			 tol, absTol, tolLS, maxIterLS);
  CGsolver.setWolfeParameters(1.0e-6, 1.0e-3, 1.0e-4);

  //CGsolver.zeroOutData(true, false, false);
  */
  model.print("BeforeCG");
  
  double v = bd.volume();
  double a = bd.area(); 
  double tc = bd.totalCurvature();
  double vred = 6.0*sqrt(M_PI)*v/std::pow(a,3.0/2.0);
  double V = bd.constraintVolume(); 
  double A = bd.constraintArea();
  double TC = bd.constraintTotalCurvature();
  double Vred = 6.0*sqrt(M_PI)*V/std::pow(A,3.0/2.0);
  std::cout << "Before CG:"<< std::endl
	    << "Volume = "<< bd.volume() << std::endl
	    << "Cons. Volume = "<< bd.constraintVolume() << std::endl
	    << "Area = "<< bd.area() << std::endl
	    << "Cons. Area = "<< bd.constraintArea() << std::endl
	    << "Reduced Volume = " << vred << std::endl
	    << "Cons. Reduced Volume = " << Vred << std::endl
	    << "Total Curvature= " << tc <<std::endl
	    << "Cons. Total Curvature= " << TC <<std::endl
	    << "Energy = " << solver.function() << std::endl
            <<  "ConstraintEnergy = " << bd.constraintEnergy() << std::endl
	    <<  "strainEnergy = " << bd.totalStrainEnergy() << std::endl;
    //	    <<  "workenergy = " << bd.work() <<std::endl;
    
  //model.checkRank(model.dof()-6,true);
  //  for(double nu=0.615424/*0.82127*/; nu > 0.6; nu-=0.1) {
//     bd.reduceVolume(nu); 


  //model.checkConsistency(true,false);
  //return 0;
  bool boundaryConditionFlag = true;

  double fixedLength;
  double dis; //extension for each deformation
  int beginDeformation;
  int endDeformation;
  double constraintEnergyRatio;
  double VredDiff;
  double tcDiff;
  double factor;
  ifstream fixinp("fixedDeformation.inp");
  fixinp >> fixedLength >> dis >> beginDeformation >> endDeformation >> constraintEnergyRatio >> VredDiff >> tcDiff>>factor;
  std::cout << "fixed Length = "  << fixedLength <<  std::endl
            << "extension of each deformation = " << dis <<std::endl
	    << "Begin deformation = " << beginDeformation <<std::endl
            << "End deformation = " <<endDeformation <<std::endl
	    << "constraint energy ratio = " << constraintEnergyRatio<<std::endl
	    << "reduced volume difference = "<< VredDiff <<std::endl
	    << "total curvature difference = "<< tcDiff <<std::endl
	    << "factor for AL = "<< factor <<std::endl;
   
  //ves11    
  vector<int> CN;
  double topZ=nodes[0]->getPoint(2);
  double endZ=nodes[nodes.size()-1]->getPoint(2);

  for (int a=0; a<nodes.size(); a++){
    double dt = std::abs(nodes[a]->getPoint(2)-topZ);
    double de = std::abs(nodes[a]->getPoint(2)-endZ);
    if (dt < fixedLength || de <fixedLength)
      CN.push_back(a);
  } 
  int NCN = CN.size();    


  cout<<"No. of fixed nodes "<<NCN<<endl; 
  //for (int i=0; i<NCN; i++){
  //cout<<CN[i]<<", ";
  //}
  
  double nu = vred;

  for(int deformationTime=beginDeformation; deformationTime < endDeformation+1; deformationTime++){
    //for (double nu=0.95; nu>0.93; nu-=0.02){
//set upper/lower bounds
  if(boundaryConditionFlag){

    blitz::Array<int,1> nbd(3*nodes.size());
    blitz::Array<double,1> lo(3*nodes.size());
    blitz::Array<double,1> hi(3*nodes.size());
    nbd = 0;
    lo = 0.0;
    hi = 0.0;


    for (int a=0; a<nodes.size(); a++){

      for (int i=0; i<NCN; i++){
	if (a==CN[i]){

	  nbd(3*a) = 2; nbd(3*a+1) = 2; nbd(3*a+2) = 2;

	    hi(3*a)=nodes[a]->getPoint(0);
	    lo(3*a)=nodes[a]->getPoint(0);  

	    hi(3*a+1)=nodes[a]->getPoint(1);
	    lo(3*a+1)=nodes[a]->getPoint(1);

	    if (nodes[a]->getPoint(2)< -0.0){
	      hi(3*a+2)=nodes[a]->getPoint(2)-dis;
	      lo(3*a+2)=nodes[a]->getPoint(2)-dis;
	    } 
	    else{
	      hi(3*a+2)=nodes[a]->getPoint(2)+dis;
	      lo(3*a+2)=nodes[a]->getPoint(2)+dis;
	    }
	}    
      }

      solver.setBounds(nbd, lo, hi);
    }
  }

  

  bd.reduceVolume(nu);
  /*tvmet::Vector<double,3> nf(0.0);

  for (double nf2=85.0; nf2 < 200.0; nf2+=5.0){
        
    nf(2)= nf2/2.0;
    topLoad->setForce(nf);
    nf(2)= -nf2/2.0;
    bottomLoad->setForce(nf);

    for(int i=0; i<nn; i++){
    nf(2)= nf2/nn/2.0;
    //neighborTopLoad(i)->setForce(nf);
    }

    for(int i=0; i<nn; i++){
    nf(2)= -nf2/nn/2.0;
    neighborBottomLoad(i)->setForce(nf);
    }*/   
  /*
    double energyDiff = 1.0;
    double CE;
    int iter=0, maxIter=10;


    bd.setCytoSpring(0.0, 0.0, 500.0);
 
    while (energyDiff>0.1 || energyDiff<0.0 && iter < maxIter ) {	
      double pCE = CE;
      solver.solve( &model );      
      CE = bd.constraintEnergy(); 
      iter++;
      std::cout<< "ConstraintEnergy = " << bd.constraintEnergy() << std::endl;
      bd.resetReference();
      energyDiff = pCE - CE;
    }

    char name[20]; sprintf(name,"f2-%i", deformationTime );
    model.print(name);
  */
  double diffVred=1.0;
  double diffTC=1.0;
  double fP=0.0;
  double fT=0.0;
  double fTC=0.0;

  bd.updateFixedForce(fP, fT, fTC);

  //bd.setCytoSpring(0.0, 0.0, 5.0e3);

  //model.checkConsistency(true,false);
  //return 0;
  //int aln=1;
    double energyRatio = 1.0;
    int iter=0, maxIter=100;

    while ( std::abs(diffVred) > VredDiff/*0.5e-6*/ || std::abs(diffTC) > tcDiff || energyRatio > constraintEnergyRatio/*5.0e-6*/ && iter < maxIter){


    //while ( energyRatio > 1.0e-5 && iter < maxIter ) {	
      solver.solve( &model );      
      energyRatio = bd.constraintEnergy() / bd.totalStrainEnergy(); 
      iter++;
      std::cout<< "ConstraintEnergy = " << bd.constraintEnergy() << std::endl;
      bd.resetReference();
      //}
    

      char name[20];sprintf(name,"ves-%i", deformationTime );
      model.print(name);
    //aln++;


//     char name[20]; sprintf(name,"nu%f-%f",nu, nf2);
    //char name[20]; sprintf(name,"f-%f", nf2);
    //model.print(name);
    
    //bd.resetReference();
    
    v = bd.volume();
    a = bd.area(); 
    tc = bd.totalCurvature();
    vred = 6.0*sqrt(M_PI)*v/std::pow(a,3.0/2.0);
    V = bd.constraintVolume();
    A = bd.constraintArea(); 
    TC = bd.constraintTotalCurvature();
    Vred = 6.0*sqrt(M_PI)*V/std::pow(A,3.0/2.0);

    
    diffVred=(vred-Vred)/Vred;
    diffTC=(tc-TC)/TC;
    fP += -penaltyVolume*(v-V)/V/V;
    fT += penaltyArea*(a-A)/A/A;
    fTC += penaltyTotalCurvature*(tc-TC)/TC/TC;
    bd.updateFixedForce(fP, fT, fTC);

    if (std::abs(diffVred) > VredDiff*factor /*1.5e-5*/){
      penaltyVolume*=2.0;
      penaltyArea*=2.0;
      penaltyTotalCurvature*=2.0;
      bd.updatePenaltyVolumeArea(penaltyVolume, penaltyArea, penaltyTotalCurvature);
    }
    std::cout  << "No. "<< deformationTime << std::endl
	       << "Volume = "<< bd.volume() << std::endl
	       << "Cons. Volume = "<< bd.constraintVolume() << std::endl
	       << "Area = "<< bd.area() << std::endl
	       << "Cons. Area = "<< bd.constraintArea() << std::endl
	       << "Reduced Volume = " << vred << std::endl
	       << "Ref. Reduced Volume = " << Vred << std::endl
	       << "Total Curvature= " << tc <<std::endl
	       << "Cons. Total Curvature= " << TC <<std::endl
	       << "Energy = " << solver.function() << std::endl
	       << "ConstraintEnergy = " << bd.constraintEnergy() << std::endl
	       << "strainEnergy = " << bd.totalStrainEnergy() << std::endl
               << "pressure = " << fP <<std::endl
	       << "tension = "  << fT <<std::endl;

    std::cout <<"top force: "<<bd.topForce()<<std::endl
	      <<"end force: "<<bd.endForce()<<std::endl;


    //model.checkRank(model.dof()-6,true);
    //}
  //model.checkRank(model.dof()-6,true);
 }

  }


 
  return 0;
}


void ioSetting(int argc, char* argv[], ifstream& ifs, string& ofn)
{

  // if no arguement or too many
  if (argc < 2 || argc > 3)
    {
      cout << "Usage: ProgramName InputFileName [OutputFileName]." << endl;
      exit(0);
    }

  string inFullName = argv[1];
  string pathName;
  //	basic_string <char>::size_type sztp = inFullName.find_last_of("/");
	
  // 	if ( (sztp) != string::npos )  // no position was found
  // 		pathName = inFullName.substr(0, sztp) + "/";
  // 	else
  // 		pathName = "";

  pathName = "";
	
  // only one arguement
  if ( argc == 2 ) ofn = pathName + "output.dat";


  // two arguements
  if ( argc == 3 ) {
    ofn = argv[2];
    ofn = pathName + ofn;
  }

  // create input stream
  ifs.open( inFullName.c_str(), ios::in);
  if (!ifs)
    {
      cout << "can not open input file: " << inFullName << endl;
      exit(0);
    }
	
}
