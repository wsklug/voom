// -*- C++ -*-
//-------------------------------------------------------------------------
//
//                     HoHai Van and William S. Klug
//                University of California Los Angeles
//                     (C) 2006 All Rights Reserved
//
// 05/30/07: Revision #14
// Added cell_type definition to the element fiber file.
//
// 03/21/07: Revision #13.
// This version of app.cc is used for simulations involving fibrosis.
//
// 03/15/07: Revision #12.
// Revised the output to char (1 byte) instead of float (4 bytes) to 
// save filesize and time.
//
// 03/12/07: Revision #11.
// Revised the nodeconnect array to be compressed saving looping.
//
// 02/10/07: Revision #10.
// Added anisotropic fiber direction.
//
// 10/31/06: Revision #9.  
// Courtemache atrial cell kinetics.
//
// 10/17/06: Revision #8.
// Added hexadradral elements. Added Runge-Kutta 2 algorithm.  Moved 
// compute function to one function to speed things up.
//
// 10/09/06: Revision #7.  Added MPI functinality using metis as the
// mesh partitioner.  
//
// 09/25/06: Revision #6.  Added output to OpenDx.  Also implemented 3D
// formulation.
//
// 09/14/06: Revision #5.  Template style implemented by William Klug 
// and changed the reading of .node and .ele file to output by
// triangle meshing program instead of abaqus.
//
// 09/3/06: Revision #4: Moved all the gating variables to a table which
// does NOT speed up calculations so this ended up commented off.
// 
// 09/1/06: Revision #3: Added OpenMP for the for loop over the nodes and 
// elements.
//
// Revision #2: Now used the CardiacPotential class which is much faster
//
//---------------------------------------------------------------------------

#include "CP.fibrosis.h"
#include "CP.fibrosis.icc"
#include "Shape.h"
#include <string>
#include <vector>
#include <tvmet/Vector.h>
#include <tvmet/Matrix.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <fstream>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef WITH_MPI
#include <mpi.h>
#endif

#include "ShapeTri6.h"
#include "TriangleQuadrature.h"
#include "ShapeQ4.h"
#include "QuadQuadrature.h"
#include "ShapeTri3.h"
#include "TetQuadrature.h"
#include "ShapeTet4CP.h"
#include "HexQuadrature.h"
#include "ShapeHex8.h"


using namespace std;
using namespace voom;


/*! To switch element type, uncomment only one of these four  sets of lines.
***************************************************************************
Remember to change the quadrature rule when constructing the body of the 
mesh, i.e., TriangleQuadrature(1).
Also, the FEM.dx file needs to be modified to define the element type:
triangles, tetrahedra, cubes
***************************************************************************
*/

// typedef CardiacPotential<QuadQuadrature,ShapeQ4,2> ElementType;
// const int NODES_PER_ELEMENT=4;

// typedef CardiacPotential<TriangleQuadrature,ShapeTri3,2> ElementType;
// const int NODES_PER_ELEMENT=3;

// typedef CardiacPotential<TriangleQuadrature,ShapeTri6,2> ElementType;
// const int NODES_PER_ELEMENT=6;

// typedef CardiacPotential<TetQuadrature,ShapeTet4,3> ElementType;
// const int NODES_PER_ELEMENT=4;

 typedef CardiacPotential<HexQuadrature,ShapeHex8,3> ElementType;
 const int NODES_PER_ELEMENT=8;


/*! To switch dimensions, i.e., 2D vs 3D, uncomment only one of the two lines. 
*****************************************************************************
Remember to change the corresponding line in CP.fibrosis.h
*****************************************************************************
*/
//This is for 3D
typedef VoltageNode<3> VoltageNodeType;

//This is for 2D
//typedef VoltageNode<2> VoltageNodeType;

/*! n_datafiles is a global variable that is used in main and printdataDx.  It
    keeps track of the number of data files or frames that will be used in the movie.
*/
int n_datafiles=0;

void printmeshDx( ElementType::VoltageNodeContainer,
               vector< tvmet::Vector<int,NODES_PER_ELEMENT> >, 
               int); 
void printheaderDx(int ,int , int , int);

ofstream ofs("./data1/AP.dat");
ofstream ofsecg("./data1/ECG.dat");

/* Define constants here.  Make sure that it is initialized in the corresponding 
 * Element header file, i.e., CP.atria.h 
 */

float ElementType::gate_table[1701][25];
float ElementType::Gna=7.8;
float ElementType::GK1=3.0*0.09;
float ElementType::Gto=1.0*0.1652;
float ElementType::Gkur=1.0*1.0;
float ElementType::GCaL=1.5*0.1238;
float ElementType::GKr=1.0*0.0294; 
float ElementType::GKs=1.0*0.129; 
float ElementType::FVi=96.4867*13668.;        
float ElementType::VuVi=0.08117647;           
float ElementType::VrVi=0.00705882 ;          
float ElementType::VrVu=0.08695652 ;          
float ElementType::RT_F=26.71283192;          
float ElementType::ENa=RT_F*log(140.0/11.2);
float ElementType::EK=RT_F*log(5.4/139.0);
float ElementType::Trpn=0.07;
float ElementType::Cmdn=0.05;
float ElementType::Csqn=10.;

float ElementType::E_kf=-87.0;
float ElementType::g_kv=0.25;
float ElementType::g_k1f=.4822;
float ElementType::K_o=5.4;
float ElementType::Na_i=8.5547;
float ElementType::v_rev=-150.0;
float ElementType::inakf_bar=2.002;
float ElementType::K_mk=1.0;
float ElementType::K_mna=11.0;
float ElementType::g_bnaf=0.0095;
float ElementType::E_na_f=54.4;


int main(int argc, char* argv[])
{  

  if(argc != 2) {
    cout << "Usage: ./app modelName" << endl;
    return 0;
  }

  string modelName = argv[1];
  string nodeFileName = "../Mesh/" +modelName + ".node";
  string elementFileName = "../Mesh/" +modelName + ".ele";
  string elem_processors = "../Mesh/" +modelName + ".mesh.epart.32";
  string node_processors = "../Mesh/" +modelName + ".mesh.npart.32";
  string output = "../Mesh/" +modelName + ".nc";
  string fiberFileName = "../Mesh/" + modelName + ".fiber";
 
  int n_elem,n_nodes;
  string dummy_str;
  int dummy_int;
  double x0,y0,z0;
  int n_dim;
  int id;

  //  Read node coordinates for all the entire mesh
  ifstream inputnode(nodeFileName.c_str());
  inputnode>>n_nodes>>n_dim>>dummy_int>>dummy_int;

  VoltageNodeType::PositionVector x;
  VoltageNodeType::Point v;

  NodeBase::DofIndexMap idx(1);
  ElementType::VoltageNodeContainer nodes;
  nodes.clear();
  idx[0]=0;
  v=-82.0;

  for (int id=0;id<n_nodes;id++) {
    inputnode>>dummy_int>>x0>>y0;
    if (n_dim==3) {inputnode>>z0;}
    if (n_dim==2) {inputnode>>dummy_int;}

    x[0]=x0*0.025; x[1]=y0*0.025;
    if (n_dim==3) {x[2]=z0*0.025;}
    nodes.push_back(new VoltageNodeType(id,idx,x,v));

  }
//  cout<<"Number of nodes:  "<<nodes.size()<<endl;
  inputnode.close();

  //This reads in the element connectivities of the mesh
  tvmet::Vector<int,NODES_PER_ELEMENT> connect;
  vector< tvmet::Vector<int,NODES_PER_ELEMENT> > Connectivity;
  Connectivity.clear();
  ifstream inputele(elementFileName.c_str());
  inputele>>n_elem>>dummy_int>>dummy_int;
 
  for (int i=0;i<n_elem;i++) {
    connect =0;
    inputele >>dummy_int;
    for (int j=0;j<NODES_PER_ELEMENT;j++) {
        inputele>>connect[j];
        }
    connect -= 1;
    Connectivity.push_back(connect);
  }

//  cout<<"Number of elements:  " << Connectivity.size()<<endl;
  inputele.close();

  int _myRank=0;
  int _nProcessors;

#ifdef WITH_MPI
  MPI_Init (&argc,&argv); 
  MPI_Comm_size( MPI_COMM_WORLD, &_nProcessors );
  MPI_Comm_rank( MPI_COMM_WORLD, &_myRank );
#endif

  //read in the file for the elem processors
  ifstream elem_p(elem_processors.c_str());
  std::vector<int> _elem_proc;
  for (int i=0;i<n_elem;i++) {
     dummy_int=0;
#ifdef WITH_MPI
     elem_p>>dummy_int;
#endif
     _elem_proc.push_back(dummy_int);   
  }
  elem_p.close();

  //read in the file for the node processors
  ifstream node_p(node_processors.c_str());
  std::vector<int> _node_proc;
  for (int i=0;i<n_nodes;i++) {
     dummy_int=0;
#ifdef WITH_MPI
     node_p>>dummy_int;
#endif
     _node_proc.push_back(dummy_int);   
  }
  node_p.close();
  
  //read in the nc (nodeconnect) file used for mpi
#ifdef WITH_MPI
  ifstream nc(output.c_str());
  std::vector<int> nodeconnect;
  std::vector< std::vector<int> > nodeconn_container;
  int num_nc,c1,c2;
  nc>>num_nc; 
  for (int i=0;i<num_nc;i++) {
     nodeconnect.clear();
     nc>>c1;
     nodeconnect.push_back(c1);
     nc>>c1;
     for (int j=0;j<c1;j++) {
       nc>>c2;
       nodeconnect.push_back(c2);
     }
     nodeconn_container.push_back(nodeconnect);
  }
  nc.close();
#endif

  // Create body of elements
  vector<ElementType*>  ElementContainer;
  HexQuadrature q(1);
  int e_c=0;
  //Fiber direction file 
  ifstream iffiber(fiberFileName.c_str());
  for (vector< tvmet::Vector<int,NODES_PER_ELEMENT> >::const_iterator c=Connectivity.begin();c!=Connectivity.end();c++) {
    // build element node container
    ElementType::VoltageNodeContainer nds;
    nds.clear();
    for (tvmet::Vector<int,NODES_PER_ELEMENT>::const_iterator n=c->begin();n!=c->end();n++) {
      nds.push_back(nodes[*n]);
    }
    //  Load fiber and diffusion data from "fiber.dat"
    float fiber[3];
    float diffusion;
    int cell_type;
    iffiber>>cell_type;
    iffiber>>diffusion;
    iffiber>>fiber[0]>>fiber[1]>>fiber[2];
    //create the element
#ifdef WITH_MPI
    if (_myRank==_elem_proc[e_c]) ElementContainer.push_back(new ElementType(q,nds,diffusion,fiber,cell_type));
#endif
#ifndef WITH_MPI
    ElementContainer.push_back(new ElementType(q,nds,diffusion,fiber,cell_type));
#endif
    e_c++;
  }
  iffiber.close();
//  cout <<"Processor "<<_myRank<<" created body of "<<ElementContainer.size()<<" elements."<<endl;

  int es=0;
  int ee=ElementContainer.size();

  int ns=0;
  int ne=nodes.size();


  //Create array of elements to be stimulated; flag=true if stimulated
  std::vector<bool> stim_elem;
  for (int i=0;i<ElementContainer.size();i++) {
    ElementType::VoltageNodeContainer nds;
    nds=ElementContainer[i]->nodes();
    bool flag=true;
    for (int j=0;j<nds.size();j++) {
      double p=nds[j]->getPosition(0);
      if (p>=0.40) {flag=false;}
    }
    stim_elem.push_back(flag);
  }

  //Create array of elements to be stimulated; flag=true if stimulated
  std::vector<bool> stim2_elem;
  for (int i=0;i<ElementContainer.size();i++) {
    ElementType::VoltageNodeContainer nds;
    nds=ElementContainer[i]->nodes();
    bool flag=true;
    for (int j=0;j<nds.size();j++) {
      double p1=nds[j]->getPosition(0);
      double p2=nds[j]->getPosition(1);
      double p3=nds[j]->getPosition(2);

//    Sinus node
      if (p1<=1.75) {flag=false;}
      if (p1>=2.25) {flag=false;}
      if (p2<=3.00) {flag=false;}
      if (p2>=3.75) {flag=false;}
      if (p3<=4.50) {flag=false;}
      if (p3>=5.00) {flag=false;}


    }
    stim2_elem.push_back(flag);
  }

  //Create array of elements to be stimulated; flag=true if stimulated
  //This is a 3D stimulus protocol
  std::vector<bool> stim3_elem;
  for (int i=0;i<ElementContainer.size();i++) {
    ElementType::VoltageNodeContainer nds;
    nds=ElementContainer[i]->nodes();
    bool flag=true;
    for (int j=0;j<nds.size();j++) {
      double p1=nds[j]->getPosition(0);
      double p2=nds[j]->getPosition(1);
      double p3=nds[j]->getPosition(2);

//      if (p1<=(161*0.025)) {flag=false;}
//      if (p1>=(179*0.025)) {flag=false;}
//      if (p2<=(181*0.025)) {flag=false;}
//      if (p2>=(188*0.025)) {flag=false;}
//      if (p3<=(4*0.025)) {flag=false;}
//      if (p3>=(31*0.025)) {flag=false;}

//    LSPV
      if (p1<=5.10) {flag=false;}
      if (p1>=5.35) {flag=false;}
      if (p2<=6.65) {flag=false;}
      if (p2>=6.95) {flag=false;}
      if (p3<=2.90) {flag=false;}
      if (p3>=3.15) {flag=false;}



//    RSPV
/*
      if (p1<=70*0.025) {flag=false;}
      if (p1>=80*0.025) {flag=false;}
      if (p2<=230*0.025) {flag=false;}
      if (p2>=240*0.025) {flag=false;}
      if (p3<=120*0.025) {flag=false;}
      if (p3>=130*0.025) {flag=false;}

      



//    Sinus node
      if (p1<=1.75) {flag=false;}
      if (p1>=2.25) {flag=false;}
      if (p2<=3.00) {flag=false;}
      if (p2>=3.75) {flag=false;}
      if (p3<=4.50) {flag=false;}
      if (p3>=5.00) {flag=false;}
*/
    }
    stim3_elem.push_back(flag);
  }


#ifdef WITH_MPI
// Initialize the damping
    for (int i=0;i<nodeconn_container.size();i++) {
       int nc=nodeconn_container[i][0];
       int c3;
       if (_myRank==_node_proc[nc]) {
          for (int j=1;j<nodeconn_container[i].size();j++) {
              c3=nodeconn_container[i][j];
              if (c3!=_myRank) {
                 float other_damping=0.0;
                 int source=c3;
                 MPI_Status status;
                 MPI_Recv(&other_damping,1,MPI_FLOAT,source,0,MPI_COMM_WORLD,&status);
                 nodes[nc]->addDamping(0,other_damping);
                 
              }          
          }
       }
       else {
          for (int j=1;j<nodeconn_container[i].size();j++) {
               c3=nodeconn_container[i][j];
               if (c3==_myRank) {
                  float damping=nodes[nc]->getDamping(0);
                  int dest=_node_proc[nc];         
                  MPI_Send(&damping,1,MPI_FLOAT,dest,0,MPI_COMM_WORLD);   
               }
           }
        }

    }
#endif

// make gate_table
  for (int i=0;i<1701;i++) {   
     double volt_quad=-90+i*0.1;
          

     double a_m;
     if (fabs(volt_quad+47.13)<=0.001) a_m=3.2;
     else a_m=0.32*(volt_quad+47.13)/(1.0-exp(-0.1*(volt_quad+47.13)));
     double b_m=0.08*exp(-volt_quad/11.0);
     double m_tau=(a_m+b_m);
     double m_inf=a_m/m_tau;
     ElementType::gate_table[i][0]=m_tau;
     ElementType::gate_table[i][1]=m_inf;
          
     double a_h,b_h,a_j,b_j;
     if (volt_quad>=-40.0) {
        a_h=0.0;
        b_h=1.0/(0.13*(1.0+exp((volt_quad+10.66)/-11.1)));
        a_j=0.0;
        b_j=0.3*exp(-0.0000002535*volt_quad)/(1.0+exp(-0.1*(volt_quad+32.0)));
        }
     else {
        a_h=0.135*exp((80.0+volt_quad)/-6.8);
        b_h=3.56*exp(0.079*volt_quad)+310000.0*exp(0.35*volt_quad);
        a_j=(-127140.0*exp(0.2444*volt_quad)-0.00003474*exp(-0.04391*volt_quad))*(volt_quad+37.78)/(1+exp(0.311*(volt_quad+79.23)));
        b_j=0.1212*exp(-0.01052*volt_quad)/(1.0+exp(-0.1378*(volt_quad+40.14)));

        }  

     double h_tau=(a_h+b_h);
     double h_inf=a_h/h_tau;
     double j_tau=(a_j+b_j);
     double j_inf=a_j/j_tau;
     ElementType::gate_table[i][2]=h_tau;
     ElementType::gate_table[i][3]=h_inf;
     ElementType::gate_table[i][4]=j_tau;
     ElementType::gate_table[i][5]=j_inf;

     double a_oa=0.65/(exp(-(volt_quad+10.0)/8.5)+exp(-(volt_quad-30.0)/59.0));
     double b_oa=0.65/(2.5+exp((volt_quad+82.)/17.0));
     double oa_inf=1.0/(1.0+exp(-(volt_quad+20.47)/17.54));
     double oa_tau=(a_oa+b_oa)*3.0;
     ElementType::gate_table[i][6]=oa_tau;
     ElementType::gate_table[i][7]=oa_inf;
     

     double a_oi=1./(18.53+exp((volt_quad+113.7)/10.95));
     double b_oi=1./(35.56+exp(-(volt_quad+1.26)/7.44));
     double oi_inf=1./(1.+exp((volt_quad+43.1)/5.3));
     double oi_tau=(a_oi+b_oi)*3.0;
     ElementType::gate_table[i][8]=oi_tau;
     ElementType::gate_table[i][9]=oi_inf;

     double a_ua=0.65/(exp(-(volt_quad+10.)/8.5)+exp(-(volt_quad-30.)/59.0));
     double b_ua=0.65/(2.5+exp((volt_quad+82.)/17.0));
     double ua_inf=1./(1+exp(-(volt_quad+30.3)/9.6));
     double ua_tau=(a_ua+b_ua)*3.0;
     ElementType::gate_table[i][10]=ua_tau;
     ElementType::gate_table[i][11]=ua_inf;
         
     double a_ui=1./(21.+exp(-(volt_quad-185.)/28.));
     double b_ui=exp((volt_quad-158.)/16.); 
     double ui_inf=1./(1.+exp((volt_quad-99.45)/27.48));
     double ui_tau=(a_ui+b_ui)*3.0;
     ElementType::gate_table[i][12]=ui_tau;
     ElementType::gate_table[i][13]=ui_inf;

     double a_xr;
     if (fabs(volt_quad+14.1)<=0.0001)
        a_xr=0.0015;
     else
        a_xr=3.e-4*(volt_quad+14.1)/(1.-exp(-(volt_quad+14.1)*.2));
     double b_xr;
     if (fabs(volt_quad-3.3328)<=0.0001)
        b_xr=7.3898e-5*5.1237;
     else
        b_xr=7.3898e-5*(volt_quad-3.3328)/(exp((volt_quad-3.3328)/5.1237)-1.);
     double xr_inf=1./(1.+exp(-(volt_quad+14.1)/6.5));
     double xr_tau=(a_xr+b_xr);
     ElementType::gate_table[i][14]=xr_tau;
     ElementType::gate_table[i][15]=xr_inf;

     double a_xs,b_xs;         
     if (fabs(volt_quad-19.9)<=1.e-4) {
        a_xs=6.8e-4;
        b_xs=3.15e-4;
        }
     else {
        a_xs=4.e-5*(volt_quad-19.9)/(1.-exp(-(volt_quad-19.9)/17.));
        b_xs=3.5e-5*(volt_quad-19.9)/(exp((volt_quad-19.9)/9.)-1.);
        } 
     double xs_inf=1./sqrt(1.+exp(-(volt_quad-19.9)/12.7));
     double xs_tau=2.*(a_xs+b_xs);
     ElementType::gate_table[i][16]=xs_tau;
     ElementType::gate_table[i][17]=xs_inf;

     double d_tau;
     if (fabs(volt_quad+10.0)<=.0001) d_tau=0.4368;
     else d_tau=(0.035*(volt_quad+10)*(1.0+exp(-(volt_quad+10.0)/6.24)))/(1.0-exp(-(volt_quad+10.0)/6.24));
     double d_inf=1.0/(1.0+exp(-(volt_quad+10.0)/8.0));
     ElementType::gate_table[i][18]=d_tau;
     ElementType::gate_table[i][19]=d_inf;


     double f_tau=(0.0197*exp(-(0.0337*0.0337*(volt_quad+10.0)*(volt_quad+10.0)))+0.02)/9;
     double f_inf=1/(1+exp((volt_quad+28.0)/6.9));
     ElementType::gate_table[i][20]=f_tau;
     ElementType::gate_table[i][21]=f_inf;

     double w_tau;
     if (fabs(volt_quad-7.9)<=1.e-4)
        w_tau=6.5/6.;
     else {
        w_tau=(volt_quad-7.9)*(1.+0.3*exp(-(volt_quad-7.9)*.2));
        w_tau/=(6.*(1.-exp(-(volt_quad-7.9)*.2)));}
     double w_inf=1.-1./(1.+exp(-(volt_quad-40.)/17.));
     ElementType::gate_table[i][22]=w_tau;
     ElementType::gate_table[i][23]=w_inf;


  }

  double curr_time=0.0;
  // The time step of the system
  float dt=0.05;
  // The stimulus strength;
  double istim=0.0;
  // These variables correspond to the when to start and end the pacing
  int    t_0,t_f,t_02,t_f2,t_03,t_f3;
  // The pulse 
  int    t_pulse=62;
  int    t_pulse_ACh=100;
  int    t_init=100;
  int    t_period=20000;
  int    t_init2=107900;
  int    t_period2=3200;
  int    t_period_ACh=10000;
  int    t_init_ACh=15000;
//  std::cout <<"Processor "<<_myRank<< " Starting time loop now" << std::endl;
  bool startflag=true;
  
  // Creates the mesh.data file to be read by OpenDx
  if (_myRank==0) {printmeshDx(nodes,Connectivity,n_dim);
    cout<<"CL=160ms"<<endl;

  }



  //Pacer 1
  t_0=t_init;
  t_f=t_init+t_pulse;

  //Pacer 2
  t_02=t_init2;
  t_f2=t_init2+t_pulse;
  
  //ACh stimulation
  t_03=t_init_ACh;
  t_f3=t_init_ACh+t_pulse_ACh;
  
  
  for  (int output_cnt=0;output_cnt<140000;output_cnt++) {

    //Increase the next stimulus time step
    if (output_cnt>t_f){
      t_0+=t_period;
      t_f+=t_period;
    }

    if (output_cnt>t_f2){
      t_02+=t_period2;
      t_f2+=t_period2;
    }

    if (output_cnt>t_f3){
      t_03+=t_period_ACh;
      t_f3+=t_period_ACh;
    }
    // Turn off pacing after a set time point
    if (curr_time>=5100.0) {
      t_0=0;
      t_f=0;
    }

    if (curr_time>=62000.0) {
      t_02=0;
      t_f2=0;
    }

    if (curr_time>=4000.0) {
      t_03=0;
      t_f3=0;
    }

//   Used for RK2 looping
//    nodes_old.clear();
//    for (int i=0;i<n_nodes;i++) {
//       if (_myRank==_node_proc[i])
//          nodes_old.push_back(nodes[i]->getPoint(0));	  
//    }

//    double _dt=dt*.5;
//    bool rk2_flag=false;
//    for (int rk2=1;rk2<2;rk2++) {

//    if (rk2==1) {rk2_flag=true;_dt=dt;}

   
    // Set the force contribution to each node to zero
#ifdef _OPENMP
#pragma omp parallel for			\
  schedule(static) default(shared) 
     for (int i=0;i<n_nodes;i++) {
	   nodes[i]->setForce(0.0);
     }
#endif

#ifdef WITH_MPI
    for (int i=0;i<nodeconn_container.size();i++) {
       int nc=nodeconn_container[i][0];
       int c3;
       if (_myRank==_node_proc[nc]) {
          for (int j=1;j<nodeconn_container[i].size();j++) {
              c3=nodeconn_container[i][j];
              if (c3!=_myRank) {
                 float volt=nodes[nc]->getPoint(0);
                 int dest=c3;         
                 MPI_Send(&volt,1,MPI_FLOAT,dest,0,MPI_COMM_WORLD);   
                 
              }          
          }
       }
       else {
          for (int j=1;j<nodeconn_container[i].size();j++) {
               c3=nodeconn_container[i][j];
               if (c3==_myRank) {
                 float volt=0.0;
                 int source=_node_proc[nc];
                 MPI_Status status;
                 MPI_Recv(&volt,1,MPI_FLOAT,source,0,MPI_COMM_WORLD,&status);
                 nodes[nc]->setPoint(0,volt);
                 nodes[nc]->setForce(0.0);
               }
           }
        }

    }

#endif


    // Loop over each element in the body
#ifdef _OPENMP
#pragma omp parallel for			\
  schedule(static) default(shared) 
#endif
    for (int elem_cnt=0;elem_cnt<ElementContainer.size();elem_cnt++) {
      istim=0.0;
      if (((output_cnt >=t_0) && (output_cnt<=t_f) && (stim2_elem[elem_cnt])))
         {istim=20.0;}
      if (((output_cnt >=t_02) && (output_cnt<=t_f2) && (stim3_elem[elem_cnt])))
         {istim=20.0;}
      float ACh=0.0;
//      if ((output_cnt >=t_03) && (output_cnt<=t_f3))
//         {ACh=0.03;}
      
      ElementContainer[elem_cnt]->compute_ion(dt,istim,true,ACh);
    }

   

#ifdef WITH_MPI
    for (int i=0;i<nodeconn_container.size();i++) {
       int nc=nodeconn_container[i][0];
       int c3;
       if (_myRank==_node_proc[nc]) {
          for (int j=1;j<nodeconn_container[i].size();j++) {
              c3=nodeconn_container[i][j];
              if (c3!=_myRank) {
                 float other_force=0.0;
                 int source=c3;
                 MPI_Status status;
                 MPI_Recv(&other_force,1,MPI_FLOAT,source,0,MPI_COMM_WORLD,&status);
                 nodes[nc]->addForce(0,other_force);
                 
              }          
          }
       }
       else {
          for (int j=1;j<nodeconn_container[i].size();j++) {
               c3=nodeconn_container[i][j];
               if (c3==_myRank) {
                  float force=nodes[nc]->getForce(0);
                  int dest=_node_proc[nc];         
                  MPI_Send(&force,1,MPI_FLOAT,dest,0,MPI_COMM_WORLD);   
               }
           }
        }

    }
#endif

    // Calculate the new voltages by the RK2 method in time and add to each node

#ifdef _OPENMP
#pragma omp parallel for			\
  schedule(static) default(shared) 
#endif
    for (int i=0;i<n_nodes;i++) {
      if (_myRank==_node_proc[i]) {
         float f=nodes[i]->getForce(0);
         float d=nodes[i]->getDamping(0);
         nodes[i]->addPoint(0,dt*f/d);
         nodes[i]->setForce(0,0.0);
      }
    }

  
//    } //end rk2 loop

    curr_time+=dt;

    //Output voltage of all the nodes at certain time step

#ifdef WITH_MPI   
    if ((output_cnt %200==0) && (curr_time>5490.0)) {
       char fname[200];
       n_datafiles++;
       double step=n_datafiles;
       double rank=_myRank;
       sprintf(fname, "/u/home/htvan/hvan/voom-cardiac-export/src/Applications/Fibrosis/data1/volt.%05.0f.data.%02.0f", step,rank);
       ofstream dataofs;
       dataofs.open(fname);
       
       for (int i=0;i<n_nodes;i++) {
          if (_node_proc[i]==_myRank) {
             double volt=nodes[i]->getPoint(0);
             int volt_int=(volt+85)*(250/70);
             if (volt_int>250) volt_int=250;
             if (volt_int<0) volt_int=0;
             char volt_out=volt_int;
	           dataofs.write( (char *) &volt_out, sizeof( volt_out ) );
             if ((i==1)) {ofs<<curr_time<<" "<<volt<<endl;}
	        }
       }
       dataofs.close();
    }
#endif

#ifdef WITH_MPI   
    if ((output_cnt %100==0) && (curr_time>1990.0)) {
       
       // Compute ECG for all the elements and output to file
       tvmet::Vector<float,3> ref(0.0);
       ref(0)=-3.0;
       ref(1)=-3.0;
       float ECG=0.0;
       float total_ECG1=0.0;
       for (int elem_cnt=0;elem_cnt<ElementContainer.size();elem_cnt++) {
          ECG+=ElementContainer[elem_cnt]->compute_ECG(ref);
       }
       MPI_Reduce(&ECG,&total_ECG1,1,MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
       ref(0)=3.0;
       ref(1)=-3.0;
       ECG=0.0;
       float total_ECG2=0.0;
       for (int elem_cnt=0;elem_cnt<ElementContainer.size();elem_cnt++) {
          ECG+=ElementContainer[elem_cnt]->compute_ECG(ref);
       }
       MPI_Reduce(&ECG,&total_ECG2,1,MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
       if (_myRank==0) {ofsecg<<curr_time<<" "<<(total_ECG1-total_ECG2)<<endl;}
    }
#endif

#ifdef _OPENMP
    if (output_cnt%200==0) {
       char fname[200];
       n_datafiles++;
       double step=n_datafiles;
       sprintf(fname, "/u/home/htvan/hvan/voom-cardiac-export/src/Applications/Fibrosis/data1/volt.%05.0f.data", step);
       ofstream dataofs;
       dataofs.open(fname);       
       for (int i=0;i<n_nodes;i++) {
         double volt=nodes[i]->getPoint(0);
         int volt_int=(volt+81)*(126/81);
         if (volt_int>126) volt_int=126;
         if (volt_int<0) volt_int=0;
         char volt_out=volt_int;
//         if ((volt_out>=60.0) || (volt_out<=-90.0)) {
//           cout<<curr_time<<" Node ="<<i<<" Vm="<<volt_out<<endl;}
	       dataofs.write( (char *) &volt_out, sizeof( volt_out ) );
         if (i==0) {cout<<curr_time<<" "<<i<<" "<<volt<<endl;}
//                              ofs<<curr_time<<" "<<volt<<endl;}
	   

       }
       dataofs.close();
    }
#endif

// Stop when wave hits edge of tissue to calculate wavespeed
//    double q1=nodes[0]->getPoint(0);
//    double start_time;
//    if ((q1>=0.0)&&(startflag)) {start_time=curr_time/*;cout<<curr_time<<endl*/;startflag=false;}
//    double q2=nodes[1]->getPoint(0);
//    if ((q2>=0.0)) {cout<<curr_time-start_time<<endl;/*curr_time=20000.1;*/}


//    output_cnt+=1;


  } //end loop over time

  if (_myRank==0) { printheaderDx(nodes.size(),n_elem,n_dim,NODES_PER_ELEMENT);}

#ifdef WITH_MPI
  MPI_Finalize();
#endif
  ofs.close();
  ofsecg.close();

  return 0;
} //end main{}


void printVTK( ElementType::VoltageNodeContainer nodes,
               vector< tvmet::Vector<int,NODES_PER_ELEMENT> > con,
	       double time ) 
{

  char fname[100];
  sprintf(fname, "C.%07.2f.vtk", time);
  ofstream ofs(fname);
     
  ofs << "# vtk DataFile Version 2.0\n"
      << "Test example" << std::endl
      << "ASCII" << std::endl
      << "DATASET POLYDATA" << std::endl
      << "POINTS  " << nodes.size() << "  double" << std::endl;
    
  //
  // output nodal postions
  for(int i=0;i<nodes.size();i++) {
    double x=nodes[i]->getPosition(0);
    double y=nodes[i]->getPosition(1);
    ofs << std::setprecision(16) 
    	<< x << "  "
	<< y << "  "
	<< 0.0 << std::endl;
  }



  /////////////////////////////////////////////////////////////////////
  //
  //    Element Section
  //
  ofs << "POLYGONS  " << con.size() << "  "
      << (1+NODES_PER_ELEMENT)*con.size() << std::endl;
  if( NODES_PER_ELEMENT == 6 ) {
    for(int e=0; e<con.size(); e++) {      
      ofs << NODES_PER_ELEMENT << "  ";
      ofs << std::setw(10) << con[e](0);
      ofs << std::setw(10) << con[e](3);
      ofs << std::setw(10) << con[e](1);
      ofs << std::setw(10) << con[e](4);
      ofs << std::setw(10) << con[e](2);
      ofs << std::setw(10) << con[e](5);
      ofs << std::endl;
    }
  } else {
    for(int e=0; e<con.size(); e++) {      
      ofs << NODES_PER_ELEMENT << "  ";
      for(int a=0; a<NODES_PER_ELEMENT; a++)
	ofs << std::setw(10) << con[e](a);
      ofs << std::endl;
    }
  } 

  ofs << endl << "POINT_DATA " << nodes.size() << endl
      << "SCALARS  voltage  double  "<< endl
      << "LOOKUP_TABLE default" << endl;
  for(int n=0; n<nodes.size(); n++) {
    double p=nodes[n]->getPoint(0);
    ofs << std::setprecision(16) 
	<< p << std::endl;
  }

  ofs.close();
  return;
}

namespace voom {
  template<>
  void CardiacPotential<TriangleQuadrature,ShapeTri6,2>::computeLumpedCapacitance() {
    double c_m=1.0;
    //Compute area
    double area=0.0;
    for(QuadPointIterator p=_quadPoints.begin();p!=_quadPoints.end(); p++){
      area += p->weight;
    }
  
    for(int a=0; a<3; a++) {
      _vNodes[a]->addDamping(0,area/12.0);
    }
    for(int a=3; a<6; a++) {
      _vNodes[a]->addDamping(0,area/4.0);
    }
  
    return;
  }
}
void printdataDx( ElementType::VoltageNodeContainer nodes) 
{
  
  char fname[100];
  n_datafiles+=1;
  double step=n_datafiles;
  sprintf(fname, "volt.%05.0f.data", step);
  ofstream ofs(fname);

  for(int n=0; n<nodes.size(); n++) {
    float p=nodes[n]->getPoint(0);
    ofs.write( (char *) &p, sizeof( p ) );
  }
  ofs.close();
  return;
}

/*!  This function prints out the mesh information for OpenDx.  The node
x,y,z coordinates are printed first, then the element connectivities.

*/
void printmeshDx( ElementType::VoltageNodeContainer nodes,
               vector< tvmet::Vector<int,NODES_PER_ELEMENT> > con, 
               int d) 
{

  ofstream ofs("./data1/mesh.data");
  // output nodal postions
  for(int i=0;i<nodes.size();i++) {
    float x=nodes[i]->getPosition(0);
    float y=nodes[i]->getPosition(1);
    ofs.write( (char *) &x, sizeof( x ) );
    ofs.write( (char *) &y, sizeof( y ) );
    

    if (d==3) {
       float z=nodes[i]->getPosition(2);
       ofs.write( (char *) &z, sizeof( z ) );
    }
  }

  //    Element Section

  if( NODES_PER_ELEMENT == 8 ) {
    for (int e=0;e<con.size();e++) {
    int h=con[e](0);
    ofs.write((char*)&h,sizeof(h));
    h=con[e](4);
    ofs.write((char*)&h,sizeof(h));
    h=con[e](3);
    ofs.write((char*)&h,sizeof(h));
    h=con[e](7);
    ofs.write((char*)&h,sizeof(h));
    h=con[e](1);
    ofs.write((char*)&h,sizeof(h));
    h=con[e](5);
    ofs.write((char*)&h,sizeof(h));
    h=con[e](2);
    ofs.write((char*)&h,sizeof(h));
    h=con[e](6);
    ofs.write((char*)&h,sizeof(h));
    }

  } else {
    for(int e=0; e<con.size(); e++) {      
      for(int a=0; a<NODES_PER_ELEMENT; a++) {
        int h=con[e](a);
        ofs.write( (char *) &h, sizeof( h ) );

       }
    }
  } 

  ofs.close();
  return;
}

/*! Makes the header file for OpenDx.  The header files contain all sequence of 
    frames that will be used in the simulation movie
*/

void printheaderDx(int s,int r, int d, int n) {
ofstream Dxofs("./data1/FEM.dx");

Dxofs << "object 1 class array type float rank 1 "
      << "shape " << d
      << " items " << s<<" ieee"<<endl
      << "data file \"./mesh.data\",0" << endl
      <<"attribute \"dep\" string \"positions\""<<endl<<endl;

Dxofs << "object 2 class array type int rank 1 "
      << "shape "<<n 
      << " items " << r<<" ieee"<<endl
      << "data file \"./mesh.data\","<<s*4*d << endl
      <<"attribute \"element type\" string \"cubes\""<<endl
      <<"attribute \"ref\" string \"positions\""<<endl<<endl;

for (int i=0;i<n_datafiles;i++) {
  char fname[100];
  double step=i+1;
  sprintf(fname, "volt.%05.0f.data", step);
  Dxofs<<"object "<<i+1+2<<" class array type unsigned byte rank 0 items "
       <<s<<" ieee"<<endl;
  Dxofs<<"data file \"./"<<fname<<"\",0"<<endl
       <<"attribute \"dep\" string \"positions\""<<endl<<endl;

  }


for (int i=0;i<n_datafiles;i++) {
  char fname[100];
  double step=i+1;
  sprintf(fname, "volt.%05.0f.data", step);
  Dxofs<<"object "<<n_datafiles+i+1+2<<" class field"<<endl
  <<"component \"positions\" value 1"<<endl
  <<"component \"connections\" value 2"<<endl
  <<"component \"data\" value "<<i+1+2<<endl<<endl; 
  }

Dxofs<<"object \"series\" class series"<<endl;
for (int i=0;i<n_datafiles;i++) {
  Dxofs<<"member "<<i<<" value "<<n_datafiles+i+1+2<<" position "<<i+1<<" "<<endl;
  }
Dxofs<<endl<<"end"<<endl;
Dxofs.close();


}


