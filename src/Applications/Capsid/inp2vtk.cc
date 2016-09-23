#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

int main(int argc, char * argv[]) {
  
  if(argc != 2 ) {
    std::cout << endl
	      << "Usage: inp2vtk modelName" << endl
	      << endl
	      << "modelName should be the prefix of the inp file you want to process (everything before the .inp extension)." 
	      << endl
	      << endl;

    return 0;
  }

  string fileName = argv[1];
    
  string inputFileName = fileName + ".inp";
  
  ifstream ifs( inputFileName.c_str() );
  if (!ifs) {
    cout << endl
	 << "Cannot open file " << inputFileName 
	 << endl
	 << endl;
    return (0);
  }
  
  cout << "Opened file " << inputFileName << endl;

  string outputFileName = fileName + ".vtk";
  ofstream ofs( outputFileName.c_str() );
 
  int N = 1;

  string line;
  int number;
  string name;

  std::vector<double> coords(0); // coordinates
  std::vector<int> connect(0); // coordinates

  int nNodes=0;
  int nElements=0;

  // look for the *Node card
  while ( getline( ifs, line ) ) {
    if ( line.substr(0,5) == "*Node" ) break;
    if ( line.substr(0,5) == "*NODE" ) break;

  }
  
  // read in nodes until *Element card is found
  while ( getline( ifs, line ) ) {



    if( line.substr(0,8) == "*Element" ) break;
    if( line.substr(0,8) == "*ELEMENT" ) break;


    // replace any commas with spaces
    for( string::size_type pos=line.find(","); 
	 pos != string::npos; 
	 pos=line.find(",") ) {
     
      line.replace(pos,1," ");

    }

    // parse line for node data
    istringstream sstr(line);
    int id;
    sstr >> id;
    for(int i=0; i<3; i++) {
      double x;
      sstr >> x;
      coords.push_back(x);
    }
    nNodes++;
  }

  // now we've found the *Element card, read in and store connectivities
  
  // first figure out what type of elements we have
  string::size_type pos = line.find("type=");
  if( pos == string::npos ) {
    pos = line.find("TYPE=");
  }

  if( pos == string::npos ) {

    cout << "Cannot determine element type.  Exiting." << endl;
    return 0; 
  }
  pos+=5;

  int vtk_cell_type=0;
  int nodes_per_element=0;
  if( line.substr(pos,2) == "S3" ) {
    // 3-node triangle
    vtk_cell_type = 5;
    nodes_per_element = 3;
  } else if( line.substr(pos,6) == "STRI65" ) {
    // 6-node triangle
    vtk_cell_type = 22;
    nodes_per_element = 6;
  } else {
    cout << "Cannot determine element type.  Quitting." << endl;
    return (0);
  }

  while ( getline( ifs, line ) ) {

    if( line[0] == '*' ) break;

    // replace any commas with spaces
    for( string::size_type pos=line.find(","); 
	 pos != string::npos; 
	 pos=line.find(",") ) {
     
      line.replace(pos,1," ");

    }
    
    istringstream sstr(line);
    int id;

    sstr >> id;

    for(int i=0; i<nodes_per_element; i++) {
      int a;
      sstr >> a;
      cout << a << " ";
      connect.push_back(a-1);
    }
    cout << endl;
    nElements++;
  }
  
  // print vtk file

  //    Node Section
  ofs << "# vtk DataFile Version 2.0\n"
      << "Test example" << std::endl
      << "ASCII" << std::endl
      << "DATASET UNSTRUCTURED_GRID" << std::endl
      << "POINTS  " << nNodes << "  double" << std::endl;

  for(int a=0; a<nNodes; a++) {
      ofs << std::setprecision(16) 
	  << coords[3*a+0] << "  "
	  << coords[3*a+1] << "  "
	  << coords[3*a+2] << endl;
  }

    //    Element Section
  ofs << "CELLS  " << nElements << " " << nElements*(1+nodes_per_element) 
      << std::endl;
  for(int f=0; f<nElements; f++) {
    ofs << nodes_per_element << "  ";
    for(int i=0; i<nodes_per_element; i++) {
      ofs << std::setw(10) << connect[nodes_per_element*f+i];
    }
    ofs << std::endl;
  }

 
  ofs << "CELL_TYPES  " << nElements << std::endl;
  for(int f=0; f<nElements; f++) {
    ofs << vtk_cell_type << endl;
  }
 
  ifs.close();
  ofs.close();

  cout << "Wrote file " << outputFileName << endl;

  return 0;

}
