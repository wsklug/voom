// Macromolecule.cpp: implementation of the Macromolecule class.
//
//////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstring>
#include <ctime>
#include <iostream>
#include <algorithm>
#include <vector>

#include "Macromolecule.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Macromolecule::Macromolecule(ifstream & ifs)
{
  // Search through ifs for ATOM records and construct atoms in vector _atoms
  string line;
  set< string > skippedRecords;
  while ( getline( ifs, line ) ) {
    string recordName = line.substr(0,6);
    if ( recordName == "ATOM  " ) {
      _atoms.push_back( Atom( line ) );
    } else {
      skippedRecords.insert( recordName );
    }
  }
  if ( !skippedRecords.empty() ) {
    cout << "Records with the following names were skipped "
	 << "because they are not yet supported." << endl;
    for ( set< string >::iterator i = skippedRecords.begin();
	  i != skippedRecords.end();
	  i++ ) {
      cout << "\t" << *i << endl;
    }
  }
 //  vector<Atom>::iterator iA=_atoms.begin();
//   unsigned int a=0;
//   for( ; iA!=_atoms.end(); iA++, a++) 
//     {
//       _atomNumbers[&(*iA)] = a;
//     }
}


Macromolecule::~Macromolecule()
{

}

//////////////////////////////////////////////////////////////////////
// Methods
//////////////////////////////////////////////////////////////////////

const tvmet::Vector<Vector3D, 2>  Macromolecule::getBoundingBox()
{
  tvmet::Vector<Vector3D, 2> v;

  if ( !_atoms.empty() ) {
	  Vector3D Min/*(3)*/;
	  Vector3D Max/*(3)*/;
	  Vector3D MinNew;
	  Vector3D MaxNew;
	
    // compute min and max coordinates
    vector<Atom>::iterator i = _atoms.begin();
    Min =i -> getPosition();
	Max = i -> getPosition();
    i++;
    for ( ; i != _atoms.end(); i++ ) {
      for ( unsigned int j=0; j < 3; j++ ) {
	Min[j] = min( Min[j], i->getPosition( j ) );
	Max[j] = max( Max[j], i->getPosition( j ) );	
	MinNew[j] = Min[j];
	MaxNew[j] = Max[j];
      }
    }

    Vector3D diag;
    diag = MaxNew - MinNew;
    
    v[0] = MinNew - 0.15*diag;
    v[1] = MaxNew + 0.15*diag;
    
  }
  return ( v );
}
	

double Macromolecule::distance( const tvmet::Vector<double, 3> & x, 
				const tvmet::Vector<double, 3> & y )
{
  if (x.size() != y.size()) 
    printf("Error detected by Macromolecule::distance: unmatched dimensions");

  tvmet::Vector<double, 3> d;
  d = x - y;
  d = d*d;
  //for (unsigned int i=0; i<x.size(); i++) a += pow(x[i] - y[i], 2);
  return( sqrt( sum(d) ) );
	
	
}

double Macromolecule::scalarProduct( const tvmet::Vector<double, 3>& x, 
				     const tvmet::Vector<double, 3>& y )
{
  if (x.size() != y.size()) 
    printf("Error detected by Macromolecule::scalarProduct: unmatched dimensions");

  //for (unsigned int i=0; i<x.size(); i++) a += pow(x[i] - y[i], 2);
  return( sum(x*y) );
}

void Macromolecule::printAtoms()
{
  vector<Atom>::iterator i = _atoms.begin();
  for ( ; i != _atoms.end(); i++ ) {
    cout << *i << endl;
  }

}


ostream& operator<< (ostream& s, const Atom& a)	
{
  return s<<'{'<<a.getPosition(0)
	  <<','<<a.getPosition(1)
	  <<','<<a.getPosition(2)<<'}';

}

void Macromolecule::locateAtomsInGrid(std::vector<int>& atomElementConnectivity,
				      std::vector< std::vector<double> >& atomCoords,
				      const tvmet::Vector<tvmet::Vector<double,3>, 2>& box,
				      tvmet::Vector<double,3>& d,
				      const int gridSize) {

  vector<Atom>::iterator a = _atoms.begin();
  for ( ; a != _atoms.end(); a++ ) {
    Vector3D pos = a->getPosition();
    std::vector<double> tempCoord(3);
    tempCoord[0] = pos[0];
    tempCoord[1] = pos[1];
    tempCoord[2] = pos[2];
    atomCoords.push_back(tempCoord);
    
    int xval = 0;
    int yval = 0;
    int zval = 0;

    xval = fabs(pos[0]-box[0][0])/d[0];
    yval = fabs(pos[1]-box[0][1])/d[1];
    zval = fabs(pos[2]-box[0][2])/d[2];

    int elemNum = 0;
    elemNum = (xval+1) + yval*gridSize + zval*gridSize*gridSize;

    atomElementConnectivity.push_back(elemNum);    
  }

  return;
}
// end of file
