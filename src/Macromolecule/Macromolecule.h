// Macromolecule.h: interface for the Macromolecule class.
//
//////////////////////////////////////////////////////////////////////

/*!  
  \file Atom.h

  \brief Macromolecule is a class that stores information about a macromolecule that is constructed from a number of atoms.
  */


#if !defined(AFX_MACROMOLECULE_H__DAEBAC11_9A78_4BC5_A549_509056537A20__INCLUDED_)
#define AFX_MACROMOLECULE_H__DAEBAC11_9A78_4BC5_A549_509056537A20__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>

#include <valarray>
#include <vector>
#include <set>
#include <map>
#include <tvmet/Vector.h>
#include <blitz/array.h>

typedef tvmet::Vector<double,3> Vector3D;
typedef tvmet::Vector<int,3> Vector3DInt;
using namespace std;

#include "../Atom/Atom.h"

class Macromolecule  
{
 public:

  //! Constructor, takes an ifstream as input, usually from .pdb file
  Macromolecule(ifstream &);  
  //! Default destructor
  virtual ~Macromolecule();
  //! Returns a bounding box that surrounds macromolecule, where the first vector, v[0] = MinPoint - 0.15*diag, and the second vector, v[1] = MaxPoint + 0.15*diag.
  const tvmet::Vector<Vector3D, 2> getBoundingBox();

  //! Returns container of all atoms in the macromolecule
  const vector<Atom >& getAtoms()  const { return _atoms; };
  //! Returns the number of atoms in the macromolecule
  unsigned int getNumberOfAtoms()  const { return( _atoms.size()  ); };

  //! Prints information about all atoms in the macromolecule
  void printAtoms();
  //! Based on coordinates of all atoms in macromolecule, a bounding box, and a grid size, the atomElementConnectivity vector is filled to represent which integer element number each atom belongs to
  void locateAtomsInGrid(std::vector<int>& atomElementConnectivity,
			 std::vector< std::vector<double> >& atomCoords,
			 const tvmet::Vector<tvmet::Vector<double,3>, 2>& box,
			 tvmet::Vector<double,3>& d,
			 const int gridSize);

 private:

  vector<Atom > _atoms;

  //! Calculates the scalar distance between the two input coordinates
  static double distance(const tvmet::Vector<double, 3>&, const tvmet::Vector<double, 3>&);
  //! Calculates the scalar product between the two input coordinates
  static double scalarProduct(const tvmet::Vector<double, 3>&, const tvmet::Vector<double, 3>&);
};

ostream& operator << (ostream&, const Atom& );

#endif // !defined(AFX_MACROMOLECULE_H__DAEBAC11_9A78_4BC5_A549_509056537A20__INCLUDED_)
