// Atom.h: interface for the Atom class.
//
//////////////////////////////////////////////////////////////////////

/*!  
  \file Atom.h

  \brief Atom is a class that stores information about a single atom that is constructed (usually) from a single line of a .pdb file.  This information includes its number, name, Cartesian coordinates, etc.
  */

#if !defined(AFX_ATOM_H__99F556F5_4288_4E5E_896F_13D5DA9099E6__INCLUDED_)
#define AFX_ATOM_H__99F556F5_4288_4E5E_896F_13D5DA9099E6__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <valarray>
#include <string>
#include <tvmet/Vector.h>

using namespace std;

class Atom  
{
 public:
	
  typedef tvmet::Vector<double,3> Vector3D;
  
  //! Default constructor
  Atom();
  //! Default destructor
  virtual ~Atom();
  //! Construct with .pdb string input
  Atom(const string &);
  //! Copy constructor
  Atom(const Atom &);

  //! Returns atom number 
  inline unsigned int getNumber() const {return _number;}
  //! Returns atom name
  inline string getName() const {return _name;}
  //! Returns residue number
  inline unsigned int getResidueNumber() const {return _residueNumber;}
  //  inline const valarray <double> &getPosition() const {return _position;}
  //! Returns position vector of atom
  inline const Vector3D  &getPosition() const {return _position;}
  //! Returns desired component of position vector (0, 1, or 2)
  double getPosition(unsigned int i) const;
	
 private:
  // Columns     Description
  unsigned int      _number;       //  7-11  Atom serial number.
  string            _name;         // 13-16 Atom name.
  char              _altLoc;       //   17  Alternate location indicator.
  string            _residueName;  // 18-20 Residue name.
  char              _chainId;      //   22  Chain identifier.
  unsigned int      _residueNumber;// 23-26 Residue sequence number.
  char              _insertionCode;//   27  Code for insertion of residues.
  Vector3D          _position;     // 31-54 Orthogonal coordinates in angstroms.	
  //       (31-38 x,39-46 y, 47-54 z)
  double            _occupancy;	 // 55-60 Occupancy.
  double            _tempFactor; // 61-66 Temperature factor.
  string            _segId;      // 73-76 Segment identifier.
  string            _element;    // 77-78 Element symbol.
  string            _charge;     // 79-80 Charge on the atom.

};

#endif // !defined(AFX_ATOM_H__99F556F5_4288_4E5E_896F_13D5DA9099E6__INCLUDED_)
