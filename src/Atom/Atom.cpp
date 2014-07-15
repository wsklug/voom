// Atom.cpp: implementation of the Atom class.
//
//////////////////////////////////////////////////////////////////////

#include "Atom.h"
#include <iostream>
#include <cstdlib>
#include <cctype>
#include <cmath>
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
Atom::Atom()
{
	_position = 0.0, 0.0, 0.0;

}
Atom::Atom(const string& line)
{
	_position = 0.0, 0.0, 0.0;
	
	unsigned int lineLength = line.size();
	if ( lineLength < 54 ) {
		// throw an exception
	}		
	else {
		_number = atoi( line.substr(7,5).c_str() );
		_name = line.substr(12,3);
		for ( string::size_type i=_name.find(" "); i!=string::npos; i=_name.find(" ") ) {
			_name.erase(i,1);
		}
		for ( unsigned int i=0; i<_name.length(); i++ ) {
			if ( isdigit(_name[i]) ) _name.erase(i,1);
		}
		_altLoc = line[15];
		_residueName = line.substr(17,3);
		_chainId = line[21];
		_residueNumber = atoi( line.substr(22,4).c_str() );
		_insertionCode = line[26];
		_position[0] = atof( line.substr(30,8).c_str() );
		_position[1] = atof( line.substr(38,8).c_str() );
		_position[2] = atof( line.substr(46,8).c_str() );
	}

	if ( lineLength >= 60 )
		_occupancy = atof( line.substr(54,6).c_str() );

	if ( lineLength >= 66 )
		_tempFactor = atof( line.substr(60,6).c_str() );

	if ( lineLength >= 76 )
		_segId = line.substr(72,4);

	if ( lineLength >= 78 )
		_element = line.substr(76,2);

	if ( lineLength >= 80 )
		_charge = line.substr(78,2);

}

Atom::Atom(const Atom &input)
{
	if ( input._position.size() != 3 ) {
	// throw an exception
	}

	_position = 0.0, 0.0, 0.0;

	_number        = input._number;
	_name          = input._name;
	_altLoc        = input._altLoc;
	_residueName   = input._residueName;
	_chainId       = input._chainId;
	_residueNumber = input._residueNumber;
	_insertionCode = input._insertionCode;
	_position      = input._position;
	_occupancy 	   = input._occupancy;
	_tempFactor    = input._tempFactor;
	_segId         = input._segId;
	_element       = input._element;
	_charge        = input._charge;
	
}

Atom::~Atom()
{

}

double Atom::getPosition(unsigned int i) const
{
	if ( i < _position.size() ) {
		return _position[i];
	} else {
		// throw exception
	}
	return ( 0 );	
}
