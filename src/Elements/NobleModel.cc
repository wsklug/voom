// -*- C++ -*-
//----------------------------------------------------------------------
//
//                     HoHai Van and William S. Klug
//                University of California Los Angeles
//                 (C) 2004-2005 All Rights Reserved
//
//----------------------------------------------------------------------
//
//----------------------------------------------------------------------

#include "IonicNoble.h"


int main {

//Read node coordinates for all the entire mesh

//Read in the connectivities for the elements

//Initialize and store the voltage and dv/dt at all the nodes

//Calculate the voltage at the quadrature points for each element and store

//Initialize and store m,n,h at the quadrature points

// Main loop over time

// Loop over each element

// Calculate the ionic current for each quadrature point in the element using m,n,h

/* Perform volume quadrature for each part of the differential equation.
   -(I_ion) + Laplacian_voltage + I_stim=dv/dt to create a the local element dv/dt vector.
*/

// Add the results of this to the global dv/dt vector

// End loop over each element

// Update voltage=old voltage + dt * dv/dt at the nodes

// Calculate the new voltage at the quadrature points

// Calculate alpha_i, beta_i for all the gating variables using new voltage at quad points

// Update m,n,h for example m=old m + dt * dm/dt at the quadrature points

// Output the new voltage for all the nodes in a file

// End loop over time



}
