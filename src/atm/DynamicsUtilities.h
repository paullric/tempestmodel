///////////////////////////////////////////////////////////////////////////////
///
///	\file    DynamicsUtilities.h
///	\author  Paul Ullrich
///	\version March 25, 2015
///
///	<remarks>
///		Copyright 2000-2010 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _DYNAMICSUTILITIES_H_
#define _DYNAMICSUTILITIES_H_

#include "DataArray1D.h"
#include "DataArray2D.h"
#include "DataArray4D.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Finite-element method (FEM) based atmospheric horizontal dynamics.
///	</summary>
class DynamicsUtilities {

private:
	///	<summary>
	///		Private constructor.
	///	</summary>
	DynamicsUtilities()
	{ }

public:
	///	<summary>
	///		Apply boundary conditions to data.
	///	</summary>
	static void ApplyBoundaryConditions(
		const Grid & grid,
		int iDataIndex
	);

public:
	///	<summary>
	///		Compute total energy.
	///	</summary>
	static double ComputeTotalEnergy(
		const Grid & grid,
		int iDataIndex
	);

	///	<summary>
	///		Compute total potential enstrophy.
	///	</summary>
	static double ComputeTotalPotentialEnstrophy(
		const Grid & grid,
		int iDataIndex
	);
};

///////////////////////////////////////////////////////////////////////////////

#endif

