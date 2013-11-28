///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridLonLat.h
///	\author  Paul Ullrich
///	\version February 25, 2013
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

#ifndef _GRIDLONLAT_H_
#define _GRIDLONLAT_H_

#include "Grid.h"
#include "GridPatch.h"

class Time;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Grid patch on the longitude-latitude grid.
///	</summary>
class GridPatchLonLat : public GridPatch {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	GridPatchLonLat(
		Grid & grid,
		int ixPatch,
		const PatchBox & box
	) :
		GridPatch(grid, ixPatch, box)
	{ }

public:
	///	<summary>
	///		Initialize data arrays.
	///	</summary>
	virtual void InitializeData();

	///	<summary>
	///		Initialize state and tracer data from a test case.
	///	</summary>
	void EvaluateTestCase(
		const TestCase & test,
		const Time & time,
		int iDataInstance = 0
	);
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		The longitude-latitude implementation of Grid.
///	</summary>
class GridLonLat : public Grid {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	GridLonLat(
		Model & model,
		int nLongitudes,
		int nLatitudes,
		int nRefinementRatio,
		int nRElements
	);
};

///////////////////////////////////////////////////////////////////////////////

#endif

