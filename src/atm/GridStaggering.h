///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridStaggering.h
///	\author  Paul Ullrich
///	\version May 22, 2013
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

#ifndef _GRIDSTAGGERING_H_
#define _GRIDSTAGGERING_H_

#include "EquationSet.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Type of horizontal staggering.
///	</summary>
enum HorizontalStaggering {
	HorizontalStaggering_Elements,
	HorizontalStaggering_AlphaEdge,
	HorizontalStaggering_BetaEdge,
	HorizontalStaggering_Corners,
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Type of vertical staggering.
///	</summary>
enum VerticalStaggering {
	VerticalStaggering_Interfaces,
	VerticalStaggering_Levels,
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A descriptor for the staggering of prognostic variable solution nodes.
///	</summary>
class GridStaggering {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	GridStaggering() :
		m_fInitialized(false)
	{ }

public:
	///	<summary>
	///		Initialize the staggering information.
	///	</summary>
	void Initialize(
		const EquationSet & eqn
	) {
		m_fInitialized = true;

		// Resize arrays
		m_vecHorizStaggeringState.resize(eqn.GetComponents());
		m_vecVertStaggeringState.resize(eqn.GetComponents());

		// Initialize to default
		for (int c = 0; c < eqn.GetComponents(); c++) {
			m_vecHorizStaggeringState[c] = HorizontalStaggering_Elements;
			m_vecVertStaggeringState[c] = VerticalStaggering_Levels;
		}

		m_vecHorizStaggeringTracer = HorizontalStaggering_Elements;
		m_vecVertStaggeringTracer  = VerticalStaggering_Levels;
	}
/*
public:
	///	<summary>
	///		Get the horizontal staggering for the specified variable.
	///	</summary>
	HorizontalStaggering GetStaggering
*/
protected:
	///	<summary>
	///		Initialization flag.
	///	</summary>
	bool m_fInitialized;

	///	<summary>
	///		Vector of horizontal staggering of state variables.
	///	</summary>
	std::vector<HorizontalStaggering> m_vecHorizStaggeringState;

	///	<summary>
	///		Vector of vertical staggering of state variables.
	///	</summary>
	std::vector<VerticalStaggering> m_vecVertStaggeringState;

	///	<summary>
	///		Horizontal staggering of tracer variables.
	///	</summary>
	HorizontalStaggering m_vecHorizStaggeringTracer;

	///	<summary>
	///		Vertical staggering of tracer variables.
	///	</summary>
	VerticalStaggering m_vecVertStaggeringTracer;


};

///////////////////////////////////////////////////////////////////////////////

#endif

