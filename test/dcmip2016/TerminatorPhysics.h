///////////////////////////////////////////////////////////////////////////////
///
///	\file    TerminatorPhysics.cpp
///	\author  Paul Ullrich
///	\version June 1, 2016
///
///	<remarks>
///		Copyright 2000-2016 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _TERMINATORPHYSICS_H_
#define _TERMINATORPHYSICS_H_

#include "Defines.h"
#include "DataArray1D.h"
#include "DataArray2D.h"
#include "DataArray3D.h"
#include "DataArray4D.h"

#include "WorkflowProcess.h"

///////////////////////////////////////////////////////////////////////////////

class GridPatch;
class Time;

///	<summary>
///		KesslerPhysics type physics forcing.
///	</summary>
class TerminatorPhysics : public WorkflowProcess {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	TerminatorPhysics(
		Model & model,
		const Time & timeFrequency);
	
public:
  	///	<summary>
	///		Initializer.
	///	</summary>
	virtual void Initialize(
		const Time & timeStart
	);
	
	///	<summary>
	///		Apply KesslerPhysics physics.
	///	</summary>
	virtual void Perform(
		const Time & time
	);
};

///////////////////////////////////////////////////////////////////////////////

#endif
