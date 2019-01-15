///////////////////////////////////////////////////////////////////////////////
///
///	\file    HeldSuarezPhysics.h
///	\author  Paul Ullrich
///	\version June 27, 2014
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

#ifndef _HELDSUAREZPHYSICS_H_
#define _HELDSUAREZPHYSICS_H_

#include "WorkflowProcess.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Held-Suarez type physics forcing.
///	</summary>
class HeldSuarezPhysics : public WorkflowProcess {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	HeldSuarezPhysics(
		Model & model,
		const Time & timeFrequency
	);

public:
	///	<summary>
	///		Apply Held-Suarez physics.
	///	</summary>
	virtual void Perform(
		const Time & time
	);
};

///////////////////////////////////////////////////////////////////////////////

#endif

