///////////////////////////////////////////////////////////////////////////////
///
///	\file    HorizontalDynamics.h
///	\author  Paul Ullrich
///	\version June 18, 2013
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

#ifndef _HORIZONTALDYNAMICS_H_
#define _HORIZONTALDYNAMICS_H_

#include "WorkflowProcess.h"

///////////////////////////////////////////////////////////////////////////////

class Model;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Interface for atmospheric horizontal dynamics.
///	</summary>
class HorizontalDynamics : public WorkflowProcess {

protected:
	///	<summary>
	///		Constructor.  Do not allow initialization of standalone
	///		HorizontalDynamics.
	///	</summary>
	HorizontalDynamics(
		Model & model
	) :
		WorkflowProcess(model)
	{ }

public:
	///	<summary>
	///		Virtual destructor.
	///	</summary>
	virtual ~HorizontalDynamics() { }

public:
	///	<summary>
	///		Get the number of halo elements needed by the model.
	///	</summary>
	virtual int GetHaloElements() const = 0;
};

///////////////////////////////////////////////////////////////////////////////

#endif

