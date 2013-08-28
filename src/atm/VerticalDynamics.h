///////////////////////////////////////////////////////////////////////////////
///
///	\file    VerticalDynamics.h
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

#ifndef _VERTICALDYNAMICS_H_
#define _VERTICALDYNAMICS_H_

#include "WorkflowProcess.h"

///////////////////////////////////////////////////////////////////////////////

class Model;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Interface for vertical dynamics in the atmospheric model.
///	</summary>
class VerticalDynamics : public WorkflowProcess {

protected:
	///	<summary>
	///		Constructor.  Do not allow initialization of standalone
	///		VerticalDynamics.
	///	</summary>
	VerticalDynamics(
		Model & model
	) :
		WorkflowProcess(model)
	{ }

public:
	///	<summary>
	///		Number of auxiliary data objects needed.
	///	</summary>
	virtual int GetAuxDataCount() const {
		return 2;
	}
};

///////////////////////////////////////////////////////////////////////////////

#endif

