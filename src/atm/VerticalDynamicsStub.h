///////////////////////////////////////////////////////////////////////////////
///
///	\file    VerticalDynamicsStub.h
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

#ifndef _VERTICALDYNAMICSSTUB_H_
#define _VERTICALDYNAMICSSTUB_H_

#include "VerticalDynamics.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Stub for atmospheric vertical dynamics.
///	</summary>
class VerticalDynamicsStub : public VerticalDynamics {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	VerticalDynamicsStub(
		Model & model
	) :
		VerticalDynamics(model)
	{ }
};

///////////////////////////////////////////////////////////////////////////////

#endif

