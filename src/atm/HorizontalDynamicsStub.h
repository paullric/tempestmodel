///////////////////////////////////////////////////////////////////////////////
///
///	\file    HorizontalDynamicsStub.h
///	\author  Paul Ullrich
///	\version August 31, 2015
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

#ifndef _HORIZONTALDYNAMICSSTUB_H_
#define _HORIZONTALDYNAMICSSTUB_H_

#include "HorizontalDynamics.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Stub for horizontal dynamics.
///	</summary>
class HorizontalDynamicsStub : public HorizontalDynamics {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	HorizontalDynamicsStub(
		Model & model
	) :
		HorizontalDynamics(model)
	{ }
};

///////////////////////////////////////////////////////////////////////////////

#endif

