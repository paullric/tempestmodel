///////////////////////////////////////////////////////////////////////////////
///
///	\file    OutputManagerChecksum.h
///	\author  Paul Ullrich
///	\version May 15, 2013
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

#ifndef _OUTPUTMANAGERCHECKSUM_H_
#define _OUTPUTMANAGERCHECKSUM_H_

#include "OutputManager.h"

class Time;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		An OutputManager which provides checksums to standard output.  This
///		type of OutputManager is used to verify conservation properties.
///	</summary>
class OutputManagerChecksum : public OutputManager {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	OutputManagerChecksum(
		Grid & grid,
		const Time & timeOutputFrequency
	);

	///	<summary>
	///		Get the name of the OutputManager.
	///	</summary>
	virtual const char * GetName() const {
		return "Checksum ";
	}

protected:
	///	<summary>
	///		Perform an output.
	///	</summary>
	void Output(
		const Time & time
	);
};

///////////////////////////////////////////////////////////////////////////////

#endif

