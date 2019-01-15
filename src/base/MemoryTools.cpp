///////////////////////////////////////////////////////////////////////////////
///
///	\file    MemoryTools.cpp
///	\author  Paul Ullrich
///	\version April 23, 2014
///
///	<summary>
///		This header file provides access to various mathematical functions
///		that are not normally exposed in the standard C++ libraries.
///	</summary>
///	<remarks>
///		Copyright 2000-2010 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#include "MemoryTools.h"
#include "Announce.h"

///////////////////////////////////////////////////////////////////////////////

void PrintMemoryLine(const char * szString) {
	rusage ruse;

	getrusage(RUSAGE_SELF, &ruse);

	if (szString == NULL) {
		Announce("MEMORY RES %lu", ruse.ru_maxrss);
	} else {
		Announce("%s : RES %lu : DATA %lu : STACK %lu",
			szString, ruse.ru_maxrss);
	}
}

///////////////////////////////////////////////////////////////////////////////

