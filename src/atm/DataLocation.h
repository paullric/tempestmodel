///////////////////////////////////////////////////////////////////////////////
///
///	\file    DataLocation.h
///	\author  Paul Ullrich
///	\version March 25, 2013
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

#ifndef _DATALOCATION_H_
#define _DATALOCATION_H_

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Discrete grid locations where data is stored.
///	</summary>
enum DataLocation {
	DataLocation_None = (-1),
	DataLocation_Node = (0),
	DataLocation_AEdge,
	DataLocation_BEdge,
	DataLocation_REdge,
	DataLocation_Count = DataLocation_REdge+1,
	DataLocation_Default = DataLocation_Node,
};

///////////////////////////////////////////////////////////////////////////////

#endif
