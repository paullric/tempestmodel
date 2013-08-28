///////////////////////////////////////////////////////////////////////////////
///
///	\file    Direction.h
///	\author  Paul Ullrich
///	\version March 24, 2013
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

#ifndef _CHECKSUMTYPE_H_
#define _CHECKSUMTYPE_H_

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Type of checksum to perform.
///	</summary>
enum ChecksumType {
	ChecksumType_Sum,
	ChecksumType_L1,
	ChecksumType_L2,
	ChecksumType_Linf
};

///////////////////////////////////////////////////////////////////////////////

#endif

