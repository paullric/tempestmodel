///////////////////////////////////////////////////////////////////////////////
///
///	\file    DataChunk.h
///	\author  Paul Ullrich
///	\version June 26, 2015
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

#ifndef _DATACHUNK_H_
#define _DATACHUNK_H_

#include <string>
#include <vector>

///////////////////////////////////////////////////////////////////////////////

class DataChunk {

public:
	///	<summary>
	///		Virtual destructor.
	///	</summary>
	virtual ~DataChunk()
	{ }

	///	<summary>
	///		Get the size of this DataChunk, in bytes.
	///	</summary>
	virtual size_t GetByteSize() const = 0;

	///	<summary>
	///		Determine if this DataChunk is attached to a data array.
	///	</summary>
	virtual bool IsAttached() const = 0;

	///	<summary>
	///		Attach this DataChunk to an array of pre-allocated data.
	///	</summary>
	virtual void AttachToData(void *) = 0;

	///	<summary>
	///		Detach data from this DataChunk.
	///	</summary>
	virtual void Detach() = 0;
};

///////////////////////////////////////////////////////////////////////////////

#endif

