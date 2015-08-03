///////////////////////////////////////////////////////////////////////////////
///
///	\file    DataContainer.h
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

#ifndef _DATACONTAINER_H_
#define _DATACONTAINER_H_

#include "DataChunk.h"

#include <vector>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class for storing a unified chunk of data.
///	</summary>
class DataContainer {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	DataContainer();

	///	<summary>
	///		Destructor.
	///	</summary>
	~DataContainer();

	///	<summary>
	///		Add a new data object to the DataContainer of specified size.
	///	</summary>
	size_t PushDataChunk(DataChunk * pChunk) {
		m_vecDataChunks.push_back(pChunk);

		return (m_vecDataChunks.size() - 1);
	}

	///	<summary>
	///		Allocate an array for all DataChunks.
	///	</summary>
	void Allocate();

	///	<summary>
	///		Deallocate all data and detach DataChunks.
	///	</summary>
	void Deallocate();

private:
	///	<summary>
	///		Pointer to allocated memory chunk.
	///	</summary>
	unsigned char * m_pAllocatedMemory;

	///	<summary>
	///		Vector of DataChunks stored in this DataContainer.
	///	</summary>
	std::vector<DataChunk *> m_vecDataChunks;
};

///////////////////////////////////////////////////////////////////////////////

#endif

