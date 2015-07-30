///////////////////////////////////////////////////////////////////////////////
///
///	\file    DataContainer.cpp
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

#include "DataContainer.h"

#include "Exception.h"

#include <cstdlib>

///////////////////////////////////////////////////////////////////////////////

DataContainer::DataContainer() :
	m_pAllocatedMemory(NULL)
{ }

///////////////////////////////////////////////////////////////////////////////

DataContainer::~DataContainer() {
	Deallocate();
}

///////////////////////////////////////////////////////////////////////////////

void DataContainer::Allocate() {

	// Check that memory has not been allocated
	if (m_pAllocatedMemory != NULL) {
		_EXCEPTIONT("Attempting to allocate already allocated DataContainer");
	}
	
	// Get the accumulated size of all DataChunks
	size_t sAccumulated = 0;

	for (size_t i = 0; i < m_vecDataChunks.size(); i++) {
		sAccumulated += m_vecDataChunks[i]->GetSize();
	}

	// Allocate memory as one contiguous chunk
	m_pAllocatedMemory =
		reinterpret_cast<unsigned char*>(malloc(sAccumulated));

	if (m_pAllocatedMemory == NULL) {
		_EXCEPTIONT("Out of memory");
	}

	// Assign memory to DataChunks
	sAccumulated = 0;

	for (size_t i = 0; i < m_vecDataChunks.size(); i++) {
		m_vecDataChunks[i]->AttachTo(
			reinterpret_cast<void *>(m_pAllocatedMemory + sAccumulated));

		sAccumulated += m_vecDataChunks[i]->GetSize();
	}
}

///////////////////////////////////////////////////////////////////////////////

void DataContainer::Deallocate() {

	for (size_t i = 0; i < m_vecDataChunks.size(); i++) {
		m_vecDataChunks[i]->Detach();
	}

	free(m_pAllocatedMemory);

	m_pAllocatedMemory = NULL;
}

///////////////////////////////////////////////////////////////////////////////

