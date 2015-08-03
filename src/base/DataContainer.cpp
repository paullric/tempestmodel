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

#include <iostream>

///////////////////////////////////////////////////////////////////////////////

DataContainer::DataContainer() :
	m_fOwnsData(true),
	m_pAllocatedMemory(NULL)
{ }

///////////////////////////////////////////////////////////////////////////////

DataContainer::~DataContainer() {
	if ((m_fOwnsData) && (m_pAllocatedMemory != NULL)) {
		free(m_pAllocatedMemory);
	}
}

///////////////////////////////////////////////////////////////////////////////

size_t DataContainer::GetTotalByteSize() {
	
	// Get the accumulated size of all DataChunks
	size_t sAccumulated = 0;

	for (size_t i = 0; i < m_vecDataChunks.size(); i++) {
		sAccumulated += m_vecDataChunks[i]->GetByteSize();
	}

	return sAccumulated;
}

///////////////////////////////////////////////////////////////////////////////

void DataContainer::Allocate() {

	// Check that memory has not been allocated
	if (m_pAllocatedMemory != NULL) {
		_EXCEPTIONT("Attempting Allocate() on attached DataContainer");
	}

	// Get the accumulated size of all DataChunks
	size_t sAccumulated = GetTotalByteSize();

	// Allocate memory as one contiguous chunk
	m_pAllocatedMemory =
		reinterpret_cast<unsigned char*>(malloc(sAccumulated));

	if (m_pAllocatedMemory == NULL) {
		_EXCEPTIONT("Out of memory");
	}

	// Initialize allocated memory to zero
	memset(m_pAllocatedMemory, 0, sAccumulated);

	// Assign memory to DataChunks
	sAccumulated = 0;

	for (size_t i = 0; i < m_vecDataChunks.size(); i++) {
		m_vecDataChunks[i]->AttachToData(
			reinterpret_cast<void *>(m_pAllocatedMemory + sAccumulated));

		sAccumulated += m_vecDataChunks[i]->GetByteSize();
	}
}

///////////////////////////////////////////////////////////////////////////////

void DataContainer::Deallocate() {
	Detach();
}

///////////////////////////////////////////////////////////////////////////////

void DataContainer::AttachTo(
	unsigned char * pAllocatedMemory
) {
	if (m_pAllocatedMemory != NULL) {
		_EXCEPTIONT("Attempting AttachTo() on attached DataContainer");
	}

	m_fOwnsData = false;
	m_pAllocatedMemory = pAllocatedMemory;

	// Assign memory to DataChunks
	size_t sAccumulated = 0;

	for (size_t i = 0; i < m_vecDataChunks.size(); i++) {
		m_vecDataChunks[i]->AttachToData(
			reinterpret_cast<void *>(m_pAllocatedMemory + sAccumulated));

		sAccumulated += m_vecDataChunks[i]->GetByteSize();
	}
}

///////////////////////////////////////////////////////////////////////////////

void DataContainer::Detach() {

	for (size_t i = 0; i < m_vecDataChunks.size(); i++) {
		m_vecDataChunks[i]->Detach();
	}

	if ((m_fOwnsData) && (m_pAllocatedMemory != NULL)) {
		free(m_pAllocatedMemory);
	}

	m_fOwnsData = true;
	m_pAllocatedMemory = NULL;
}

///////////////////////////////////////////////////////////////////////////////

