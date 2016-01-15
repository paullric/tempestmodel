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
#include <cstring>

#include <iostream>

///////////////////////////////////////////////////////////////////////////////

static const size_t int_aligned_size = sizeof(size_t);
static const size_t float_aligned_size = sizeof(size_t);
static const size_t double_aligned_size = sizeof(double);
static const size_t size_t_aligned_size = sizeof(size_t);
static const size_t bool_aligned_size = sizeof(size_t);

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

size_t DataContainer::GetTotalByteSize() const {
	
	// Get the accumulated size of all DataChunks
	size_t sAccumulated = 0;

	for (size_t i = 0; i < m_vecDataChunks.size(); i++) {
		if (m_vecType[i] == PRIMTYPE_DATACHUNK) {
			DataChunk * pDataChunk =
				reinterpret_cast<DataChunk*>(m_vecDataChunks[i]);

			// Verify byte alignment
			size_t sByteSize = pDataChunk->GetByteSize();
			if (sByteSize % sizeof(size_t) != 0) {
				_EXCEPTIONT("Misaligned array detected in DataContainer");
			}
			sAccumulated += sByteSize;

		} else if (m_vecType[i] == PRIMTYPE_INT) {
			sAccumulated += int_aligned_size;
		} else if (m_vecType[i] == PRIMTYPE_FLOAT) {
			sAccumulated += float_aligned_size;
		} else if (m_vecType[i] == PRIMTYPE_DOUBLE) {
			sAccumulated += double_aligned_size;
		} else if (m_vecType[i] == PRIMTYPE_SIZET) {
			sAccumulated += size_t_aligned_size;
		} else if (m_vecType[i] == PRIMTYPE_BOOL) {
			sAccumulated += bool_aligned_size;
		} else {
			_EXCEPTIONT("Invalid PRIMTYPE");
		}
	}

	return sAccumulated;
}

///////////////////////////////////////////////////////////////////////////////

void DataContainer::Allocate() {

	// Check that memory has not been allocated
	if (m_pAllocatedMemory != NULL) {
		_EXCEPTIONT("Attempting Allocate() on attached DataContainer");
	}

	// Allocate memory as one contiguous chunk
	size_t sTotalByteSize = GetTotalByteSize();

	m_pAllocatedMemory =
		reinterpret_cast<unsigned char*>(malloc(sTotalByteSize));

	if (m_pAllocatedMemory == NULL) {
		_EXCEPTIONT("Out of memory");
	}

	// Initialize allocated memory to zero
	memset(m_pAllocatedMemory, 0, sTotalByteSize);

	// Assign memory to DataChunks
	unsigned char * pAccumulated = m_pAllocatedMemory;

	for (size_t i = 0; i < m_vecDataChunks.size(); i++) {

		// DataChunk type
		if (m_vecType[i] == PRIMTYPE_DATACHUNK) {
			DataChunk * pDataChunk =
				reinterpret_cast<DataChunk*>(m_vecDataChunks[i]);

			pDataChunk->AttachToData(
				reinterpret_cast<void *>(pAccumulated));

#pragma message "Alignment may be an issue here on 32-bit systems"
			pAccumulated += pDataChunk->GetByteSize();

		// int type
		} else if (m_vecType[i] == PRIMTYPE_INT) {
			pAccumulated += int_aligned_size;

		// float type
		} else if (m_vecType[i] == PRIMTYPE_FLOAT) {
			pAccumulated += float_aligned_size;

		// double type
		} else if (m_vecType[i] == PRIMTYPE_DOUBLE) {
			pAccumulated += double_aligned_size;

		// size_t type
		} else if (m_vecType[i] == PRIMTYPE_SIZET) {
			pAccumulated += size_t_aligned_size;

		// bool type
		} else if (m_vecType[i] == PRIMTYPE_BOOL) {
			pAccumulated += bool_aligned_size;

		// Invalid type
		} else {
			_EXCEPTIONT("Invalid PRIMTYPE");
		}
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
	unsigned char * pAccumulated = m_pAllocatedMemory;

	for (size_t i = 0; i < m_vecDataChunks.size(); i++) {

		// DataChunk type
		if (m_vecType[i] == PRIMTYPE_DATACHUNK) {
			DataChunk * pDataChunk =
				reinterpret_cast<DataChunk *>(m_vecDataChunks[i]);
			pDataChunk->AttachToData(
				reinterpret_cast<void *>(pAccumulated));

			pAccumulated += pDataChunk->GetByteSize();

		// int type
		} else if (m_vecType[i] == PRIMTYPE_INT) {
			int * pInt = reinterpret_cast<int *>(m_vecDataChunks[i]);
			(*pInt) = *(reinterpret_cast<int *>(pAccumulated));

			pAccumulated += int_aligned_size;

		// float type
		} else if (m_vecType[i] == PRIMTYPE_FLOAT) {
			float * pFloat = reinterpret_cast<float *>(m_vecDataChunks[i]);
			(*pFloat) = *(reinterpret_cast<float *>(pAccumulated));

			pAccumulated += float_aligned_size;

		// double type
		} else if (m_vecType[i] == PRIMTYPE_DOUBLE) {
			double * pDouble = reinterpret_cast<double *>(m_vecDataChunks[i]);
			(*pDouble) = *(reinterpret_cast<size_t *>(pAccumulated));

			pAccumulated += double_aligned_size;

		// size_t type
		} else if (m_vecType[i] == PRIMTYPE_SIZET) {
			size_t * pSizeT = reinterpret_cast<size_t *>(m_vecDataChunks[i]);
			(*pSizeT) = *(reinterpret_cast<size_t *>(pAccumulated));

			pAccumulated += size_t_aligned_size;

		// bool type
		} else if (m_vecType[i] == PRIMTYPE_BOOL) {
			bool * pBool = reinterpret_cast<bool *>(m_vecDataChunks[i]);
			(*pBool) = *(reinterpret_cast<bool *>(pAccumulated));

			pAccumulated += bool_aligned_size;

		// Invalid type
		} else {
			_EXCEPTIONT("Invalid PRIMTYPE");
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void DataContainer::Detach() {

	for (size_t i = 0; i < m_vecDataChunks.size(); i++) {
		if (m_vecType[i] == PRIMTYPE_DATACHUNK) {
			DataChunk * pDataChunk =
				reinterpret_cast<DataChunk *>(m_vecDataChunks[i]);

			pDataChunk->Detach();
		}
	}

	if ((m_fOwnsData) && (m_pAllocatedMemory != NULL)) {
		free(m_pAllocatedMemory);
	}

	m_fOwnsData = true;
	m_pAllocatedMemory = NULL;
}

///////////////////////////////////////////////////////////////////////////////

