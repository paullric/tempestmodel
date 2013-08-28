///////////////////////////////////////////////////////////////////////////////
///
///	\file    WorkspaceMemory.cpp
///	\author  Paul Ullrich
///	\version January 18, 2012
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

#include "WorkspaceMemory.h"

#include "Exception.h"

///////////////////////////////////////////////////////////////////////////////
// MemoryAllotment
///////////////////////////////////////////////////////////////////////////////

MemoryAllotment::~MemoryAllotment() {
	m_aWorkspace.Free(m_iAllotmentId);
}

///////////////////////////////////////////////////////////////////////////////
// WorkspaceMemory
///////////////////////////////////////////////////////////////////////////////

void WorkspaceMemory::Initialize(
	int nMemorySize
) {
	// Check for validity of parameter
	if (nMemorySize <= 2) {
		_EXCEPTIONT("Workspace memory must be at least 3 bytes.");
	}

	// Already initialized
	if (m_nMemorySize != 0) {
		_EXCEPTIONT("Workspace memory already initialized.");
	}

	// Initialize
	m_nNextPointer = 1;
	m_nMemorySize = nMemorySize;

	m_pMemory = new char[nMemorySize];
	if (m_pMemory == NULL) {
		_EXCEPTIONT("Memory allocation error.");
	}

	// Set left block buffer byte
	m_pMemory[0] = 0xFF;
}

///////////////////////////////////////////////////////////////////////////////

void WorkspaceMemory::Get(
	int nSize,
	MemoryAllotment & aAllotment
) {

	// Verify sufficient memory is available
	if (m_nNextPointer + nSize + 1 > m_nMemorySize) {
		_EXCEPTIONT("Insufficient memory available.");
	}

	// Verify allotment is uninitialized
	if (aAllotment.m_iAllotmentId != (-1)) {
		_EXCEPTIONT("MemoryAllotment already allocated.");
	}

	// Allocate a new allotment
	aAllotment.m_iAllotmentId = static_cast<int>(m_vecAllotments.size());

	aAllotment.m_nMemorySize = nSize;

	aAllotment.m_pMemory = m_pMemory + m_nNextPointer;

	// Adjust the current pointer
	m_nNextPointer += nSize + 1;

	// Set the block buffer byte
	m_pMemory[m_nNextPointer] = 0xFF;

	// Store the allotment
	m_vecAllotments.push_back(&aAllotment);
}

///////////////////////////////////////////////////////////////////////////////

void WorkspaceMemory::Free(
	int iAllotmentId
) {
	// Verify the allotment id is legitimate
	if ((iAllotmentId < 0) || (iAllotmentId > m_vecAllotments.size())) {
		_EXCEPTIONT("Invalid allotment id.");
	}

	// Remove allotment
	m_vecAllotments.erase(m_vecAllotments.begin() + iAllotmentId);

	// If no allotments are left, reset
	if (m_vecAllotments.size() == 0) {
		m_nNextPointer = 1;
	}
}

///////////////////////////////////////////////////////////////////////////////

void WorkspaceMemory::CheckBounds() const {
	
	int nLocation = 0;

	if (m_pMemory[nLocation] != 0xFF) {
		_EXCEPTIONT("CheckBounds failure: Memory overflow detected.");
	}

	for (int i = 0; i < m_vecAllotments.size(); i++) {
		nLocation += m_vecAllotments[i]->m_nMemorySize + 1;

		if (m_pMemory[nLocation] != 0xFF) {
			_EXCEPTIONT("CheckBounds failure: Memory overflow detected.");
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

