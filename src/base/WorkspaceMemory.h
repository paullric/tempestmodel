///////////////////////////////////////////////////////////////////////////////
///
///	\file    WorkspaceMemory.h
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

#include "Exception.h"

#include <vector>

///////////////////////////////////////////////////////////////////////////////

class WorkspaceMemory;

///////////////////////////////////////////////////////////////////////////////

class MemoryAllotment {

friend class WorkspaceMemory;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	MemoryAllotment(
		WorkspaceMemory & aWorkspace
	) :
		m_aWorkspace(aWorkspace),
		m_iAllotmentId(-1),
		m_nMemorySize(0),
		m_pMemory(NULL)
	{ }

private:
	///	<summary>
	///		Copy constructor (not allowed).
	///	</summary>
	MemoryAllotment(const MemoryAllotment & aAllotment) :
		m_aWorkspace(aAllotment.m_aWorkspace)
	{
		_EXCEPTION();
	}

public:
	///	<summary>
	///		Destructor.
	///	</summary>
	~MemoryAllotment();

private:
	///	<summary>
	///		Assignment operator.
	///	</summary>
	MemoryAllotment & operator=(const MemoryAllotment & aAllotment) {
		_EXCEPTION();
	}

private:
	///	<summary>
	///		Parent workspace memory.
	///	</summary>
	WorkspaceMemory & m_aWorkspace;

	///	<summary>
	///		Identifier of this memory allotment.
	///	</summary>
	int m_iAllotmentId;

	///	<summary>
	///		Allotment size.
	///	</summary>
	int m_nMemorySize;

	///	<summary>
	///		Pointer to memory allotment.
	///	</summary>
	char * m_pMemory;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class used for managing a chunk of scratch workspace data.
///	</summary>
class WorkspaceMemory {

friend class MemoryAllotment;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	WorkspaceMemory() :
		m_nMemorySize(0),
		m_nNextPointer(0),
		m_pMemory(NULL)
	{ }

public:
	///	<summary>
	///		Destructor.
	///	</summary>
	~WorkspaceMemory();

public:
	///	<summary>
	///		Initialize this workspace memory.
	///	</summary>
	void Initialize(int nMemorySize);

public:
	///	<summary>
	///		Request an allotment of memory of a specified size.
	///	</summary>
	void Get(int nSize, MemoryAllotment & aAllotment);

private:
	///	<summary>
	///		Free an allotment of memory.
	///	</summary>
	void Free(int iAllotmentId);

public:
	///	<summary>
	///		Check memory bounds to ensure no overwrites have occurred.
	///	</summary>
	void CheckBounds() const;

private:
	///	<summary>
	///		Next available pointer.
	///	</summary>
	int m_nNextPointer;

	///	<summary>
	///		Size of the currently allocated chunk of memory.
	///	</summary>
	int m_nMemorySize;

	///	<summary>
	///		Currently allocated chunk of memory.
	///	</summary>
	char * m_pMemory;

	///	<summary>
	///		Set of all memory allotments.
	///	</summary>
	std::vector<MemoryAllotment *> m_vecAllotments;
};

