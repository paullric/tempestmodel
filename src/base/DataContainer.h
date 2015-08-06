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
	///		Primtive types stored in the DataContainer.
	///	</summary>
	enum PRIMTYPE {
		PRIMTYPE_DATACHUNK,
		PRIMTYPE_INT,
		PRIMTYPE_FLOAT,
		PRIMTYPE_DOUBLE,
		PRIMTYPE_SIZET,
		PRIMTYPE_BOOL
	};

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
		m_vecType.push_back(PRIMTYPE_DATACHUNK);

		return (m_vecDataChunks.size() - 1);
	}

	///	<summary>
	///		Add a new integer to the DataContainer.
	///	</summary>
	size_t PushInteger(int * pInt) {
		m_vecDataChunks.push_back(pInt);
		m_vecType.push_back(PRIMTYPE_INT);

		return (m_vecDataChunks.size() - 1);
	}

	///	<summary>
	///		Add a new float to the DataContainer.
	///	</summary>
	size_t PushFloat(float * pFloat) {
		m_vecDataChunks.push_back(pFloat);
		m_vecType.push_back(PRIMTYPE_FLOAT);

		return (m_vecDataChunks.size() - 1);
	}

	///	<summary>
	///		Add a new double to the DataContainer.
	///	</summary>
	size_t PushDouble(double * pDouble) {
		m_vecDataChunks.push_back(pDouble);
		m_vecType.push_back(PRIMTYPE_DOUBLE);

		return (m_vecDataChunks.size() - 1);
	}

	///	<summary>
	///		Add a new size_t to the DataContainer.
	///	</summary>
	size_t PushSizeT(size_t * pSizeT) {
		m_vecDataChunks.push_back(pSizeT);
		m_vecType.push_back(PRIMTYPE_SIZET);

		return (m_vecDataChunks.size() - 1);
	}

	///	<summary>
	///		Add a new size_t to the DataContainer.
	///	</summary>
	size_t PushBool(bool * pBool) {
		m_vecDataChunks.push_back(pBool);
		m_vecType.push_back(PRIMTYPE_BOOL);

		return (m_vecDataChunks.size() - 1);
	}

	///	<summary>
	///		Get the total size of the DataContainer (in bytes).
	///	</summary>
	size_t GetTotalByteSize();

	///	<summary>
	///		Allocate an array for all DataChunks.
	///	</summary>
	void Allocate();

	///	<summary>
	///		Deallocate all data and detach DataChunks.
	///	</summary>
	void Deallocate();

	///	<summary>
	///		Attach to the specified pointer.
	///	</summary>
	void AttachTo(unsigned char * pAllocatedMemory);

	///	<summary>
	///		Detach from the specified pointer.
	///	</summary>
	void Detach();

private:
	///	<summary>
	///		Flag indicating that this DataContainer owns is memory space.
	///	</summary>
	bool m_fOwnsData;

	///	<summary>
	///		Pointer to allocated memory chunk.
	///	</summary>
	unsigned char * m_pAllocatedMemory;

	///	<summary>
	///		Vector of DataChunks stored in this DataContainer.
	///	</summary>
	std::vector<void *> m_vecDataChunks;

	///	<summary>
	///		Vector of primitive types.
	///	</summary>
	std::vector<PRIMTYPE> m_vecType;
};

///////////////////////////////////////////////////////////////////////////////

#endif

