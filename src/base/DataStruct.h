///////////////////////////////////////////////////////////////////////////////
///
///	\file    DataStruct.h
///	\author  Paul Ullrich
///	\version January 11, 2016
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

#ifndef _DATASTRUCT_H_
#define _DATASTRUCT_H_

///////////////////////////////////////////////////////////////////////////////

#include "Exception.h"
#include "DataChunk.h"

#include <cstdlib>
#include <cstring>

template <typename T>
class DataStruct : public DataChunk {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	DataStruct() :
		m_fOwnsData(true),
		m_data(NULL)
	{ }

	///	<summary>
	///		Copy constructor.
	///	</summary>
	DataStruct(const DataStruct<T> & da) :
		m_fOwnsData(true),
		m_data(NULL)
	{
		Assign(da);
	}

	///	<summary>
	///		Destructor.
	///	</summary>
	virtual ~DataStruct() {
		Detach();
	}

	///	<summary>
	///		Get the size of the data, in bytes.
	///	</summary>
	virtual size_t GetByteSize() const {

		// Verify data aligns on word boundaries
		if (sizeof(T) % sizeof(size_t) == 0) {
			return (sizeof(T));
		} else {
			return (sizeof(T) / sizeof(size_t) + 1) * sizeof(size_t);
		}
	}

	///	<summary>
	///		Allocate data in this DataStruct.
	///	</summary>
	void Allocate() {
		if (!m_fOwnsData) {
			_EXCEPTIONT("Attempting to Allocate() on attached DataStruct");
		}

		if (m_data == NULL) {
			Detach();

			m_data = reinterpret_cast<T *>(malloc(GetByteSize()));
		}

		Zero();

		m_fOwnsData = true;
	}

public:
	///	<summary>
	///		Determine if this DataChunk is attached to a data array.
	///	</summary>
	virtual bool IsAttached() const {
		return (m_data != NULL);
	}

	///	<summary>
	///		Attach this DataChunk to an array of pre-allocated data.
	///	</summary>
	virtual void AttachToData(void * ptr) {
		if (IsAttached()) {
			_EXCEPTIONT("Attempting to attach already attached DataStruct");
		}

		m_data = reinterpret_cast<T *>(ptr);
		m_fOwnsData = false;
	}

	///	<summary>
	///		Detach data from this DataChunk.
	///	</summary>
	virtual void Detach() {
		if ((m_fOwnsData) && (m_data != NULL)) {
			delete[] m_data;
		}
		m_fOwnsData = true;
		m_data = NULL;
	}

	///	<summary>
	///		Deallocate data from this DataChunk.
	///	</summary>
	void Deallocate() {
		if (!m_fOwnsData) {
			_EXCEPTIONT("Attempting to Deallocate an attached DataStruct");
		}

		Detach();
	}

public:
	///	<summary>
	///		Get the size of the data.
	///	</summary>
	size_t GetTotalSize() const {
		return (1);
	}

public:
	///	<summary>
	///		Assignment operator.
	///	</summary>
	void Assign(const DataStruct<T> & da) {

		// Verify source array existence
		if (!da.IsAttached()) {
			if (!IsAttached()) {
				return;
			}

			_EXCEPTIONT("Attempting to assign unattached DataStruct\n"
				"to attached DataStruct (undefined behavior)");
		}

		// Allocate if necessary
		if (!IsAttached()) {
			Allocate();
		}

		// Copy data
		(*m_data) = (*(da.m_data));
	}

	///	<summary>
	///		Assignment operator.
	///	</summary>
	DataStruct<T> & operator= (const DataStruct<T> & da) {
		Assign(da);
		return (*this);
	}

public:
	///	<summary>
	///		Zero the data content of this object.
	///	</summary>
	void Zero() {

		// Check that this DataStruct is attached to a data object
		if (!IsAttached()) {
			_EXCEPTIONT("Attempted operation on uninitialized DataStruct");
		}

		// Set content to zero
		memset(m_data, 0, sizeof(T));
	}

	///	<summary>
	///		Scale data by a given constant.
	///	</summary>
	void Scale(const T & x) {

		// Check that this DataStruct is attached to a data object
		if (!IsAttached()) {
			_EXCEPTIONT("Attempted operation on unattached DataStruct");
		}

		// Scale data values
		(*m_data) *= x;
	}

	///	<summary>
	///		Add a factor of the given DataStruct to this DataStruct.
	///	</summary>
	void AddProduct(
		const DataStruct<T> & da,
		const T & x
	) {
		// Check that this DataStruct is attached to a data object
		if (!IsAttached()) {
			_EXCEPTIONT("Attempted operation on unattached DataStruct");
		}
		if (!da.IsAttached()) {
			_EXCEPTIONT("Attempted operation on unattached DataStruct");
		}

		// Scale data values
		(*m_data) += x * (*(da.m_data));
	}

public:
	///	<summary>
	///		Cast to a primitive.
	///	</summary>
	inline operator T() const {
		return (*m_data);
	}

	///	<summary>
	///		Assignment operator.
	///	</summary>
	inline DataStruct<T> & operator=(const T & t) {
		if (!IsAttached()) {
			_EXCEPTIONT("Attempted operation on unattached DataStruct");
		}
		(*m_data) = t;
		return (*this);
	}

private:
	///	<summary>
	///		A flag indicating this array owns its data.
	///	</summary>
	bool m_fOwnsData;

	///	<summary>
	///		A pointer to the data for this DataStruct.
	///	</summary>
	T * m_data;
};

///////////////////////////////////////////////////////////////////////////////

#endif

