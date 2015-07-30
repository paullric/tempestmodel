///////////////////////////////////////////////////////////////////////////////
///
///	\file    DataArray1D.h
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

#ifndef _DATAARRAY1D_H_
#define _DATAARRAY1D_H_

///////////////////////////////////////////////////////////////////////////////

#include "Exception.h"
#include "DataChunk.h"
#include "DataType.h"

#include <sstream>
#include <iostream>
#include <cstdlib>
#include <cstring>

template <typename T>
class DataArray1D : public DataChunk {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	DataArray1D(
		DataType eDataType = DataType_Default,
		DataLocation eDataLocation = DataLocation_Default
	) :
		m_sRows(0),
		m_eDataType(eDataType),
		m_eDataLocation(eDataLocation),
		m_data(NULL)
	{ }

	///	<summary>
	///		Constructor with rows.
	///	</summary>
	DataArray1D(
		size_t sRows,
		DataType eDataType = DataType_Default,
		DataLocation eDataLocation = DataLocation_Default
	) :
		m_sRows(sRows),
		m_eDataType(eDataType),
		m_eDataLocation(eDataLocation),
		m_data(NULL)
	{ }

	///	<summary>
	///		Destructor.
	///	</summary>
	virtual ~DataArray1D() {
	}

public:
	///	<summary>
	///		Get the size of this DataChunk, in bytes.
	///	</summary>
	virtual size_t GetSize() const {
		return (m_sRows * sizeof(T));
	}

	///	<summary>
	///		Determine if this DataChunk is attached to a data array.
	///	</summary>
	virtual bool IsAttached() const {
		return (m_data != NULL);
	}

	///	<summary>
	///		Attach this DataChunk to an array of pre-allocated data.
	///	</summary>
	virtual void AttachTo(void * ptr) {
		m_data = reinterpret_cast<T *>(ptr);
	}

	///	<summary>
	///		Detach data from this DataChunk.
	///	</summary>
	virtual void Detach() {
		m_data = NULL;
	}

public:
	///	<summary>
	///		Get the number of elements in this DataArray1D.
	///	</summary>
	inline size_t GetRows() const {
		return m_sRows;
	}

	///	<summary>
	///		Set the number of elements in this DataArray1D.
	///	</summary>
	inline void SetRows(size_t sRows) const {
		if (IsAttached()) {
			_EXCEPTIONT("Attempting to SetRows on attached DataArray1D");
		}

		m_sRows = sRows;
	}

public:
	///	<summary>
	///		Get the DataLocation.
	///	</summary>
	inline DataLocation GetDataLocation() const {
		return m_eDataLocation;
	}

	///	<summary>
	///		Get the DataType.
	///	</summary>
	inline DataType GetDataType() const {
		return m_eDataType;
	}

public:
	///	<summary>
	///		Assignment operator.
	///	</summary>
	void Assign(const DataArray1D<T> & da) {

		// Check initialization status
		if (!da.IsAttached()) {
			_EXCEPTIONT("Attempting to assign unattached DataArray1D");
		}
		if (da.GetRows() != GetRows()) {
			_EXCEPTIONT("Size mismatch in assignment of DataArray1D");
		}
		if (da.m_eDataType != m_eDataType) {
			_EXCEPTIONT("DataType mismatch in assignment of DataArray1D");
		}
		if (da.m_eDataLocation != m_eDataLocation) {
			_EXCEPTIONT("DataLocation mismatch in assignment of DataArray1D");
		}

		// Copy data
		memcpy(m_data, da.m_data, m_sRows * sizeof(T));
	}

	///	<summary>
	///		Assignment operator.
	///	</summary>
	DataArray1D<T> & operator= (const DataArray1D<T> & da) {
		Assign(da);
		return (*this);
	}

public:
	///	<summary>
	///		Zero the data content of this object.
	///	</summary>
	void Zero() {

		// Check that this DataArray1D is attached to a data object
		if (!IsAttached()) {
			_EXCEPTIONT("Attempted operation on uninitialized DataArray1D");
		}

		// Set content to zero
		memset(m_data, 0, m_sRows * sizeof(T));
	}

	///	<summary>
	///		Scale data by a given constant.
	///	</summary>
	void Scale(const T & x) {

		// Check that this DataArray1D is attached to a data object
		if (!IsAttached()) {
			_EXCEPTIONT("Attempted operation on uninitialized DataArray1D");
		}

		// Scale data values
		for (size_t i = 0; i < m_sRows; i++) {
			m_data[i] *= x;
		}
	}

	///	<summary>
	///		Add a factor of the given DataArray1D to this DataArray1D.
	///	</summary>
	void AddProduct(
		const DataArray1D<T> & darr,
		const T & x
	) {

		// Check that this DataArray1D is attached to a data object
		if (!IsAttached()) {
			_EXCEPTIONT("Attempted operation on uninitialized DataArray1D");
		}
		if (!darr.IsAttached()) {
			_EXCEPTIONT("Attempted operation on uninitialized DataArray1D");
		}
		if (GetRows() != darr.GetRows()) {
			_EXCEPTIONT("Mismatch in DataArray1D size");
		}

		// Scale data values
		for (size_t i = 0; i < m_sRows; i++) {
			m_data[i] += x * darr.m_data[i];
		}
	}

public:
	///	<summary>
	///		Cast to an array.
	///	</summary>
	inline operator (T *)() {
		return m_data;
	}

	///	<summary>
	///		Cast to an array.
	///	</summary>
	inline operator const (T *)() const {
		return m_data;
	}

private:
	///	<summary>
	///		The number of rows in this DataArray1D.
	///	</summary>
	size_t m_sRows;

	///	<summary>
	///		The type of data stored in this DataArray1D.
	///	</summary>
	DataType m_eDataType;

	///	<summary>
	///		The location of data stored in this DataArray1D.
	///	</summary>
	DataLocation m_eDataLocation;

	///	<summary>
	///		A pointer to the data for this DataArray1D.
	///	</summary>
	T * m_data;

};

///////////////////////////////////////////////////////////////////////////////

#endif

