///////////////////////////////////////////////////////////////////////////////
///
///	\file    DataArray3D.h
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

#ifndef _DATAARRAY3D_H_
#define _DATAARRAY3D_H_

///////////////////////////////////////////////////////////////////////////////

#include "Exception.h"
#include "DataChunk.h"
#include "DataType.h"
#include "DataLocation.h"

#include <cstdlib>
#include <cstring>

template <typename T>
class DataArray3D : public DataChunk {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	DataArray3D(
		DataType eDataType = DataType_Default,
		DataLocation eDataLocation = DataLocation_Default
	) :
		m_fOwnsData(true),
		m_fOwnsPointerTree(true),
		m_eDataType(eDataType),
		m_eDataLocation(eDataLocation),
		m_data(NULL),
		m_data1D(NULL)
	{
		m_sSize[0] = 0;
		m_sSize[1] = 0;
		m_sSize[2] = 0;
	}

	///	<summary>
	///		Constructor allowing specification of size.
	///	</summary>
	DataArray3D(
		size_t sSize0,
		size_t sSize1,
		size_t sSize2,
		DataType eDataType = DataType_Default,
		DataLocation eDataLocation = DataLocation_Default,
		bool fAllocate = true
	) :
		m_fOwnsData(true),
		m_fOwnsPointerTree(true),
		m_eDataType(eDataType),
		m_eDataLocation(eDataLocation),
		m_data(NULL),
		m_data1D(NULL)
	{
		m_sSize[0] = sSize0;
		m_sSize[1] = sSize1;
		m_sSize[2] = sSize2;

		if (fAllocate) {
			Allocate();
		}
	}

	///	<summary>
	///		Copy constructor.
	///	</summary>
	DataArray3D(const DataArray3D<T> & da) :
		m_fOwnsData(true),
		m_fOwnsPointerTree(true),
		m_eDataType(DataType_Default),
		m_eDataLocation(DataLocation_Default),
		m_data(NULL)
	{
		m_sSize[0] = 0;
		m_sSize[1] = 0;
		m_sSize[2] = 0;

		Assign(da);
	}

	///	<summary>
	///		Destructor.
	///	</summary>
	virtual ~DataArray3D() {
		Detach();
		DeletePointerTree();
	}

private:
	///	<summary>
	///		Build the pointer tree for the specified data.
	///	</summary>
	void BuildPointerTree() {
		if (m_data != NULL) {
			_EXCEPTIONT("Attempting to rebuild existing pointer tree");
		}
		if (!m_fOwnsPointerTree) {
			_EXCEPTIONT("Logic error");
		}

		m_data = new T**[m_sSize[0]];
		for (size_t i = 0; i < m_sSize[0]; i++) {
			m_data[i] = new T*[m_sSize[1]];
			for (size_t j = 0; j < m_sSize[1]; j++) {
				m_data[i][j] = m_data1D
					+ (i * m_sSize[1] + j) * m_sSize[2];
			}
		}
	}

	///	<summary>
	///		Delete the pointer tree.
	///	</summary>
	void DeletePointerTree() {
		if (!m_fOwnsPointerTree) {
			m_data = NULL;
			m_fOwnsPointerTree = true;
			return;
		}
		if (m_data == NULL) {
			return;
		}

		for (int i = 0; i < m_sSize[0]; i++) {
			delete[] m_data[i];
		}
		delete[] m_data;

		m_data = NULL;
	}

	///	<summary>
	///		Attach pointer tree to 1D data.
	///	</summary>
	void AttachPointerTree() {
		if (!m_fOwnsPointerTree) {
			_EXCEPTIONT("Attempting to modify attached pointer tree");
		}
		for (size_t i = 0; i < m_sSize[0]; i++) {
			for (size_t j = 0; j < m_sSize[1]; j++) {
				m_data[i][j] = m_data1D
					+ (i * m_sSize[1] + j) * m_sSize[2];
			}
		}
	}

public:
	///	<summary>
	///		Get the size of this DataChunk, in bytes.
	///	</summary>
	virtual size_t GetByteSize() const {
		return (m_sSize[0] * m_sSize[1] * m_sSize[2] * sizeof(T));
	}

	///	<summary>
	///		Allocate data in this DataArray3D.
	///	</summary>
	void Allocate(
		size_t sSize0 = 0,
		size_t sSize1 = 0,
		size_t sSize2 = 0
	) {
		if ((!m_fOwnsData) || (!m_fOwnsPointerTree)) {
			_EXCEPTIONT("Attempting to Allocate() on attached DataArray3D");
		}

		if (sSize0 == 0) {
			sSize0 = m_sSize[0];
		}
		if (sSize1 == 0) {
			sSize1 = m_sSize[1];
		}
		if (sSize2 == 0) {
			sSize2 = m_sSize[2];
		}
		if ((sSize0 == 0) || (sSize1 == 0) || (sSize2 == 0)) {
			_EXCEPTIONT("Attempting to Allocate() zero-size DataArray3D");
		}
		if ((m_data1D == NULL) ||
		    (m_sSize[0] != sSize0) ||
		    (m_sSize[1] != sSize1) ||
		    (m_sSize[2] != sSize2)
		) {
			Detach();
			DeletePointerTree();

			m_sSize[0] = sSize0;
			m_sSize[1] = sSize1;
			m_sSize[2] = sSize2;

			m_data1D = reinterpret_cast<T *>(malloc(GetByteSize()));

			BuildPointerTree();
		}

		Zero();

		m_fOwnsData = true;
	}

	///	<summary>
	///		Set the dimension sizes of this DataArray3D.
	///	</summary>
	inline void SetSize(
		size_t sSize0,
		size_t sSize1,
		size_t sSize2
	) {
		if (IsAttached()) {
			_EXCEPTIONT("Attempting SetSize() on attached DataArray3D");
		}

		if (m_data != NULL) {
			if ((m_sSize[0] != sSize0) ||
			    (m_sSize[1] != sSize1) ||
				(m_sSize[2] != sSize2)
			) {
				DeletePointerTree();
			}
		}

		m_sSize[0] = sSize0;
		m_sSize[1] = sSize1;
		m_sSize[2] = sSize2;
	}

public:
	///	<summary>
	///		Determine if this DataChunk is attached to a data array.
	///	</summary>
	virtual bool IsAttached() const {
		return (m_data1D != NULL);
	}

	///	<summary>
	///		Attach this DataChunk to an array of pre-allocated data.
	///	</summary>
	virtual void AttachToData(void * ptr) {
		if (IsAttached()) {
			_EXCEPTIONT("Attempting AttachToData() on attached DataArray3D");
		}
		if (!m_fOwnsPointerTree) {
			_EXCEPTIONT("Attempting AttachToData() on attached DataArray3D");
		}

		m_data1D = reinterpret_cast<T *>(ptr);
		m_fOwnsData = false;

		if (m_data == NULL) {
			BuildPointerTree();
		} else {
			AttachPointerTree();
		}
	}

	///	<summary>
	///		Attach this DataArray3D to an array of pre-allocated data.
	///	</summary>
	void AttachTo(T *** ptr) {
		if (IsAttached()) {
			_EXCEPTIONT("Attempting AttachTo() on attached DataArray3D");
		}
		if (!m_fOwnsPointerTree) {
			_EXCEPTIONT("Attempting AttachTo() on attached DataArray3D");
		}

		m_data = ptr;
		m_data1D = ptr[0][0];
		m_fOwnsData = false;
		m_fOwnsPointerTree = false;
	}

	///	<summary>
	///		Detach data from this DataChunk.
	///	</summary>
	virtual void Detach() {
		if (!m_fOwnsPointerTree) {
			m_fOwnsPointerTree = true;
			m_data = NULL;
		}
		if ((m_fOwnsData) && (m_data1D != NULL)) {
			delete[] m_data1D;
		}
		m_fOwnsData = true;
		m_data1D = NULL;
	}

	///	<summary>
	///		Deallocate data from this DataChunk.
	///	</summary>
	void Deallocate() {
		if (!m_fOwnsData) {
			_EXCEPTIONT("Attempting Deallocate() on attached DataArray3D");
		}

		Detach();
		DeletePointerTree();
	}

public:
	///	<summary>
	///		Get the size of the data.
	///	</summary>
	size_t GetTotalSize() const {
		return (m_sSize[0] * m_sSize[1] * m_sSize[2]);
	}

	///	<summary>
	///		Get the size of the specified dimension.
	///	</summary>
	inline size_t GetSize(int dim) const {
		return m_sSize[dim];
	}

	///	<summary>
	///		Get the number of rows in this DataArray3D.
	///	</summary>
	inline size_t GetRows() const {
		return m_sSize[0];
	}

	///	<summary>
	///		Get the number of columns in this DataArray3D.
	///	</summary>
	inline size_t GetColumns() const {
		return m_sSize[1];
	}

	///	<summary>
	///		Get the number of subcolumns in this DataArray3D.
	///	</summary>
	inline size_t GetSubColumns() const {
		return m_sSize[2];
	}

public:
	///	<summary>
	///		Set the DataLocation.
	///	</summary>
	inline void SetDataLocation(const DataLocation & eDataLocation) {
		m_eDataLocation = eDataLocation;
	}

	///	<summary>
	///		Set the DataType.
	///	</summary>
	inline void SetDataType(const DataType & eDataType) {
		m_eDataType = eDataType;
	}

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
	void Assign(const DataArray3D<T> & da) {

		// Verify source array existence
		if (!da.IsAttached()) {
			if (!IsAttached()) {
				m_sSize[0] = da.m_sSize[0];
				m_sSize[1] = da.m_sSize[1];
				m_sSize[2] = da.m_sSize[2];

				m_eDataType = da.m_eDataType;
				m_eDataLocation = da.m_eDataLocation;
				return;
			}

			_EXCEPTIONT("Attempting to assign unattached DataArray3D\n"
				"to attached DataArray3D (undefined behavior)");
		}

		// Allocate if necessary
		if (!IsAttached()) {
			Allocate(da.m_sSize[0], da.m_sSize[1], da.m_sSize[2]);
			m_eDataType = da.m_eDataType;
			m_eDataLocation = da.m_eDataLocation;
		}
		if (IsAttached() && m_fOwnsData) {
			if ((m_sSize[0] != da.m_sSize[0]) ||
			    (m_sSize[1] != da.m_sSize[1]) ||
			    (m_sSize[2] != da.m_sSize[2])
			) {
				Deallocate();
				Allocate(
					da.m_sSize[0],
					da.m_sSize[1],
					da.m_sSize[2]);
			}

			m_eDataType = da.m_eDataType;
			m_eDataLocation = da.m_eDataLocation;
		}

		// Check initialization status
		if (da.GetRows() != GetRows()) {
			_EXCEPTIONT("Rows mismatch in assignment of DataArray3D");
		}
		if (da.GetColumns() != GetColumns()) {
			_EXCEPTIONT("Columns mismatch in assignment of DataArray3D");
		}
		if (da.GetSubColumns() != GetSubColumns()) {
			_EXCEPTIONT("Subcolumns mismatch in assignment of DataArray3D");
		}
		if (da.m_eDataType != m_eDataType) {
			_EXCEPTIONT("DataType mismatch in assignment of DataArray3D");
		}
		if (da.m_eDataLocation != m_eDataLocation) {
			_EXCEPTIONT("DataLocation mismatch in assignment of DataArray3D");
		}

		// Copy data
		memcpy(m_data1D, da.m_data1D, GetByteSize());
	}

	///	<summary>
	///		Assignment operator.
	///	</summary>
	DataArray3D<T> & operator= (const DataArray3D<T> & da) {
		Assign(da);
		return (*this);
	}

public:
	///	<summary>
	///		Zero the data content of this object.
	///	</summary>
	void Zero() {

		// Check that this DataArray3D is attached to a data object
		if (!IsAttached()) {
			_EXCEPTIONT("Attempted operation on unattached DataArray3D");
		}

		// Set content to zero
		memset(m_data1D, 0, GetByteSize());
	}

	///	<summary>
	///		Scale data by a given constant.
	///	</summary>
	void Scale(const T & x) {

		// Check that this DataArray3D is attached to a data object
		if (!IsAttached()) {
			_EXCEPTIONT("Attempted operation on unattached DataArray3D");
		}

		// Scale data values
		size_t sTotalSize = GetTotalSize();

		for (size_t i = 0; i < sTotalSize; i++) {
			m_data1D[i] *= x;
		}
	}

	///	<summary>
	///		Add a factor of the given DataArray3D to this DataArray3D.
	///	</summary>
	void AddProduct(
		const DataArray3D<T> & da,
		const T & x
	) {
		// Check that this DataArray3D is attached to a data object
		if (!IsAttached()) {
			_EXCEPTIONT("Attempted operation on unattached DataArray3D");
		}
		if (!da.IsAttached()) {
			_EXCEPTIONT("Attempted operation on unattached DataArray3D");
		}
		if (da.GetRows() != GetRows()) {
			_EXCEPTIONT("Rows mismatch in DataArray3D");
		}
		if (da.GetColumns() != GetColumns()) {
			_EXCEPTIONT("Columns mismatch in DataArray3D");
		}
		if (da.GetSubColumns() != GetSubColumns()) {
			_EXCEPTIONT("SubColumns mismatch in DataArray3D");
		}

		// Scale data values
		size_t sTotalSize = GetTotalSize();

		for (size_t i = 0; i < sTotalSize; i++) {
			m_data1D[i] += x * da.m_data1D[i];
		}
	}

public:
	///	<summary>
	///		Get the data.
	///	</summary>
	inline T*** GetData() const {
		return m_data;
	}

	///	<summary>
	///		Cast to an array.
	///	</summary>
	inline operator T***() const {
		return m_data;
	}

private:
	///	<summary>
	///		A flag indicating this array owns its data.
	///	</summary>
	bool m_fOwnsData;

	///	<summary>
	///		A flag indicating this array owns its pointer tree.
	///	</summary>
	bool m_fOwnsPointerTree;

	///	<summary>
	///		The size of each dimension of this DataArray3D.
	///	</summary>
	size_t m_sSize[3];

	///	<summary>
	///		The type of data stored in this DataArray3D.
	///	</summary>
	DataType m_eDataType;

	///	<summary>
	///		The location of data stored in this DataArray3D.
	///	</summary>
	DataLocation m_eDataLocation;

	///	<summary>
	///		The top level pointer in pointer tree.
	///	</summary>
	T *** m_data;

	///	<summary>
	///		A pointer to the data for this DataArray3D.
	///	</summary>
	T * m_data1D;
};

///////////////////////////////////////////////////////////////////////////////

#endif

