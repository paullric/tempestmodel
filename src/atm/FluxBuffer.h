///////////////////////////////////////////////////////////////////////////////
///
///	\file    FluxBuffer.h
///	\author  Paul Ullrich
///	\version February 19, 2014
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

#ifndef _FLUXBUFFER_H_
#define _FLUXBUFFER_H_

#include "DataMatrix3D.h"
#include "DataType.h"
#include "Direction.h"

#include <vector>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Class for storing horizontal fluxes through an edge.
///	</summary>
class FluxBuffer {

public:
	///	<summary>
	///		Default constructor.
	///	</summary>
	FluxBuffer();

	///	<summary>
	///		Destructor.
	///	</summary>
	virtual ~FluxBuffer() {
		Deinitialize();
	}

public:
	///	<summary>
	///		Initializer.
	///	</summary>
	void Initialize(
		DataType eDataType,
		int nComponents,
		int nRElements,
		Direction dirBlockEdge,
		int iABegin,
		int iAEnd
	);

	///	<summary>
	///		Deinitialize this data object.
	///	</summary>
	void Deinitialize();

public:
	///	<summary>
	///		Zero operator.
	///	</summary>
	inline void Zero() {
		if (!m_fInitialized) {
			_EXCEPTIONT("Attempting to operate on uninitialized data.");
		}

		m_data.Zero();
	}

public:
	///	<summary>
	///		Bracket accessor.
	///	</summary>
	inline double** operator[](int n) const {
		return m_data[n];
	}

	///	<summary>
	///		Conversion to nested array.
	///	</summary>
	inline operator double***() {
		return (double ***)(m_data);
	}

	///	<summary>
	///		Get the DataType.
	///	</summary>
	inline DataType GetDataType() const {
		return m_eDataType;
	}

	///	<summary>
	///		Get a constant reference to the data matrix.
	///	</summary>
	inline const DataMatrix3D<double> & GetDataMatrix() const {
		return m_data;
	}

	///	<summary>
	///		Get the total number of elements in the data matrix.
	///	</summary>
	inline int GetTotalElements() const {
		return 
			  m_data.GetRows()
			* m_data.GetColumns()
			* m_data.GetSubColumns();
	}

	///	<summary>
	///		Get the number of components in the data matrix.
	///	</summary>
	inline int GetComponents() const {
		return m_data.GetRows();
	}

	///	<summary>
	///		Get the number of vertical elements in the data matrix.
	///	</summary>
	inline int GetRElements() const {
		return m_data.GetColumns();
	}

	///	<summary>
	///		Get the number of horizontal elements in the data matrix.
	///	</summary>
	inline int GetAElements() const {
		return m_data.GetSubColumns();
	}

protected:
	///	<summary>
	///		Initialization flag.
	///	</summary>
	bool m_fInitialized;

	///	<summary>
	///		The type of data stored in this flux buffer.
	///	</summary>
	DataType m_eDataType;

	///	<summary>
	///		Edge of the block associated with this flux buffer.
	///	</summary>
	Direction m_dirBlockEdge;

	///	<summary>
	///		First index associated with this flux buffer.
	///	</summary>
	int m_iABegin;

	///	<summary>
	///		End index of the flux buffer.
	///	</summary>
	int m_iAEnd;

	///	<summary>
	///		The data for this flux buffer.
	///	</summary>
	DataMatrix3D<double> m_data;
};

///////////////////////////////////////////////////////////////////////////////

class FluxBufferVector : public std::vector<FluxBuffer> {

public:
	///	<summary>
	///		Deinitialize the vector: Delete all data and revert vector
	///		to empty.
	///	</summary>
	void Deinitialize() {
		for (int n = 0; n < size(); n++) {
			(*this)[n].Deinitialize();
		}
		resize(0);
	}
};

///////////////////////////////////////////////////////////////////////////////

#endif
