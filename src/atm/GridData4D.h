///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridData4D.h
///	\author  Paul Ullrich
///	\version February 25, 2013
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

#ifndef _GRIDDATA4D_H_
#define _GRIDDATA4D_H_

#include "DataMatrix4D.h"
#include "DataType.h"
#include "DataLocation.h"
#include "GridData3D.h"

#include <vector>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		4D data array storing prognostic variables in space.  This class
///		implements important arithmatic operations on the data.
///	</summary>
class GridData4D {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	GridData4D() :
		m_fInitialized(false),
		m_eDataType(DataType_Default),
		m_eDataLocation(DataLocation_Default)
	{ }

	///	<summary>
	///		Copy constructor.
	///	</summary>
	GridData4D(
		const GridData4D & statedata
	) {
		Duplicate(statedata);
	}

	///	<summary>
	///		Destructor.
	///	</summary>
	virtual ~GridData4D() {
		Deinitialize();
	}

public:
	///	<summary>
	///		Initializer.
	///	</summary>
	void Initialize(
		DataType eDataType,
		DataLocation eDataLocation,
		int nComponents,
		int nRElements,
		int nAElements,
		int nBElements,
		int nHaloElements
	);

	///	<summary>
	///		Deinitialize this data object.
	///	</summary>
	void Deinitialize() {
		m_fInitialized = false;
		m_eDataType = DataType_Default;
		m_eDataLocation = DataLocation_Default;
		m_data.Deinitialize();
	}

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

	///	<summary>
	///		Copy operator.
	///	</summary>
	inline GridData4D & operator=(const GridData4D & statedata) {
		Duplicate(statedata);

		return *this;
	}

	///	<summary>
	///		Duplicate data in statedata.
	///	</summary>
	inline void Duplicate(const GridData4D & statedata) {
		m_fInitialized = statedata.m_fInitialized;
		m_eDataType = statedata.m_eDataType;
		m_data = statedata.m_data;
	}

	///	<summary>
	///		Scale data by given constant.
	///	</summary>
	void Scale(
		double dFactor
	);

	///	<summary>
	///		Add a factor of the given state data to this state data.
	///	</summary>
	void AddProduct(
		const GridData4D & data,
		double dFactor
	);

public:
	///	<summary>
	///		Check if this object is initialized.
	///	</summary>
	bool IsInitialized() const {
		return m_fInitialized;
	}

	///	<summary>
	///		Bracket accessor.
	///	</summary>
	inline double*** operator[](int n) const {
		return m_data[n];
	}

	///	<summary>
	///		GridData3D accessor.
	///	</summary>
	inline void GetAsGridData3D(
		int n,
		GridData3D & data
	) const {
		data.Attach(
			m_eDataType,
			m_eDataLocation,
			GetRElements(),
			GetAElements(),
			GetBElements(),
			m_nHaloElements,
			m_data[n]);
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

	///	<summary>
	///		Get a constant reference to the data matrix.
	///	</summary>
	inline const DataMatrix4D<double> & GetDataMatrix() const {
		return m_data;
	}

	///	<summary>
	///		Get the total number of elements in the data matrix.
	///	</summary>
	inline int GetTotalElements() const {
		return 
			  m_data.GetSize(0)
			* m_data.GetSize(1)
			* m_data.GetSize(2)
			* m_data.GetSize(3);
	}

	///	<summary>
	///		Get the number of components per grid node in the data matrix.
	///	</summary>
	inline int GetComponents() const {
		return m_data.GetSize(0);
	}

	///	<summary>
	///		Get the number of vertical elements in the data matrix.
	///	</summary>
	inline int GetRElements() const {
		return m_data.GetSize(1);
	}

	///	<summary>
	///		Get the number of horizontal elements in the first coordinate
	///		direction in the data matrix.
	///	</summary>
	inline int GetAElements() const {
		return m_data.GetSize(2);
	}

	///	<summary>
	///		Get the number of horizontal elements in the second coordinate
	///		direction in the data matrix.
	///	</summary>
	inline int GetBElements() const {
		return m_data.GetSize(3);
	}

	///	<summary>
	///		Get the size in the specified coordinate in the data matrix.
	///	</summary>
	inline int GetSize(int dim) const {
		return m_data.GetSize(dim);
	}

	///	<summary>
	///		Get the number of halo elements in this GridData.
	///	</summary>
	inline int GetHaloElements() const {
		return m_nHaloElements;
	}

protected:
	///	<summary>
	///		Initialization flag.
	///	</summary>
	bool m_fInitialized;

	///	<summary>
	///		Location of data.
	///	</summary>
	DataLocation m_eDataLocation;

	///	<summary>
	///		The type of data stored in this data object.
	///	</summary>
	DataType m_eDataType;

	///	<summary>
	///		Number of halo elements in this data.
	///	</summary>
	int m_nHaloElements;

	///	<summary>
	///		The data for this system.
	///	</summary>
	DataMatrix4D<double> m_data;
};

///////////////////////////////////////////////////////////////////////////////

class GridData4DVector : public std::vector<GridData4D> {

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

