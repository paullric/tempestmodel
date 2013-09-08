///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridData3D.h
///	\author  Paul Ullrich
///	\version August 5, 2013
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

#ifndef _GRIDDATA3D_H_
#define _GRIDDATA3D_H_

#include "DataMatrix3D.h"
#include "DataType.h"
#include "DataLocation.h"

#include <vector>

///////////////////////////////////////////////////////////////////////////////

class GridData4D;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		3D data array storing variables in space.  This class
///		implements important arithmatic operations on the data.
///	</summary>
class GridData3D {
	friend class GridData4D;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	GridData3D() :
		m_fInitialized(false),
		m_eDataType(DataType_Default),
		m_eDataLocation(DataLocation_Default)
	{ }

	///	<summary>
	///		Copy constructor.
	///	</summary>
	GridData3D(
		const GridData3D & statedata
	) {
		Duplicate(statedata);
	}

	///	<summary>
	///		Destructor.
	///	</summary>
	virtual ~GridData3D() {
		Deinitialize();
	}

public:
	///	<summary>
	///		Initializer.
	///	</summary>
	void Initialize(
		DataType eDataType,
		DataLocation eDataLocation,
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

private:
	///	<summary>
	///		Attach this GridData3D to an existing array.  This function is
	///		typically called from GridData4D.
	///	</summary>
	void Attach(
		DataType eDataType,
		DataLocation eDataLocation,
		int nRElements,
		int nAElements,
		int nBElements,
		int nHaloElements,
		double *** data
	);

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
	inline GridData3D & operator=(const GridData3D & statedata) {
		Duplicate(statedata);

		return *this;
	}

	///	<summary>
	///		Duplicate data in statedata.
	///	</summary>
	inline void Duplicate(const GridData3D & statedata) {
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
		const GridData3D & data,
		double dFactor
	);

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
	///		Get the number of vertical elements in the data matrix.
	///	</summary>
	inline int GetRElements() const {
		return m_data.GetRows();
	}

	///	<summary>
	///		Get the number of horizontal elements in the first coordinate
	///		direction in the data matrix.
	///	</summary>
	inline int GetAElements() const {
		return m_data.GetColumns();
	}

	///	<summary>
	///		Get the number of horizontal elements in the second coordinate
	///		direction in the data matrix.
	///	</summary>
	inline int GetBElements() const {
		return m_data.GetSubColumns();
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
	///		The type of data stored in this data object.
	///	</summary>
	DataType m_eDataType;

	///	<summary>
	///		Location of data.
	///	</summary>
	DataLocation m_eDataLocation;

	///	<summary>
	///		Number of halo elements in this data.
	///	</summary>
	int m_nHaloElements;

	///	<summary>
	///		The data for this system.
	///	</summary>
	DataMatrix3D<double> m_data;
};

///////////////////////////////////////////////////////////////////////////////

class GridData3DVector : public std::vector<GridData3D> {

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

class GridData3DVectorVector : public std::vector<GridData3DVector> {

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

