///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridData3D.cpp
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

#include "GridData3D.h"

///////////////////////////////////////////////////////////////////////////////

void GridData3D::Initialize(
	DataType eDataType,
	DataLocation eDataLocation,
	int nRElements,
	int nAElements,
	int nBElements,
	int nHaloElements
) {
	m_eDataType = eDataType;
	if (m_eDataType == DataType_All) {
		_EXCEPTIONT("A specific DataType must be used for data objects.");
	}

	m_eDataLocation = eDataLocation;
	if (m_eDataLocation == DataLocation_None) {
		_EXCEPTIONT("A specific DataLocation must be used for data objects.");
	}

	m_nHaloElements = nHaloElements;

	if (m_eDataLocation == DataLocation_Node) {
		m_data.Initialize(nRElements, nAElements, nBElements);
	} else if (m_eDataLocation == DataLocation_AEdge) {
		m_data.Initialize(nRElements, nAElements+1, nBElements);
	} else if (m_eDataLocation == DataLocation_BEdge) {
		m_data.Initialize(nRElements, nAElements, nBElements+1);
	} else if (m_eDataLocation == DataLocation_REdge) {
		m_data.Initialize(nRElements+1, nAElements, nBElements);
	} else {
		_EXCEPTIONT("Invalid DataLocation");
	}

	m_fInitialized = true;
}

///////////////////////////////////////////////////////////////////////////////

void GridData3D::Attach(
	DataType eDataType,
	DataLocation eDataLocation,
	int nRElements,
	int nAElements,
	int nBElements,
	int nHaloElements,
	double *** data
) {
	// Deinitialize existing data
	Deinitialize();

	m_fInitialized = true;
	m_eDataLocation = eDataLocation;
	m_eDataType = eDataType;
	m_nHaloElements = nHaloElements;
	m_data.Attach(nRElements, nAElements, nBElements, data);
}

///////////////////////////////////////////////////////////////////////////////

void GridData3D::Scale(
	double dFactor
) {
	int k;
	int i;
	int j;

	// Scale interior data
	for (k = 0; k < m_data.GetSize(1); k++) {
	for (i = 0; i < m_data.GetSize(2); i++) {
	for (j = 0; j < m_data.GetSize(3); j++) {
		m_data[k][i][j] *= dFactor;
	}
	}
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridData3D::AddProduct(
	const GridData3D & data,
	double dFactor
) {
	int k;
	int i;
	int j;

	// Check sizes
	if ((m_data.GetSize(0) != data.GetSize(0)) ||
		(m_data.GetSize(1) != data.GetSize(1)) ||
		(m_data.GetSize(2) != data.GetSize(2))
	) {
		_EXCEPTIONT("Incompatible GridData3D objects.");
	}

	// Scale interior data
	for (k = 0; k < m_data.GetSize(0); k++) {
	for (i = 0; i < m_data.GetSize(1); i++) {
	for (j = 0; j < m_data.GetSize(2); j++) {
		m_data[k][i][j] +=
			dFactor * data[k][i][j];
	}
	}
	}
}

///////////////////////////////////////////////////////////////////////////////

