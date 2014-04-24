///////////////////////////////////////////////////////////////////////////////
///
///	\file    FluxBuffer.cpp
///	\author  Paul Ullrich
///	\version February 19, 2013
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

#include "FluxBuffer.h"

///////////////////////////////////////////////////////////////////////////////

FluxBuffer::FluxBuffer() :
	m_fInitialized(false),
	m_eDataType(DataType_Default),
	m_dirBlockEdge(Direction_Unreachable),
	m_iABegin(0),
	m_iAEnd(0)
{ }

///////////////////////////////////////////////////////////////////////////////

void FluxBuffer::Deinitialize() {
	m_fInitialized = false;
	m_eDataType = DataType_Default;
	m_dirBlockEdge = Direction_Unreachable;
	m_iABegin = 0;
	m_iAEnd = 0;
	m_data.Deinitialize();
}

///////////////////////////////////////////////////////////////////////////////

void FluxBuffer::Initialize(
	DataType eDataType,
	int nComponents,
	int nRElements,
	Direction dirBlockEdge,
	int iABegin,
	int iAEnd
) {
	m_eDataType = eDataType;
	if (m_eDataType == DataType_All) {
		_EXCEPTIONT("A specific DataType must be used for data objects.");
	}

	m_dirBlockEdge = dirBlockEdge;

	m_iABegin = iABegin;
	m_iAEnd = iAEnd;

	m_data.Initialize(nComponents, nRElements, iAEnd - iABegin);

	m_fInitialized = true;
}

///////////////////////////////////////////////////////////////////////////////

