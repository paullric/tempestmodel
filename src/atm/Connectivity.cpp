///////////////////////////////////////////////////////////////////////////////
///
///	\file    Connectivity.h
///	\author  Paul Ullrich
///	\version January 10, 2017
///
///	<remarks>
///		Copyright 2000-2017 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#include "Connectivity.h"

#include "Grid.h"
#include "GridPatch.h"
#include "Model.h"
#include "EquationSet.h"

///////////////////////////////////////////////////////////////////////////////
// ExchangeBuffer
///////////////////////////////////////////////////////////////////////////////

void ExchangeBuffer::Reset() {
	if (m_dSendBuffer.GetByteSize() < sizeof(MessageHeader)) {
		_EXCEPTION1("Invalid ExchangeBuffer send buffer (%i)",
			m_dSendBuffer.GetByteSize());
	}
	if (sizeof(MessageHeader) % sizeof(double) != 0) {
		_EXCEPTIONT("sizeof(MessageHeader) % sizeof(double) != 0");
	}

	// Store ExchangeBuffer::MessageHeader
	GetSendMessageHeader((MessageHeader *)(&(m_dSendBuffer[0])));

	// Set read and write indices to just beyond the MessageHeader
	m_ixRecvBuffer = sizeof(MessageHeader) / sizeof(double);
	m_ixSendBuffer = sizeof(MessageHeader) / sizeof(double);
}

///////////////////////////////////////////////////////////////////////////////

void ExchangeBuffer::Pack(
	const DataArray3D<double> & data
) {
	const size_t sRElements = data.GetSize(0);
	const size_t sAElements = data.GetSize(1);
	const size_t sBElements = data.GetSize(2);

	// Check matrix bounds
	if (((m_dir == Direction_Right) || (m_dir == Direction_Left)) &&
		(m_ixSecond > sBElements)
	) {
		_EXCEPTIONT("GridData / ExteriorNeighbor inconsistency.");
	}
	if (((m_dir == Direction_Top) || (m_dir == Direction_Bottom)) &&
		(m_ixSecond > sAElements)
	) {
		_EXCEPTIONT("GridData / ExteriorNeighbor inconsistency.");
	}

	// Index for halo elements along boundary
	int ixBoundaryBegin;
	int ixBoundaryEnd;

	int ixABoundaryBegin;
	int ixABoundaryEnd;
	int ixBBoundaryBegin;
	int ixBBoundaryEnd;

	// Maximum index for radial elements
	int nVarRElements;

	// Pack data to send right
	if (m_dir == Direction_Right) {
		ixBoundaryBegin = sAElements - 2 * m_sHaloElements;
		ixBoundaryEnd   = sAElements - m_sHaloElements;

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixBoundaryEnd - ixBoundaryBegin)
			* (m_ixSecond - m_ixFirst);

		if (m_ixSendBuffer + nTotalValues > m_dSendBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in SendBuffer for operation.");
		}

		// Pack the SendBuffer
		if (m_fReverseDirection) {
			for (int k = 0; k < sRElements; k++) {
			for (int i = ixBoundaryBegin; i < ixBoundaryEnd; i++) {
			for (int j = m_ixSecond-1; j >= m_ixFirst; j--) {
				m_dSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}

		} else {
			for (int k = 0; k < sRElements; k++) {
			for (int i = ixBoundaryBegin; i < ixBoundaryEnd; i++) {
			for (int j = m_ixFirst; j < m_ixSecond; j++) {
				m_dSendBuffer[m_ixSendBuffer++] = data[k][i][j];	
			}
			}
			}
		}
	
	// Pack data to send topward
	} else if (m_dir == Direction_Top) {
		ixBoundaryBegin = sBElements - 2 * m_sHaloElements;
		ixBoundaryEnd   = sBElements - m_sHaloElements;

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixBoundaryEnd - ixBoundaryBegin)
			* (m_ixSecond - m_ixFirst);

		if (m_ixSendBuffer + nTotalValues > m_dSendBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in SendBuffer for operation.");
		}

		// Pack the SendBuffer
		if (m_fReverseDirection) {
			for (int k = 0; k < sRElements; k++) {
			for (int j = ixBoundaryBegin; j < ixBoundaryEnd; j++) {
			for (int i = m_ixSecond-1; i >= m_ixFirst; i--) {
				m_dSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}

		} else {
			for (int k = 0; k < sRElements; k++) {
			for (int j = ixBoundaryBegin; j < ixBoundaryEnd; j++) {
			for (int i = m_ixFirst; i < m_ixSecond; i++) {
				m_dSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}
		}

	// Pack data to send left
	} else if (m_dir == Direction_Left) {
		ixBoundaryBegin = m_sHaloElements;
		ixBoundaryEnd   = 2 * m_sHaloElements;

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixBoundaryEnd - ixBoundaryBegin)
			* (m_ixSecond - m_ixFirst);

		if (m_ixSendBuffer + nTotalValues > m_dSendBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in SendBuffer for operation.");
		}

		// Pack the SendBuffer
		if (m_fReverseDirection) {
			for (int k = 0; k < sRElements; k++) {
			for (int i = ixBoundaryEnd-1; i >= ixBoundaryBegin; i--) {
			for (int j = m_ixSecond-1; j >= m_ixFirst; j--) {
				m_dSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}

		} else {
			for (int k = 0; k < sRElements; k++) {
			for (int i = ixBoundaryEnd-1; i >= ixBoundaryBegin; i--) {
			for (int j = m_ixFirst; j < m_ixSecond; j++) {
				m_dSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}
		}

	// Pack data to send bottomward
	} else if (m_dir == Direction_Bottom) {
		ixBoundaryBegin = m_sHaloElements;
		ixBoundaryEnd   = 2 * m_sHaloElements;

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixBoundaryEnd - ixBoundaryBegin)
			* (m_ixSecond - m_ixFirst);

		if (m_ixSendBuffer + nTotalValues > m_dSendBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in SendBuffer for operation.");
		}

		// Pack the SendBuffer
		if (m_fReverseDirection) {
			for (int k = 0; k < sRElements; k++) {
			for (int j = ixBoundaryEnd-1; j >= ixBoundaryBegin; j--) {
			for (int i = m_ixSecond-1; i >= m_ixFirst; i--) {
				m_dSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}

		} else {
			for (int k = 0; k < sRElements; k++) {
			for (int j = ixBoundaryEnd-1; j >= ixBoundaryBegin; j--) {
			for (int i = m_ixFirst; i < m_ixSecond; i++) {
				m_dSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}
		}

	// Pack data to send toprightward 
	} else if (m_dir == Direction_TopRight) {
		ixABoundaryBegin = m_ixFirst - m_sHaloElements + 1;
		ixABoundaryEnd   = m_ixFirst + 1;
		ixBBoundaryBegin = m_ixSecond - m_sHaloElements + 1;
		ixBBoundaryEnd   = m_ixSecond + 1;
 
		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixABoundaryEnd - ixABoundaryBegin)
			* (ixBBoundaryEnd - ixBBoundaryBegin);

		if (m_ixSendBuffer + nTotalValues > m_dSendBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in SendBuffer for operation.");
		}

		// Pack the SendBuffer
		if (m_fReverseDirection) {
			for (int k = 0; k < sRElements; k++) {
			for (int i = ixABoundaryBegin; i < ixABoundaryEnd; i++) {
			for (int j = ixBBoundaryBegin; j < ixBBoundaryEnd; j++) {
				m_dSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}

		} else {
			for (int k = 0; k < sRElements; k++) {
			for (int j = ixBBoundaryBegin; j < ixBBoundaryEnd; j++) {
			for (int i = ixABoundaryBegin; i < ixABoundaryEnd; i++) {
				m_dSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}
		}

	// Pack data to send topleftward 
	} else if (m_dir == Direction_TopLeft) {
		ixABoundaryBegin = m_ixFirst;
		ixABoundaryEnd   = m_ixFirst + m_sHaloElements;
		ixBBoundaryBegin = m_ixSecond - m_sHaloElements + 1;
		ixBBoundaryEnd   = m_ixSecond + 1;
 
		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixABoundaryEnd - ixABoundaryBegin)
			* (ixBBoundaryEnd - ixBBoundaryBegin);

		if (m_ixSendBuffer + nTotalValues > m_dSendBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in SendBuffer for operation.");
		}

		// Pack the SendBuffer
		if (m_fReverseDirection) {
			for (int k = 0; k < sRElements; k++) {
			for (int i = ixABoundaryEnd-1; i >= ixABoundaryBegin; i--) {
			for (int j = ixBBoundaryBegin; j < ixBBoundaryEnd; j++) {
				m_dSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}

		} else {
			for (int k = 0; k < sRElements; k++) {
			for (int j = ixBBoundaryBegin; j < ixBBoundaryEnd; j++) {
			for (int i = ixABoundaryEnd-1; i >= ixABoundaryBegin; i--) {
				m_dSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}
		}

	// Pack data to send bottomleftward 
	} else if (m_dir == Direction_BottomLeft) {
		ixABoundaryBegin = m_ixFirst;
		ixABoundaryEnd   = m_ixFirst + m_sHaloElements;
		ixBBoundaryBegin = m_ixSecond;
		ixBBoundaryEnd   = m_ixSecond + m_sHaloElements;
 
		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixABoundaryEnd - ixABoundaryBegin)
			* (ixBBoundaryEnd - ixBBoundaryBegin);

		if (m_ixSendBuffer + nTotalValues > m_dSendBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in SendBuffer for operation.");
		}

		// Pack the SendBuffer
		if (m_fReverseDirection) {
			for (int k = 0; k < sRElements; k++) {
			for (int i = ixABoundaryEnd-1; i >= ixABoundaryBegin; i--) {
			for (int j = ixBBoundaryEnd-1; j >= ixBBoundaryBegin; j--) {
				m_dSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}

		} else {
			for (int k = 0; k < sRElements; k++) {
			for (int j = ixBBoundaryEnd-1; j >= ixBBoundaryBegin; j--) {
			for (int i = ixABoundaryEnd-1; i >= ixABoundaryBegin; i--) {
				m_dSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}
		}

	// Pack data to send bottomrightward 
	} else if (m_dir == Direction_BottomRight) {
		ixABoundaryBegin = m_ixFirst - m_sHaloElements + 1;
		ixABoundaryEnd   = m_ixFirst + 1;
		ixBBoundaryBegin = m_ixSecond;
		ixBBoundaryEnd   = m_ixSecond + m_sHaloElements;
 
		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixABoundaryEnd - ixABoundaryBegin)
			* (ixBBoundaryEnd - ixBBoundaryBegin);

		if (m_ixSendBuffer + nTotalValues > m_dSendBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in SendBuffer for operation.");
		}

		// Pack the SendBuffer
		if (m_fReverseDirection) {
			for (int k = 0; k < sRElements; k++) {
			for (int i = ixABoundaryBegin; i < ixABoundaryEnd; i++) {
			for (int j = ixBBoundaryEnd-1; j >= ixBBoundaryBegin; j--) {
				m_dSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}

		} else {
			for (int k = 0; k < sRElements; k++) {
			for (int j = ixBBoundaryEnd-1; j >= ixBBoundaryBegin; j--) {
			for (int i = ixABoundaryBegin; i < ixABoundaryEnd; i++) {
				m_dSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}
		}

	// Invalid direction
	} else {
		_EXCEPTIONT("Invalid direction");
	}
}

///////////////////////////////////////////////////////////////////////////////

void ExchangeBuffer::Pack(
	const Grid & grid,
	const DataArray4D<double> & data
) {
	// Number of components in data
	size_t sComponents = data.GetSize(0);

	// 3D Grid Data
	DataArray3D<double> data3D;
	data3D.SetSize(
		data.GetSize(1),
		data.GetSize(2),
		data.GetSize(3));

	// For state data exclude non-collacted data points
	if (data.GetDataType() == DataType_State) {
		// List of variable indices to send
		// - exclude variables which are not-collocated with this data structure
		for (int c = 0; c < sComponents; c++) {
			if (grid.GetVarLocation(c) != data.GetDataLocation()) {
				continue;
			}
			data3D.AttachToData(const_cast<double*>(&(data[c][0][0][0])));
			Pack(data3D);
			data3D.Detach();
		}

	// Send everything
	} else {
		for (int c = 0; c < sComponents; c++) {
			data3D.AttachToData(const_cast<double*>(&(data[c][0][0][0])));
			Pack(data3D);
			data3D.Detach();
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void ExchangeBuffer::Unpack(
	DataArray3D<double> & data
) {
	const size_t sRElements = data.GetSize(0);
	const size_t sAElements = data.GetSize(1);
	const size_t sBElements = data.GetSize(2);

	// Index for halo elements along boundary
	int ixBoundaryBegin;
	int ixBoundaryEnd;

	int ixABoundaryBegin;
	int ixABoundaryEnd;

	int ixBBoundaryBegin;
	int ixBBoundaryEnd;

	// Unpack data from right
	if (m_dir == Direction_Right) {
		ixBoundaryBegin = sAElements - m_sHaloElements;
		ixBoundaryEnd   = sAElements;

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixBoundaryEnd - ixBoundaryBegin)
			* (m_ixSecond - m_ixFirst);

		if (m_ixRecvBuffer + nTotalValues > m_dRecvBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in RecvBuffer for operation.");
		}

		// Unpack data
		for (int k = 0; k < sRElements; k++) {
		for (int i = ixBoundaryEnd-1; i >= ixBoundaryBegin; i--) {
		for (int j = m_ixFirst; j < m_ixSecond; j++) {
			data[k][i][j] = m_dRecvBuffer[m_ixRecvBuffer++];
		}
		}
		}

	// Unpack data from top
	} else if (m_dir == Direction_Top) {
		ixBoundaryBegin = sBElements - m_sHaloElements;
		ixBoundaryEnd   = sBElements;

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixBoundaryEnd - ixBoundaryBegin)
			* (m_ixSecond - m_ixFirst);

		if (m_ixRecvBuffer + nTotalValues > m_dRecvBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in RecvBuffer for operation.");
		}

		// Unpack data
		for (int k = 0; k < sRElements; k++) {
		for (int j = ixBoundaryEnd-1; j >= ixBoundaryBegin; j--) {
		for (int i = m_ixFirst; i < m_ixSecond; i++) {
			data[k][i][j] = m_dRecvBuffer[m_ixRecvBuffer++];
		}
		}
		}

	// Unpack data from left
	} else if (m_dir == Direction_Left) {
		ixBoundaryBegin = 0;
		ixBoundaryEnd   = m_sHaloElements;

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixBoundaryEnd - ixBoundaryBegin)
			* (m_ixSecond - m_ixFirst);

		if (m_ixRecvBuffer + nTotalValues > m_dRecvBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in RecvBuffer for operation.");
		}

		// Unpack data
		for (int k = 0; k < sRElements; k++) {
		for (int i = ixBoundaryBegin; i < ixBoundaryEnd; i++) {
		for (int j = m_ixFirst; j < m_ixSecond; j++) {
			data[k][i][j] = m_dRecvBuffer[m_ixRecvBuffer++];
		}
		}
		}

	// Unpack data from bottom
	} else if (m_dir == Direction_Bottom) {
		ixBoundaryBegin = 0;
		ixBoundaryEnd   = m_sHaloElements;

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixBoundaryEnd - ixBoundaryBegin)
			* (m_ixSecond - m_ixFirst);

		if (m_ixRecvBuffer + nTotalValues > m_dRecvBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in RecvBuffer for operation.");
		}

		// Unpack data
		for (int k = 0; k < sRElements; k++) {
		for (int j = ixBoundaryBegin; j < ixBoundaryEnd; j++) {
		for (int i = m_ixFirst; i < m_ixSecond; i++) {
			data[k][i][j] = m_dRecvBuffer[m_ixRecvBuffer++];
		}
		}
		}

	// Unpack data from top-right
	} else if (m_dir == Direction_TopRight) {
		ixABoundaryBegin = m_ixFirst + 1;
		ixABoundaryEnd   = m_ixFirst + m_sHaloElements + 1;
		ixBBoundaryBegin = m_ixSecond + 1;
		ixBBoundaryEnd   = m_ixSecond + m_sHaloElements + 1;

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixABoundaryEnd - ixABoundaryBegin)
			* (ixBBoundaryEnd - ixBBoundaryBegin);

		if (m_ixRecvBuffer + nTotalValues > m_dRecvBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in RecvBuffer for operation.");
		}

		// Unpack data
		for (int k = 0; k < sRElements; k++) {
		for (int j = ixBBoundaryEnd-1; j >= ixBBoundaryBegin; j--) {
		for (int i = ixABoundaryEnd-1; i >= ixABoundaryBegin; i--) {
			data[k][i][j] = m_dRecvBuffer[m_ixRecvBuffer++];
		}
		}
		}

	// Unpack data from top-left
	} else if (m_dir == Direction_TopLeft) {
		ixABoundaryBegin = m_ixFirst - m_sHaloElements;
		ixABoundaryEnd   = m_ixFirst;
		ixBBoundaryBegin = m_ixSecond + 1;
		ixBBoundaryEnd   = m_ixSecond + m_sHaloElements + 1;

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixABoundaryEnd - ixABoundaryBegin)
			* (ixBBoundaryEnd - ixBBoundaryBegin);

		if (m_ixRecvBuffer + nTotalValues > m_dRecvBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in RecvBuffer for operation.");
		}

		// Unpack data
		for (int k = 0; k < sRElements; k++) {
		for (int j = ixBBoundaryEnd-1; j >= ixBBoundaryBegin; j--) {
		for (int i = ixABoundaryBegin; i < ixABoundaryEnd; i++) {
			data[k][i][j] = m_dRecvBuffer[m_ixRecvBuffer++];
		}
		}
		}

	// Unpack data from bottom-left
	} else if (m_dir == Direction_BottomLeft) {
		ixABoundaryBegin = m_ixFirst - m_sHaloElements;
		ixABoundaryEnd   = m_ixFirst;
		ixBBoundaryBegin = m_ixSecond - m_sHaloElements;
		ixBBoundaryEnd   = m_ixSecond;

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixABoundaryEnd - ixABoundaryBegin)
			* (ixBBoundaryEnd - ixBBoundaryBegin);

		if (m_ixRecvBuffer + nTotalValues > m_dRecvBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in RecvBuffer for operation.");
		}

		// Unpack data
		for (int k = 0; k < sRElements; k++) {
		for (int j = ixBBoundaryBegin; j < ixBBoundaryEnd; j++) {
		for (int i = ixABoundaryBegin; i < ixABoundaryEnd; i++) {
			data[k][i][j] = m_dRecvBuffer[m_ixRecvBuffer++];
		}
		}
		}

	// Unpack data from bottom-right
	} else if (m_dir == Direction_BottomRight) {
		ixABoundaryBegin = m_ixFirst + 1;
		ixABoundaryEnd   = m_ixFirst + m_sHaloElements + 1;
		ixBBoundaryBegin = m_ixSecond - m_sHaloElements;
		ixBBoundaryEnd   = m_ixSecond;

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixABoundaryEnd - ixABoundaryBegin)
			* (ixBBoundaryEnd - ixBBoundaryBegin);

		if (m_ixRecvBuffer + nTotalValues > m_dRecvBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in RecvBuffer for operation.");
		}

		// Unpack data
		for (int k = 0; k < sRElements; k++) {
		for (int j = ixBBoundaryBegin; j < ixBBoundaryEnd; j++) {
		for (int i = ixABoundaryEnd-1; i >= ixABoundaryBegin; i--) {
			data[k][i][j] = m_dRecvBuffer[m_ixRecvBuffer++];
		}
		}
		}

	// Invalid direction
	} else {
		_EXCEPTIONT("Invalid direction");
	}
}

///////////////////////////////////////////////////////////////////////////////

void ExchangeBuffer::Unpack(
	const Grid & grid,
	DataArray4D<double> & data
) {
	// Number of components in data
	size_t sComponents = data.GetSize(0);

	// 3D Grid Data
	DataArray3D<double> data3D;
	data3D.SetSize(
		data.GetSize(1),
		data.GetSize(2),
		data.GetSize(3));

	// List of variable indices to receive
	// - exclude variables which are not-collocated with this data structure
	if (data.GetDataType() == DataType_State) {
		for (int c = 0; c < sComponents; c++) {
			if (grid.GetVarLocation(c) != data.GetDataLocation()) {
				continue;
			}
			data3D.AttachToData(&(data[c][0][0][0]));
			Unpack(data3D);
			data3D.Detach();
		}

	// Unpack all variables
	} else {
		for (int c = 0; c < sComponents; c++) {
			data3D.AttachToData(&(data[c][0][0][0]));
			Unpack(data3D);
			data3D.Detach();
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
// ExchangeBufferRegistry
///////////////////////////////////////////////////////////////////////////////

ExchangeBufferRegistry::ExchangeBufferRegistry()
{ }

///////////////////////////////////////////////////////////////////////////////

ExchangeBufferRegistry::~ExchangeBufferRegistry() {
	for (int i = 0; i < m_vecRecvBuffers.size(); i++) {
		if (m_vecRecvBuffers[i] != NULL) {
			delete[] m_vecRecvBuffers[i];
		}
	}
	for (int i = 0; i < m_vecSendBuffers.size(); i++) {
		if (m_vecSendBuffers[i] != NULL) {
			delete[] m_vecSendBuffers[i];
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void ExchangeBufferRegistry::Register(
	const ExchangeBuffer & exbuf
) {
#ifdef TEMPEST_MPIOMP
	int nRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

	if (exbuf.m_ixSourceProcessor != nRank) {
		_EXCEPTIONT("ExchangeBuffer must be local");
	}
#endif

	// Add the ExchangeBuffer to the registry
	m_vecRegistry.push_back(exbuf);
}

///////////////////////////////////////////////////////////////////////////////

void ExchangeBufferRegistry::Allocate() {

	// Allocate already called
	if (m_vecProcessors.size() != 0) {
		_EXCEPTIONT("Allocate called on already allocated object");
	}

	// A map between processor index and exchange buffer size
	std::map<int, int> mapProcessorToBufferSize;

	// Calculate the number of bytes needed per processor
	{
		for (int m = 0; m < m_vecRegistry.size(); m++) {
			int nExchangeBufferByteSize =
				m_vecRegistry[m].GetMessageSize();

			if (m_vecRegistry[m].m_ixTargetProcessor == (-1)) {
				_EXCEPTIONT("Invalid processor for ExchangeBuffer");
			}

			std::map<int, int>::iterator iterProcs =
				mapProcessorToBufferSize.find(
					m_vecRegistry[m].m_ixTargetProcessor);

			if (iterProcs == mapProcessorToBufferSize.end()) {
				iterProcs =
					mapProcessorToBufferSize.insert(
						std::pair<int, int>(
							m_vecRegistry[m].m_ixTargetProcessor, 0)).first;
			}

			iterProcs->second += nExchangeBufferByteSize;
		}
	}

	// Allocate space per processor
	{
		std::map<int, int>::const_iterator iterProcs =
			mapProcessorToBufferSize.begin();
		for (; iterProcs != mapProcessorToBufferSize.end(); iterProcs++) {
			m_vecProcessors.push_back(iterProcs->first);
			
			char * pRecvBuffer = new char[iterProcs->second];
			char * pSendBuffer = new char[iterProcs->second];

			if (pRecvBuffer == NULL) {
				_EXCEPTIONT("Out of memory");
			}
			if (pSendBuffer == NULL) {
				_EXCEPTIONT("Out of memory");
			}

			m_vecBufferSize.push_back(iterProcs->second);
			m_vecRecvBuffers.push_back(pRecvBuffer);
			m_vecSendBuffers.push_back(pSendBuffer);
			m_vecAreRecvBuffersAttached.push_back(false);
		}
	}

	// Build lookup table for ExchangeBuffers
	{
		m_vecRegistryByProcessor.resize(m_vecProcessors.size());
		for (int m = 0; m < m_vecRegistry.size(); m++) {
			int p = 0;
			for (; p < m_vecProcessors.size(); p++) {
				const int ixTargetProcessor =
					m_vecRegistry[m].m_ixTargetProcessor;

				if (ixTargetProcessor == m_vecProcessors[p]) {
					m_vecRegistryByProcessor[p].push_back(
						&(m_vecRegistry[m]));
					break;
				}
			}
			if (p ==  m_vecProcessors.size()) {
				_EXCEPTIONT("Lookup table logic error");
			}
		}
	}

	// Assign space to send buffers for each ExchangeBuffer (receive
	// buffers are attached when the first message arrives)
	{
		std::vector<int> vecSendBufferPosition;
		vecSendBufferPosition.resize(m_vecProcessors.size());

		for (int m = 0; m < m_vecRegistry.size(); m++) {
			int nExchangeBufferByteSize =
				m_vecRegistry[m].GetMessageSize();

			if (nExchangeBufferByteSize % sizeof(double) != 0) {
				_EXCEPTIONT("Message size must be aligned at "
					"double boundaries");
			}

			int p = 0;
			for (; p < m_vecProcessors.size(); p++) {
				int ixTargetProc = m_vecRegistry[m].m_ixTargetProcessor;
				if (m_vecProcessors[p] == ixTargetProc) {
					break;
				}
			}
			if (p == m_vecProcessors.size()) {
				_EXCEPTIONT("Logic error");
			}

			m_vecRegistry[m].m_dSendBuffer.SetSize(
				nExchangeBufferByteSize / sizeof(double));
			m_vecRegistry[m].m_dSendBuffer.AttachToData(
				m_vecSendBuffers[p] + vecSendBufferPosition[p]);

			vecSendBufferPosition[p] +=
				nExchangeBufferByteSize;
		}
	}

	// Allocate number of receive and send requests
	m_vecRecvRequest.resize(m_vecProcessors.size());
	m_vecSendRequest.resize(m_vecProcessors.size());

	m_vecMessageReceived.resize(m_vecProcessors.size());
}

///////////////////////////////////////////////////////////////////////////////

void ExchangeBufferRegistry::PrepareExchange() {

	// Reset all ExchangeBuffers
	for (int r = 0; r < m_vecRegistry.size(); r++) {
		m_vecRegistry[r].Reset();
	}

#ifdef TEMPEST_MPIOMP
	// Set up asynchornous receives
	for (int p = 0; p < m_vecProcessors.size(); p++) {
		m_vecMessageReceived[p] = false;

#pragma message "Check for Exchange between processors of the same rank"
		MPI_Irecv(
			m_vecRecvBuffers[p],
			m_vecBufferSize[p],
			MPI_BYTE,
			m_vecProcessors[p],
			0,
			MPI_COMM_WORLD,
			&(m_vecRecvRequest[p]));
/*
		int nRank;
		MPI_Comm_rank(MPI_COMM_WORLD, &nRank);
		printf("On %i IRecv %i bytes from %i\n",
			nRank,
			m_vecBufferSize[p],
			m_vecProcessors[p]);
*/
	}
#endif
}

///////////////////////////////////////////////////////////////////////////////

void ExchangeBufferRegistry::Send() {

#ifdef TEMPEST_MPIOMP
	for (int p = 0; p < m_vecProcessors.size(); p++) {
		MPI_Isend(
			m_vecSendBuffers[p],
			m_vecBufferSize[p],
			MPI_BYTE,
			m_vecProcessors[p],
			0,
			MPI_COMM_WORLD,
			&(m_vecSendRequest[p]));

/*
		int nRank;
		MPI_Comm_rank(MPI_COMM_WORLD, &nRank);
		printf("On %i sending %i bytes to %i\n",
			nRank,
			m_vecBufferSize[p],
			m_vecProcessors[p]);
*/
	}
#endif
}

///////////////////////////////////////////////////////////////////////////////

void ExchangeBufferRegistry::AttachRecvBuffers(int p) {

	// Find the array of ExchangeBuffers relevant to this processor
	if ((p < 0) || (p >= m_vecRegistryByProcessor.size())) {
		_EXCEPTIONT("Processor index out of range");
	}
	std::vector<ExchangeBuffer *> & vecExchangeBufs =
		m_vecRegistryByProcessor[p];

	// Number of ExchangeBuffers needed
	int nProcExchangeBuffers = vecExchangeBufs.size();

	// Assign buffers in the correct order
	int iPosition = 0;

	ExchangeBuffer::MessageHeader * msghead =
		(ExchangeBuffer::MessageHeader *)(m_vecRecvBuffers[p]);

	int nProcExchangeBuffersAssigned = 0;
	for (;;) {
		// Search for ExchangeBuffers with a header that
		// matches the header in the message
		int b = 0;
		for (; b < nProcExchangeBuffers; b++) {
			ExchangeBuffer * pExBuf = vecExchangeBufs[b];

			ExchangeBuffer::MessageHeader exbufhead;
			pExBuf->GetRecvMessageHeader(&exbufhead);

			if (exbufhead != (*msghead)) {
				continue;
			}

			pExBuf->m_dRecvBuffer.SetSize(
				pExBuf->GetMessageSize());
			pExBuf->m_dRecvBuffer.AttachToData(
				m_vecRecvBuffers[p] + iPosition);

			iPosition += pExBuf->GetMessageSize();
			msghead = (ExchangeBuffer::MessageHeader *)
				(m_vecRecvBuffers[p] + iPosition);
			break;
		}
		if (b == nProcExchangeBuffers) {
			_EXCEPTIONT("Corresponding ExchangeBuffer not found");
		}

		nProcExchangeBuffersAssigned++;

		if (iPosition > m_vecBufferSize[p]) {
			_EXCEPTIONT("Message length does not match buffer size");
		}
		if (iPosition == m_vecBufferSize[p]) {
			break;
		}
	}

	if (nProcExchangeBuffersAssigned != nProcExchangeBuffers) {
		_EXCEPTIONT("Insufficient ExchangeBuffers in message");
	}
}

///////////////////////////////////////////////////////////////////////////////

const std::vector<ExchangeBuffer *> * ExchangeBufferRegistry::WaitReceive() {

#ifdef TEMPEST_MPIOMP
	// Receive data from exterior neighbors
	int nRecvMessageCount = 0;

	while (nRecvMessageCount != m_vecProcessors.size()) {

		nRecvMessageCount = 0;
		for (int p = 0; p < m_vecProcessors.size(); p++) {

			// Check if this message has already been received
			if (m_vecMessageReceived[p]) {
				nRecvMessageCount++;
				continue;
			}

			// Check if a message is waiting
			int fRecvWaiting;
			MPI_Status status;
			MPI_Test(&(m_vecRecvRequest[p]), &fRecvWaiting, &status);
			if (!fRecvWaiting) {
				continue;
			}

			// Message received
			m_vecMessageReceived[p] = true;
/*
			int nRank;
			MPI_Comm_rank(MPI_COMM_WORLD, &nRank);
			printf("Message received on proc %i from proc %i\n", nRank, m_vecProcessors[p]);
*/
			// Attach Recv buffers
			if (!m_vecAreRecvBuffersAttached[p]) {
				AttachRecvBuffers(p);
				m_vecAreRecvBuffersAttached[p] = true;
			}

			// Return the array of ExchangeBuffers that have been filled
			return &(m_vecRegistryByProcessor[p]);
		}
	}
#endif

	return (NULL);
}

///////////////////////////////////////////////////////////////////////////////

void ExchangeBufferRegistry::WaitSend() {
#ifdef TEMPEST_MPIOMP
	for (int r = 0; r < m_vecSendRequest.size(); r++) {
		MPI_Status status;
		MPI_Wait(&m_vecSendRequest[r], &status);
	}
#endif
}

///////////////////////////////////////////////////////////////////////////////

