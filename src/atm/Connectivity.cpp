///////////////////////////////////////////////////////////////////////////////
///
///	\file    Connectivity.h
///	\author  Paul Ullrich
///	\version March 24, 2013
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

#include "Connectivity.h"

#include "Grid.h"
#include "GridPatch.h"
#include "Model.h"
#include "EquationSet.h"

///////////////////////////////////////////////////////////////////////////////
// Neighbor
///////////////////////////////////////////////////////////////////////////////

void Neighbor::AllocateBuffers() {
	m_vecSendBuffer.Allocate(
		  m_meta.sBoundarySize
		* m_meta.sMaxRElements
		* m_meta.sHaloElements
		* m_meta.sComponents);

	m_vecRecvBuffer.Allocate(
		  m_meta.sBoundarySize
		* m_meta.sMaxRElements
		* m_meta.sHaloElements
		* m_meta.sComponents);
}

///////////////////////////////////////////////////////////////////////////////

size_t Neighbor::GetByteSize() const {
	return (
		m_meta.sBoundarySize
		* m_meta.sMaxRElements
		* m_meta.sHaloElements
		* m_meta.sComponents
		* sizeof(double));
}

///////////////////////////////////////////////////////////////////////////////

bool Neighbor::CheckReceive() {

	// Check if message already received and processed
	if (m_fComplete) {
		return false;
	}

#ifdef USE_MPI
	// Test receive request
	int fRecvWaiting;
	MPI_Status status;
	MPI_Test(&m_reqRecv, &fRecvWaiting, &status);
	if (!fRecvWaiting) {
		return false;
	}
#endif

	// A message has been received
	return true;
}

///////////////////////////////////////////////////////////////////////////////
// ExteriorNeighbor
///////////////////////////////////////////////////////////////////////////////

void ExteriorNeighbor::PrepareExchange(
	const Grid & grid
) {
	// Call up the stack
	Neighbor::PrepareExchange(grid);

#ifdef USE_MPI
	// Information for receive
	int iProcessor = grid.GetPatchProcessor(m_key.ixTargetPatch);

	// Tag of the received message
	int nTag = (m_key.ixTargetPatch << 16) + (m_key.ixSourcePatch << 4) + (int)(m_key.dir);

	// Prepare an asynchronous receive
	MPI_Irecv(
		&(m_vecRecvBuffer[0]),
		m_vecRecvBuffer.GetRows(),
		MPI_DOUBLE,
		iProcessor,
		nTag,
		MPI_COMM_WORLD,
		&m_reqRecv);
#endif
}

///////////////////////////////////////////////////////////////////////////////

void ExteriorNeighbor::Pack(
	const Grid & grid,
	const DataArray3D<double> & data
) {
	const size_t sRElements = data.GetSize(0);
	const size_t sAElements = data.GetSize(1);
	const size_t sBElements = data.GetSize(2);

	// Check matrix bounds
	if (((m_key.dir == Direction_Right) || (m_key.dir == Direction_Left)) &&
		(m_meta.ixSecond > sBElements)
	) {
		_EXCEPTIONT("GridData / ExteriorNeighbor inconsistency.");
	}
	if (((m_key.dir == Direction_Top) || (m_key.dir == Direction_Bottom)) &&
		(m_meta.ixSecond > sAElements)
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
	if (m_key.dir == Direction_Right) {
		ixBoundaryBegin = sAElements - 2 * m_meta.sHaloElements;
		ixBoundaryEnd   = sAElements - m_meta.sHaloElements;

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixBoundaryEnd - ixBoundaryBegin)
			* (m_meta.ixSecond - m_meta.ixFirst);

		// Pack the SendBuffer
		if (m_fReverseDirection) {
			for (int k = 0; k < sRElements; k++) {
			for (int i = ixBoundaryBegin; i < ixBoundaryEnd; i++) {
			for (int j = m_meta.ixSecond-1; j >= m_meta.ixFirst; j--) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}

		} else {
			for (int k = 0; k < sRElements; k++) {
			for (int i = ixBoundaryBegin; i < ixBoundaryEnd; i++) {
			for (int j = m_meta.ixFirst; j < m_meta.ixSecond; j++) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];	
			}
			}
			}
		}
	
	// Pack data to send topward
	} else if (m_key.dir == Direction_Top) {
		ixBoundaryBegin = sBElements - 2 * m_meta.sHaloElements;
		ixBoundaryEnd   = sBElements - m_meta.sHaloElements;

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixBoundaryEnd - ixBoundaryBegin)
			* (m_meta.ixSecond - m_meta.ixFirst);

		if (m_ixSendBuffer + nTotalValues > m_vecSendBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in SendBuffer for operation.");
		}

		// Pack the SendBuffer
		if (m_fReverseDirection) {
			for (int k = 0; k < sRElements; k++) {
			for (int j = ixBoundaryBegin; j < ixBoundaryEnd; j++) {
			for (int i = m_meta.ixSecond-1; i >= m_meta.ixFirst; i--) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}

		} else {
			for (int k = 0; k < sRElements; k++) {
			for (int j = ixBoundaryBegin; j < ixBoundaryEnd; j++) {
			for (int i = m_meta.ixFirst; i < m_meta.ixSecond; i++) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}
		}

	// Pack data to send left
	} else if (m_key.dir == Direction_Left) {
		ixBoundaryBegin = m_meta.sHaloElements;
		ixBoundaryEnd   = 2 * m_meta.sHaloElements;

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixBoundaryEnd - ixBoundaryBegin)
			* (m_meta.ixSecond - m_meta.ixFirst);

		if (m_ixSendBuffer + nTotalValues > m_vecSendBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in SendBuffer for operation.");
		}

		// Pack the SendBuffer
		if (m_fReverseDirection) {
			for (int k = 0; k < sRElements; k++) {
			for (int i = ixBoundaryEnd-1; i >= ixBoundaryBegin; i--) {
			for (int j = m_meta.ixSecond-1; j >= m_meta.ixFirst; j--) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}

		} else {
			for (int k = 0; k < sRElements; k++) {
			for (int i = ixBoundaryEnd-1; i >= ixBoundaryBegin; i--) {
			for (int j = m_meta.ixFirst; j < m_meta.ixSecond; j++) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}
		}

	// Pack data to send bottomward
	} else if (m_key.dir == Direction_Bottom) {
		ixBoundaryBegin = m_meta.sHaloElements;
		ixBoundaryEnd   = 2 * m_meta.sHaloElements;

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixBoundaryEnd - ixBoundaryBegin)
			* (m_meta.ixSecond - m_meta.ixFirst);

		if (m_ixSendBuffer + nTotalValues > m_vecSendBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in SendBuffer for operation.");
		}

		// Pack the SendBuffer
		if (m_fReverseDirection) {
			for (int k = 0; k < sRElements; k++) {
			for (int j = ixBoundaryEnd-1; j >= ixBoundaryBegin; j--) {
			for (int i = m_meta.ixSecond-1; i >= m_meta.ixFirst; i--) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}

		} else {
			for (int k = 0; k < sRElements; k++) {
			for (int j = ixBoundaryEnd-1; j >= ixBoundaryBegin; j--) {
			for (int i = m_meta.ixFirst; i < m_meta.ixSecond; i++) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}
		}

	// Pack data to send toprightward 
	} else if (m_key.dir == Direction_TopRight) {
		ixABoundaryBegin = m_meta.ixFirst - m_meta.sHaloElements + 1;
		ixABoundaryEnd   = m_meta.ixFirst + 1;
		ixBBoundaryBegin = m_meta.ixSecond - m_meta.sHaloElements + 1;
		ixBBoundaryEnd   = m_meta.ixSecond + 1;
 
		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixABoundaryEnd - ixABoundaryBegin)
			* (ixBBoundaryEnd - ixBBoundaryBegin);

		if (m_ixSendBuffer + nTotalValues > m_vecSendBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in SendBuffer for operation.");
		}

		// Pack the SendBuffer
		if (m_fReverseDirection) {
			for (int k = 0; k < sRElements; k++) {
			for (int i = ixABoundaryBegin; i < ixABoundaryEnd; i++) {
			for (int j = ixBBoundaryBegin; j < ixBBoundaryEnd; j++) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}

		} else {
			for (int k = 0; k < sRElements; k++) {
			for (int j = ixBBoundaryBegin; j < ixBBoundaryEnd; j++) {
			for (int i = ixABoundaryBegin; i < ixABoundaryEnd; i++) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}
		}

	// Pack data to send topleftward 
	} else if (m_key.dir == Direction_TopLeft) {
		ixABoundaryBegin = m_meta.ixFirst;
		ixABoundaryEnd   = m_meta.ixFirst + m_meta.sHaloElements;
		ixBBoundaryBegin = m_meta.ixSecond - m_meta.sHaloElements + 1;
		ixBBoundaryEnd   = m_meta.ixSecond + 1;
 
		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixABoundaryEnd - ixABoundaryBegin)
			* (ixBBoundaryEnd - ixBBoundaryBegin);

		if (m_ixSendBuffer + nTotalValues > m_vecSendBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in SendBuffer for operation.");
		}

		// Pack the SendBuffer
		if (m_fReverseDirection) {
			for (int k = 0; k < sRElements; k++) {
			for (int i = ixABoundaryEnd-1; i >= ixABoundaryBegin; i--) {
			for (int j = ixBBoundaryBegin; j < ixBBoundaryEnd; j++) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}

		} else {
			for (int k = 0; k < sRElements; k++) {
			for (int j = ixBBoundaryBegin; j < ixBBoundaryEnd; j++) {
			for (int i = ixABoundaryEnd-1; i >= ixABoundaryBegin; i--) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}
		}

	// Pack data to send bottomleftward 
	} else if (m_key.dir == Direction_BottomLeft) {
		ixABoundaryBegin = m_meta.ixFirst;
		ixABoundaryEnd   = m_meta.ixFirst + m_meta.sHaloElements;
		ixBBoundaryBegin = m_meta.ixSecond;
		ixBBoundaryEnd   = m_meta.ixSecond + m_meta.sHaloElements;
 
		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixABoundaryEnd - ixABoundaryBegin)
			* (ixBBoundaryEnd - ixBBoundaryBegin);

		if (m_ixSendBuffer + nTotalValues > m_vecSendBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in SendBuffer for operation.");
		}

		// Pack the SendBuffer
		if (m_fReverseDirection) {
			for (int k = 0; k < sRElements; k++) {
			for (int i = ixABoundaryEnd-1; i >= ixABoundaryBegin; i--) {
			for (int j = ixBBoundaryEnd-1; j >= ixBBoundaryBegin; j--) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}

		} else {
			for (int k = 0; k < sRElements; k++) {
			for (int j = ixBBoundaryEnd-1; j >= ixBBoundaryBegin; j--) {
			for (int i = ixABoundaryEnd-1; i >= ixABoundaryBegin; i--) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}
		}

	// Pack data to send bottomrightward 
	} else if (m_key.dir == Direction_BottomRight) {
		ixABoundaryBegin = m_meta.ixFirst - m_meta.sHaloElements + 1;
		ixABoundaryEnd   = m_meta.ixFirst + 1;
		ixBBoundaryBegin = m_meta.ixSecond;
		ixBBoundaryEnd   = m_meta.ixSecond + m_meta.sHaloElements;
 
		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixABoundaryEnd - ixABoundaryBegin)
			* (ixBBoundaryEnd - ixBBoundaryBegin);

		if (m_ixSendBuffer + nTotalValues > m_vecSendBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in SendBuffer for operation.");
		}

		// Pack the SendBuffer
		if (m_fReverseDirection) {
			for (int k = 0; k < sRElements; k++) {
			for (int i = ixABoundaryBegin; i < ixABoundaryEnd; i++) {
			for (int j = ixBBoundaryEnd-1; j >= ixBBoundaryBegin; j--) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}

		} else {
			for (int k = 0; k < sRElements; k++) {
			for (int j = ixBBoundaryEnd-1; j >= ixBBoundaryBegin; j--) {
			for (int i = ixABoundaryBegin; i < ixABoundaryEnd; i++) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
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

void ExteriorNeighbor::Pack(
	const Grid & grid,
	const DataArray4D<double> & data
) {
	// Number of components in data
	size_t sComponents = data.GetSize(0);

	// 3D Grid Data
	DataArray3D<double> data3D(
		data.GetSize(1),
		data.GetSize(2),
		data.GetSize(3),
		data.GetDataType(),
		data.GetDataLocation(),
		false);

	// For state data exclude non-collacted data points
	if (data.GetDataType() == DataType_State) {
		// List of variable indices to send
		// - exclude variables which are not-collocated with this data structure
		for (int c = 0; c < sComponents; c++) {
			if (grid.GetVarLocation(c) != data.GetDataLocation()) {
				continue;
			}
			data3D.AttachTo(data[c]);
			Pack(grid, data3D);
			data3D.Detach();
		}

	// Send everything
	} else {
		for (int c = 0; c < sComponents; c++) {
			data3D.AttachTo(data[c]);
			Pack(grid, data3D);
			data3D.Detach();
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void ExteriorNeighbor::Send(
	const Grid & grid
) {

#ifdef USE_MPI
	// Send the data
	int iProcessor = grid.GetPatchProcessor(m_key.ixTargetPatch);

	// Tag
	if (m_key.ixTargetPatch >= (2 << 12)) {
		_EXCEPTIONT("Maximum neighbor index exceeded.");
	}

	int nTag = (m_key.ixSourcePatch << 16) + (m_key.ixTargetPatch << 4) + (int)(m_meta.dirOpposing);

	MPI_Isend(
		&(m_vecSendBuffer[0]),
		m_ixSendBuffer,
		MPI_DOUBLE,
		iProcessor,
		nTag,
		MPI_COMM_WORLD,
		&m_reqSend);
/*
	MPI_Send(
		&(m_vecSendBuffer[0]),
		m_ixSendBuffer,
		MPI_DOUBLE,
		iProcessor,
		nTag,
		MPI_COMM_WORLD);
*/
#endif
}

///////////////////////////////////////////////////////////////////////////////

void ExteriorNeighbor::Unpack(
	const Grid & grid,
	DataArray3D<double> & data
) {
/*
	// Check receive status
	int fRecvWaiting;
	MPI_Status status;
	MPI_Test(&m_reqRecv, &fRecvWaiting, &status);
	if (!fRecvWaiting) {
		_EXCEPTIONT("Receive buffer access error");
	}
	std::cout << "Test: " << status.MPI_TAG << std::endl;
*/
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
	if (m_key.dir == Direction_Right) {
		ixBoundaryBegin = sAElements - m_meta.sHaloElements;
		ixBoundaryEnd   = sAElements;

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixBoundaryEnd - ixBoundaryBegin)
			* (m_meta.ixSecond - m_meta.ixFirst);

		if (m_ixRecvBuffer + nTotalValues > m_vecRecvBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in RecvBuffer for operation.");
		}

		// Unpack data
		for (int k = 0; k < sRElements; k++) {
		for (int i = ixBoundaryEnd-1; i >= ixBoundaryBegin; i--) {
		for (int j = m_meta.ixFirst; j < m_meta.ixSecond; j++) {
			data[k][i][j] = m_vecRecvBuffer[m_ixRecvBuffer++];
		}
		}
		}

	// Unpack data from top
	} else if (m_key.dir == Direction_Top) {
		ixBoundaryBegin = sBElements - m_meta.sHaloElements;
		ixBoundaryEnd   = sBElements;

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixBoundaryEnd - ixBoundaryBegin)
			* (m_meta.ixSecond - m_meta.ixFirst);

		if (m_ixRecvBuffer + nTotalValues > m_vecRecvBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in RecvBuffer for operation.");
		}

		// Unpack data
		for (int k = 0; k < sRElements; k++) {
		for (int j = ixBoundaryEnd-1; j >= ixBoundaryBegin; j--) {
		for (int i = m_meta.ixFirst; i < m_meta.ixSecond; i++) {
			data[k][i][j] = m_vecRecvBuffer[m_ixRecvBuffer++];
		}
		}
		}

	// Unpack data from left
	} else if (m_key.dir == Direction_Left) {
		ixBoundaryBegin = 0;
		ixBoundaryEnd   = m_meta.sHaloElements;

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixBoundaryEnd - ixBoundaryBegin)
			* (m_meta.ixSecond - m_meta.ixFirst);

		if (m_ixRecvBuffer + nTotalValues > m_vecRecvBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in RecvBuffer for operation.");
		}

		// Unpack data
		for (int k = 0; k < sRElements; k++) {
		for (int i = ixBoundaryBegin; i < ixBoundaryEnd; i++) {
		for (int j = m_meta.ixFirst; j < m_meta.ixSecond; j++) {
			data[k][i][j] = m_vecRecvBuffer[m_ixRecvBuffer++];
		}
		}
		}

	// Unpack data from bottom
	} else if (m_key.dir == Direction_Bottom) {
		ixBoundaryBegin = 0;
		ixBoundaryEnd   = m_meta.sHaloElements;

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixBoundaryEnd - ixBoundaryBegin)
			* (m_meta.ixSecond - m_meta.ixFirst);

		if (m_ixRecvBuffer + nTotalValues > m_vecRecvBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in RecvBuffer for operation.");
		}

		// Unpack data
		for (int k = 0; k < sRElements; k++) {
		for (int j = ixBoundaryBegin; j < ixBoundaryEnd; j++) {
		for (int i = m_meta.ixFirst; i < m_meta.ixSecond; i++) {
			data[k][i][j] = m_vecRecvBuffer[m_ixRecvBuffer++];
		}
		}
		}

	// Unpack data from top-right
	} else if (m_key.dir == Direction_TopRight) {
		ixABoundaryBegin = m_meta.ixFirst + 1;
		ixABoundaryEnd   = m_meta.ixFirst + m_meta.sHaloElements + 1;
		ixBBoundaryBegin = m_meta.ixSecond + 1;
		ixBBoundaryEnd   = m_meta.ixSecond + m_meta.sHaloElements + 1;

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixABoundaryEnd - ixABoundaryBegin)
			* (ixBBoundaryEnd - ixBBoundaryBegin);

		// Unpack data
		for (int k = 0; k < sRElements; k++) {
		for (int j = ixBBoundaryEnd-1; j >= ixBBoundaryBegin; j--) {
		for (int i = ixABoundaryEnd-1; i >= ixABoundaryBegin; i--) {
			data[k][i][j] = m_vecRecvBuffer[m_ixRecvBuffer++];
		}
		}
		}

	// Unpack data from top-left
	} else if (m_key.dir == Direction_TopLeft) {
		ixABoundaryBegin = m_meta.ixFirst - m_meta.sHaloElements;
		ixABoundaryEnd   = m_meta.ixFirst;
		ixBBoundaryBegin = m_meta.ixSecond + 1;
		ixBBoundaryEnd   = m_meta.ixSecond + m_meta.sHaloElements + 1;

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixABoundaryEnd - ixABoundaryBegin)
			* (ixBBoundaryEnd - ixBBoundaryBegin);

		// Unpack data
		for (int k = 0; k < sRElements; k++) {
		for (int j = ixBBoundaryEnd-1; j >= ixBBoundaryBegin; j--) {
		for (int i = ixABoundaryBegin; i < ixABoundaryEnd; i++) {
			data[k][i][j] = m_vecRecvBuffer[m_ixRecvBuffer++];
		}
		}
		}

	// Unpack data from bottom-left
	} else if (m_key.dir == Direction_BottomLeft) {
		ixABoundaryBegin = m_meta.ixFirst - m_meta.sHaloElements;
		ixABoundaryEnd   = m_meta.ixFirst;
		ixBBoundaryBegin = m_meta.ixSecond - m_meta.sHaloElements;
		ixBBoundaryEnd   = m_meta.ixSecond;

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixABoundaryEnd - ixABoundaryBegin)
			* (ixBBoundaryEnd - ixBBoundaryBegin);

		// Unpack data
		for (int k = 0; k < sRElements; k++) {
		for (int j = ixBBoundaryBegin; j < ixBBoundaryEnd; j++) {
		for (int i = ixABoundaryBegin; i < ixABoundaryEnd; i++) {
			data[k][i][j] = m_vecRecvBuffer[m_ixRecvBuffer++];
		}
		}
		}

	// Unpack data from bottom-right
	} else if (m_key.dir == Direction_BottomRight) {
		ixABoundaryBegin = m_meta.ixFirst + 1;
		ixABoundaryEnd   = m_meta.ixFirst + m_meta.sHaloElements + 1;
		ixBBoundaryBegin = m_meta.ixSecond - m_meta.sHaloElements;
		ixBBoundaryEnd   = m_meta.ixSecond;

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixABoundaryEnd - ixABoundaryBegin)
			* (ixBBoundaryEnd - ixBBoundaryBegin);

		// Unpack data
		for (int k = 0; k < sRElements; k++) {
		for (int j = ixBBoundaryBegin; j < ixBBoundaryEnd; j++) {
		for (int i = ixABoundaryEnd-1; i >= ixABoundaryBegin; i--) {
			data[k][i][j] = m_vecRecvBuffer[m_ixRecvBuffer++];
		}
		}
		}

	// Invalid direction
	} else {
		_EXCEPTIONT("Invalid direction");
	}
}

///////////////////////////////////////////////////////////////////////////////

void ExteriorNeighbor::Unpack(
	const Grid & grid,
	DataArray4D<double> & data
) {
	// Number of components in data
	size_t sComponents = data.GetSize(0);

	// 3D Grid Data
	DataArray3D<double> data3D(
		data.GetSize(1),
		data.GetSize(2),
		data.GetSize(3),
		data.GetDataType(),
		data.GetDataLocation(),
		false);

	// List of variable indices to receive
	// - exclude variables which are not-collocated with this data structure
	if (data.GetDataType() == DataType_State) {
		for (int c = 0; c < sComponents; c++) {
			if (grid.GetVarLocation(c) != data.GetDataLocation()) {
				continue;
			}
			data3D.AttachTo(data[c]);;
			Unpack(grid, data3D);
			data3D.Detach();
		}

	// Unpack all variables
	} else {
		for (int c = 0; c < sComponents; c++) {
			data3D.AttachTo(data[c]);
			Unpack(grid, data3D);
			data3D.Detach();
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void ExteriorNeighbor::WaitSend() {

#ifdef USE_MPI
	MPI_Status status;
	MPI_Wait(&m_reqSend, &status);
#endif

}

///////////////////////////////////////////////////////////////////////////////
// Connectivity
///////////////////////////////////////////////////////////////////////////////

Connectivity::~Connectivity() {
/*
	for (int n = 0; n < m_vecExteriorNeighbors.size(); n++) {
		delete m_vecExteriorNeighbors[n];
	}
*/
/*
	for (int n = 0; n < m_vecNestedNeighbors; n++) {
		delete m_vecNestedNeighbors[n];
	}
*/
}

///////////////////////////////////////////////////////////////////////////////

void Connectivity::BuildFluxConnectivity() {
	_EXCEPTIONT("Disabled; should be called once all ExteriorNeighbors are set");
/*
	// PatchBox
	const PatchBox & box = m_patch.GetPatchBox();

	// Allocate space for edges
	m_vecExteriorEdge[0].resize(box.GetBTotalWidth(), NULL);
	m_vecExteriorEdge[1].resize(box.GetATotalWidth(), NULL);
	m_vecExteriorEdge[2].resize(box.GetBTotalWidth(), NULL);
	m_vecExteriorEdge[3].resize(box.GetATotalWidth(), NULL);
	m_vecExteriorEdge[4].resize(1, NULL);
	m_vecExteriorEdge[5].resize(1, NULL);
	m_vecExteriorEdge[6].resize(1, NULL);
	m_vecExteriorEdge[7].resize(1, NULL);

	// Loop over all exterior edges, hook up pointers to ExteriorNeighbors
	for (int n = 0; n < m_vecExteriorNeighbors.size(); n++) {

		ExteriorNeighbor * pNeighbor = m_vecExteriorNeighbors[n];

		int iDir = static_cast<int>(pNeighbor->m_key.dir);

		// Right or Left side of this PatchBox
		if ((pNeighbor->m_key.dir == Direction_Right) || 
		    (pNeighbor->m_key.dir == Direction_Left)
		) {
			int j = pNeighbor->m_meta.ixFirst;
			for (; j < pNeighbor->m_meta.ixSecond; j++) {
				if ((j < box.GetBInteriorBegin()) ||
					(j >= box.GetBInteriorEnd())
				) {
					_EXCEPTIONT("Edge index out of range");
				}

				m_vecExteriorEdge[iDir][j] = pNeighbor;
			}

		// Top or Bottom side of this PatchBox
		} else if (
		    (pNeighbor->m_key.dir == Direction_Top) ||
		    (pNeighbor->m_key.dir == Direction_Bottom)
		) {
			int i = pNeighbor->m_meta.ixFirst;
			for (; i < pNeighbor->m_meta.ixSecond; i++) {
				if ((i < box.GetAInteriorBegin()) ||
					(i >= box.GetAInteriorEnd())
				) {
					_EXCEPTIONT("Edge index out of range");
				}

				m_vecExteriorEdge[iDir][i] = pNeighbor;
			}

		// Corners of this PatchBox
		} else if (
		    (pNeighbor->m_key.dir == Direction_TopRight) ||
		    (pNeighbor->m_key.dir == Direction_TopLeft) ||
		    (pNeighbor->m_key.dir == Direction_BottomRight) ||
		    (pNeighbor->m_key.dir == Direction_BottomLeft)
		) {
			if (m_vecExteriorEdge[iDir][0] != NULL) {
				_EXCEPTION1("Corner patch %i already set", iDir);
			}
			m_vecExteriorEdge[iDir][0] = pNeighbor;

		} else {
			_EXCEPTIONT("Invalid direction");
		}
	}
*/
}

///////////////////////////////////////////////////////////////////////////////

void Connectivity::PrepareExchange() {

	// Prepare for asynchronous receives from each neighbor
	for (int m = 0; m < m_vecExteriorNeighbors.size(); m++) {
		m_vecExteriorNeighbors[m]->PrepareExchange(m_patch.GetGrid());
	}
}

///////////////////////////////////////////////////////////////////////////////

void Connectivity::Pack(
	const DataArray3D<double> & data
) {
	// Pack data to exterior neighbors
	for (int m = 0; m < m_vecExteriorNeighbors.size(); m++) {
		m_vecExteriorNeighbors[m]->Pack(m_patch.GetGrid(), data);
	}
}

///////////////////////////////////////////////////////////////////////////////

void Connectivity::Pack(
	const DataArray4D<double> & data
) {
	// Pack data to exterior neighbors
	for (int m = 0; m < m_vecExteriorNeighbors.size(); m++) {
		m_vecExteriorNeighbors[m]->Pack(m_patch.GetGrid(), data);
	}
}

///////////////////////////////////////////////////////////////////////////////

void Connectivity::Send() {

	// Send data to exterior neighbors
	for (int m = 0; m < m_vecExteriorNeighbors.size(); m++) {
		m_vecExteriorNeighbors[m]->Send(m_patch.GetGrid());
	}
}

///////////////////////////////////////////////////////////////////////////////

void Connectivity::SendBuffers() {

	// Send data to exterior neighbors
	for (int m = 0; m < m_vecExteriorNeighbors.size(); m++) {
		m_vecExteriorNeighbors[m]->ResetSendBufferSize();
		m_vecExteriorNeighbors[m]->Send(m_patch.GetGrid());
	}
}

///////////////////////////////////////////////////////////////////////////////

Neighbor * Connectivity::WaitReceive() {

	// Receive data from exterior neighbors
	int nRecvMessageCount = 0;

	while (nRecvMessageCount != m_vecExteriorNeighbors.size()) {

		nRecvMessageCount = 0;
		for (int m = 0; m < m_vecExteriorNeighbors.size(); m++) {
			if (m_vecExteriorNeighbors[m]->IsComplete()) {
				nRecvMessageCount++;
				continue;
			}
			if (m_vecExteriorNeighbors[m]->CheckReceive()) {
				m_vecExteriorNeighbors[m]->SetComplete();
				return m_vecExteriorNeighbors[m];
			}
		}
	}

	return NULL;
}

///////////////////////////////////////////////////////////////////////////////

void Connectivity::WaitSend() {

	// Wait for all asynchronous send requests to complete
	for (int m = 0; m < m_vecExteriorNeighbors.size(); m++) {
		m_vecExteriorNeighbors[m]->WaitSend();
	}
}

///////////////////////////////////////////////////////////////////////////////

int Connectivity::GetExpectedMessageCount() const {
	return (int)(m_vecExteriorNeighbors.size());// + m_vecNestedNeighbors.size());
}

///////////////////////////////////////////////////////////////////////////////

