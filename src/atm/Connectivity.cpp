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

void Neighbor::InitializeBuffers(
	size_t sRElements,
	size_t sHaloElements,
	size_t sComponents
) {

	// Store buffer size
	m_sMaxRElements = sRElements + 1;
	m_sHaloElements = sHaloElements;
	m_sComponents = sComponents;

	// Resize buffers
	m_vecSendBuffer.Allocate(
		  m_sBoundarySize
		* m_sMaxRElements
		* m_sHaloElements
		* m_sComponents);

	m_vecRecvBuffer.Allocate(
		  m_sBoundarySize
		* m_sMaxRElements
		* m_sHaloElements
		* m_sComponents);

}

///////////////////////////////////////////////////////////////////////////////

bool Neighbor::CheckReceive() {

	// Check if message already received and processed
	if (m_fComplete) {
		return false;
	}

	// Test receive request
	int fRecvWaiting;
	MPI_Status status;
	MPI_Test(&m_reqRecv, &fRecvWaiting, &status);
	if (!fRecvWaiting) {
		return false;
	}

	// A message has been received
	return true;
}

///////////////////////////////////////////////////////////////////////////////
// ExteriorNeighbor
///////////////////////////////////////////////////////////////////////////////

void ExteriorNeighbor::PrepareExchange() {

	// Call up the stack
	Neighbor::PrepareExchange();

	// Information for receive
	const GridPatch & patch = m_pConnect->GetGridPatch();

	const Grid & grid = patch.GetGrid();

	int ixPatch = patch.GetPatchIndex();

	int iProcessor = grid.GetPatch(m_ixNeighbor)->GetProcessor();

	// Tag of the received message
	int nTag = (m_ixNeighbor << 16) + (ixPatch << 4) + (int)(m_dir);

	// Prepare an asynchronous receive
	MPI_Irecv(
		&(m_vecRecvBuffer[0]),
		m_vecRecvBuffer.GetRows(),
		MPI_DOUBLE,
		iProcessor,
		nTag,
		MPI_COMM_WORLD,
		&m_reqRecv);
}

///////////////////////////////////////////////////////////////////////////////

void ExteriorNeighbor::Pack(
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

	// Model grid
	const Grid & grid = m_pConnect->GetGridPatch().GetGrid();

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

		// Pack the SendBuffer
		if (m_fReverseDirection) {
			for (int k = 0; k < sRElements; k++) {
			for (int i = ixBoundaryBegin; i < ixBoundaryEnd; i++) {
			for (int j = m_ixSecond-1; j >= m_ixFirst; j--) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}

		} else {
			for (int k = 0; k < sRElements; k++) {
			for (int i = ixBoundaryBegin; i < ixBoundaryEnd; i++) {
			for (int j = m_ixFirst; j < m_ixSecond; j++) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];	
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

		if (m_ixSendBuffer + nTotalValues > m_vecSendBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in SendBuffer for operation.");
		}

		// Pack the SendBuffer
		if (m_fReverseDirection) {
			for (int k = 0; k < sRElements; k++) {
			for (int j = ixBoundaryBegin; j < ixBoundaryEnd; j++) {
			for (int i = m_ixSecond-1; i >= m_ixFirst; i--) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}

		} else {
			for (int k = 0; k < sRElements; k++) {
			for (int j = ixBoundaryBegin; j < ixBoundaryEnd; j++) {
			for (int i = m_ixFirst; i < m_ixSecond; i++) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
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

		if (m_ixSendBuffer + nTotalValues > m_vecSendBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in SendBuffer for operation.");
		}

		// Pack the SendBuffer
		if (m_fReverseDirection) {
			for (int k = 0; k < sRElements; k++) {
			for (int i = ixBoundaryEnd-1; i >= ixBoundaryBegin; i--) {
			for (int j = m_ixSecond-1; j >= m_ixFirst; j--) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}

		} else {
			for (int k = 0; k < sRElements; k++) {
			for (int i = ixBoundaryEnd-1; i >= ixBoundaryBegin; i--) {
			for (int j = m_ixFirst; j < m_ixSecond; j++) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
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

		if (m_ixSendBuffer + nTotalValues > m_vecSendBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in SendBuffer for operation.");
		}

		// Pack the SendBuffer
		if (m_fReverseDirection) {
			for (int k = 0; k < sRElements; k++) {
			for (int j = ixBoundaryEnd-1; j >= ixBoundaryBegin; j--) {
			for (int i = m_ixSecond-1; i >= m_ixFirst; i--) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}

		} else {
			for (int k = 0; k < sRElements; k++) {
			for (int j = ixBoundaryEnd-1; j >= ixBoundaryBegin; j--) {
			for (int i = m_ixFirst; i < m_ixSecond; i++) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
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
	const DataArray4D<double> & data
) {
	// Model grid
	const Grid & grid = m_pConnect->GetGridPatch().GetGrid();

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
			data3D.AttachTo(data[c][0][0]);
			Pack(data3D);
			data3D.Detach();
		}

	// Send everything
	} else {
		for (int c = 0; c < sComponents; c++) {
			data3D.AttachTo(data[c][0][0]);
			Pack(data3D);
			data3D.Detach();
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void ExteriorNeighbor::Send() {

	// Send the data
	const GridPatch & patch = m_pConnect->GetGridPatch();

	const Grid & grid = patch.GetGrid();

	int ixPatch = patch.GetPatchIndex();

	int iProcessor = grid.GetPatch(m_ixNeighbor)->GetProcessor();

	// Tag
	if (m_ixNeighbor >= (2 << 12)) {
		_EXCEPTIONT("Maximum neighbor index exceeded.");
	}

	int nTag = (ixPatch << 16) + (m_ixNeighbor << 4) + (int)(m_dirOpposite);

#pragma message "Move MPI_TAG processing to Connectivity"

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
}

///////////////////////////////////////////////////////////////////////////////

void ExteriorNeighbor::Unpack(
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

	// Model grid
	const Grid & grid = m_pConnect->GetGridPatch().GetGrid();

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

		if (m_ixRecvBuffer + nTotalValues > m_vecRecvBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in RecvBuffer for operation.");
		}

		// Unpack data
		for (int k = 0; k < sRElements; k++) {
		for (int i = ixBoundaryEnd-1; i >= ixBoundaryBegin; i--) {
		for (int j = m_ixFirst; j < m_ixSecond; j++) {
			data[k][i][j] = m_vecRecvBuffer[m_ixRecvBuffer++];
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

		if (m_ixRecvBuffer + nTotalValues > m_vecRecvBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in RecvBuffer for operation.");
		}

		// Unpack data
		for (int k = 0; k < sRElements; k++) {
		for (int j = ixBoundaryEnd-1; j >= ixBoundaryBegin; j--) {
		for (int i = m_ixFirst; i < m_ixSecond; i++) {
			data[k][i][j] = m_vecRecvBuffer[m_ixRecvBuffer++];
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

		if (m_ixRecvBuffer + nTotalValues > m_vecRecvBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in RecvBuffer for operation.");
		}

		// Unpack data
		for (int k = 0; k < sRElements; k++) {
		for (int i = ixBoundaryBegin; i < ixBoundaryEnd; i++) {
		for (int j = m_ixFirst; j < m_ixSecond; j++) {
			data[k][i][j] = m_vecRecvBuffer[m_ixRecvBuffer++];
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

		if (m_ixRecvBuffer + nTotalValues > m_vecRecvBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in RecvBuffer for operation.");
		}

		// Unpack data
		for (int k = 0; k < sRElements; k++) {
		for (int j = ixBoundaryBegin; j < ixBoundaryEnd; j++) {
		for (int i = m_ixFirst; i < m_ixSecond; i++) {
			data[k][i][j] = m_vecRecvBuffer[m_ixRecvBuffer++];
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

		// Unpack data
		for (int k = 0; k < sRElements; k++) {
		for (int j = ixBBoundaryEnd-1; j >= ixBBoundaryBegin; j--) {
		for (int i = ixABoundaryEnd-1; i >= ixABoundaryBegin; i--) {
			data[k][i][j] = m_vecRecvBuffer[m_ixRecvBuffer++];
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

		// Unpack data
		for (int k = 0; k < sRElements; k++) {
		for (int j = ixBBoundaryEnd-1; j >= ixBBoundaryBegin; j--) {
		for (int i = ixABoundaryBegin; i < ixABoundaryEnd; i++) {
			data[k][i][j] = m_vecRecvBuffer[m_ixRecvBuffer++];
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

		// Unpack data
		for (int k = 0; k < sRElements; k++) {
		for (int j = ixBBoundaryBegin; j < ixBBoundaryEnd; j++) {
		for (int i = ixABoundaryBegin; i < ixABoundaryEnd; i++) {
			data[k][i][j] = m_vecRecvBuffer[m_ixRecvBuffer++];
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
	DataArray4D<double> & data
) {
	// Model grid
	const Grid & grid = m_pConnect->GetGridPatch().GetGrid();

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
			data3D.AttachTo(data[c][0][0]);;
			Unpack(data3D);
			data3D.Detach();
		}

	// Unpack all variables
	} else {
		for (int c = 0; c < sComponents; c++) {
			data3D.AttachTo(data[c][0][0]);;
			Unpack(data3D);
			data3D.Detach();
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void ExteriorNeighbor::WaitSend() {

	MPI_Status status;
	MPI_Wait(&m_reqSend, &status);

}

///////////////////////////////////////////////////////////////////////////////
// Connectivity
///////////////////////////////////////////////////////////////////////////////

Connectivity::~Connectivity() {
	for (int n = 0; n < m_vecExteriorNeighbors.size(); n++) {
		delete m_vecExteriorNeighbors[n];
	}
/*
	for (int n = 0; n < m_vecNestedNeighbors; n++) {
		delete m_vecNestedNeighbors[n];
	}
*/
}

///////////////////////////////////////////////////////////////////////////////

void Connectivity::BuildFluxConnectivity() {

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

		int iDir = static_cast<int>(pNeighbor->m_dir);

		// Right or Left side of this PatchBox
		if ((pNeighbor->m_dir == Direction_Right) || 
		    (pNeighbor->m_dir == Direction_Left)
		) {
			int j = pNeighbor->m_ixFirst;
			for (; j < pNeighbor->m_ixSecond; j++) {
				if ((j < box.GetBInteriorBegin()) ||
					(j >= box.GetBInteriorEnd())
				) {
					_EXCEPTIONT("Edge index out of range");
				}

				m_vecExteriorEdge[iDir][j] = pNeighbor;
			}

		// Top or Bottom side of this PatchBox
		} else if (
		    (pNeighbor->m_dir == Direction_Top) ||
		    (pNeighbor->m_dir == Direction_Bottom)
		) {
			int i = pNeighbor->m_ixFirst;
			for (; i < pNeighbor->m_ixSecond; i++) {
				if ((i < box.GetAInteriorBegin()) ||
					(i >= box.GetAInteriorEnd())
				) {
					_EXCEPTIONT("Edge index out of range");
				}

				m_vecExteriorEdge[iDir][i] = pNeighbor;
			}

		// Corners of this PatchBox
		} else if (
		    (pNeighbor->m_dir == Direction_TopRight) ||
		    (pNeighbor->m_dir == Direction_TopLeft) ||
		    (pNeighbor->m_dir == Direction_BottomRight) ||
		    (pNeighbor->m_dir == Direction_BottomLeft)
		) {
			if (m_vecExteriorEdge[iDir][0] != NULL) {
				_EXCEPTION1("Corner patch %i already set", iDir);
			}
			m_vecExteriorEdge[iDir][0] = pNeighbor;

		} else {
			_EXCEPTIONT("Invalid direction");
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void Connectivity::PrepareExchange() {

	// Prepare for asynchronous receives from each neighbor
	for (int m = 0; m < m_vecExteriorNeighbors.size(); m++) {
		m_vecExteriorNeighbors[m]->PrepareExchange();
	}
}

///////////////////////////////////////////////////////////////////////////////

void Connectivity::Pack(
	const DataArray3D<double> & data
) {
	// Pack data to exterior neighbors
	for (int m = 0; m < m_vecExteriorNeighbors.size(); m++) {
		m_vecExteriorNeighbors[m]->Pack(data);
	}
}

///////////////////////////////////////////////////////////////////////////////

void Connectivity::Pack(
	const DataArray4D<double> & data
) {
	// Pack data to exterior neighbors
	for (int m = 0; m < m_vecExteriorNeighbors.size(); m++) {
		m_vecExteriorNeighbors[m]->Pack(data);
	}
}

///////////////////////////////////////////////////////////////////////////////

void Connectivity::Send() {

	// Send data to exterior neighbors
	for (int m = 0; m < m_vecExteriorNeighbors.size(); m++) {
		m_vecExteriorNeighbors[m]->Send();
	}
}

///////////////////////////////////////////////////////////////////////////////

void Connectivity::SendBuffers() {

	// Send data to exterior neighbors
	for (int m = 0; m < m_vecExteriorNeighbors.size(); m++) {
		m_vecExteriorNeighbors[m]->ResetSendBufferSize();
		m_vecExteriorNeighbors[m]->Send();
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

