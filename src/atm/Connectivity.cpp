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

///////////////////////////////////////////////////////////////////////////////
// Neighbor
///////////////////////////////////////////////////////////////////////////////

void Neighbor::InitializeBuffers(
	int nRElements,
	int nHaloElements,
	int nVariables
) {
	m_vecSendBuffer.Initialize(
		m_nBoundarySize * (nRElements+1) * nHaloElements * nVariables);
	m_vecRecvBuffer.Initialize(
		m_nBoundarySize * (nRElements+1) * nHaloElements * nVariables);
}

///////////////////////////////////////////////////////////////////////////////

bool Neighbor::CheckReceive() const {

	// Check if message already received and processed
	if (m_fComplete) {
		return false;
	}

	// Check receive status
	int fRecvWaiting;
	MPI_Status status;
	MPI_Request_get_status(m_reqRecv, &fRecvWaiting, &status);
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
	const GridData3D & data
) {
	// Check matrix bounds
	if (((m_dir == Direction_Right) || (m_dir == Direction_Left)) &&
		(m_ixNodeEnd > data.GetBElements())
	) {
		_EXCEPTIONT("GridData / ExteriorNeighbor inconsistency.");
	}
	if (((m_dir == Direction_Top) || (m_dir == Direction_Bottom)) &&
		(m_ixNodeEnd > data.GetAElements())
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
		ixBoundaryBegin = data.GetAElements() - 2 * data.GetHaloElements();
		ixBoundaryEnd   = data.GetAElements() - data.GetHaloElements();

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  data.GetRElements()
			* (ixBoundaryEnd - ixBoundaryBegin)
			* (m_ixNodeEnd - m_ixNodeBegin);

		if (m_ixSendBuffer + nTotalValues > m_vecSendBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in SendBuffer for operation.");
		}

		// Pack the SendBuffer
		if (m_fReverseDirection) {
			for (int k = 0; k < data.GetRElements(); k++) {
			for (int i = ixBoundaryBegin; i < ixBoundaryEnd; i++) {
			for (int j = m_ixNodeEnd-1; j >= m_ixNodeBegin; j--) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}

		} else {
			for (int k = 0; k < data.GetRElements(); k++) {
			for (int i = ixBoundaryBegin; i < ixBoundaryEnd; i++) {
			for (int j = m_ixNodeBegin; j < m_ixNodeEnd; j++) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];	
			}
			}
			}
		}
	
	// Pack data to send topward
	} else if (m_dir == Direction_Top) {
		ixBoundaryBegin = data.GetBElements() - 2 * data.GetHaloElements();
		ixBoundaryEnd   = data.GetBElements() - data.GetHaloElements();

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  data.GetRElements()
			* (ixBoundaryEnd - ixBoundaryBegin)
			* (m_ixNodeEnd - m_ixNodeBegin);

		if (m_ixSendBuffer + nTotalValues > m_vecSendBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in SendBuffer for operation.");
		}

		// Pack the SendBuffer
		if (m_fReverseDirection) {
			for (int k = 0; k < data.GetRElements(); k++) {
			for (int j = ixBoundaryBegin; j < ixBoundaryEnd; j++) {
			for (int i = m_ixNodeEnd-1; i >= m_ixNodeBegin; i--) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}

		} else {
			for (int k = 0; k < data.GetRElements(); k++) {
			for (int j = ixBoundaryBegin; j < ixBoundaryEnd; j++) {
			for (int i = m_ixNodeBegin; i < m_ixNodeEnd; i++) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}
		}

	// Pack data to send left
	} else if (m_dir == Direction_Left) {
		ixBoundaryBegin = data.GetHaloElements();
		ixBoundaryEnd   = 2 * data.GetHaloElements();

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  data.GetRElements()
			* (ixBoundaryEnd - ixBoundaryBegin)
			* (m_ixNodeEnd - m_ixNodeBegin);

		if (m_ixSendBuffer + nTotalValues > m_vecSendBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in SendBuffer for operation.");
		}

		// Pack the SendBuffer
		if (m_fReverseDirection) {
			for (int k = 0; k < data.GetRElements(); k++) {
			for (int i = ixBoundaryEnd-1; i >= ixBoundaryBegin; i--) {
			for (int j = m_ixNodeEnd-1; j >= m_ixNodeBegin; j--) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}

		} else {
			for (int k = 0; k < data.GetRElements(); k++) {
			for (int i = ixBoundaryEnd-1; i >= ixBoundaryBegin; i--) {
			for (int j = m_ixNodeBegin; j < m_ixNodeEnd; j++) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}
		}

	// Pack data to send bottomward
	} else if (m_dir == Direction_Bottom) {
		ixBoundaryBegin = data.GetHaloElements();
		ixBoundaryEnd   = 2 * data.GetHaloElements();

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  data.GetRElements()
			* (ixBoundaryEnd - ixBoundaryBegin)
			* (m_ixNodeEnd - m_ixNodeBegin);

		if (m_ixSendBuffer + nTotalValues > m_vecSendBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in SendBuffer for operation.");
		}

		// Pack the SendBuffer
		if (m_fReverseDirection) {
			for (int k = 0; k < data.GetRElements(); k++) {
			for (int j = ixBoundaryEnd-1; j >= ixBoundaryBegin; j--) {
			for (int i = m_ixNodeEnd-1; i >= m_ixNodeBegin; i--) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}

		} else {
			for (int k = 0; k < data.GetRElements(); k++) {
			for (int j = ixBoundaryEnd-1; j >= ixBoundaryBegin; j--) {
			for (int i = m_ixNodeBegin; i < m_ixNodeEnd; i++) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}
		}

	// Pack data to send toprightward 
	} else if (m_dir == Direction_TopRight) {
		ixABoundaryBegin = data.GetAElements() - 2 * data.GetHaloElements();
		ixABoundaryEnd   = data.GetAElements() - data.GetHaloElements();
		ixBBoundaryBegin = data.GetBElements() - 2 * data.GetHaloElements();
		ixBBoundaryEnd   = data.GetBElements() - data.GetHaloElements();
 
		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  data.GetRElements()
			* (ixABoundaryEnd - ixABoundaryBegin)
			* (ixBBoundaryEnd - ixBBoundaryBegin);

		if (m_ixSendBuffer + nTotalValues > m_vecSendBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in SendBuffer for operation.");
		}

		// Pack the SendBuffer
		if (m_fReverseDirection) {
			for (int k = 0; k < data.GetRElements(); k++) {
			for (int i = ixABoundaryBegin; i < ixABoundaryEnd; i++) {
			for (int j = ixBBoundaryBegin; j < ixBBoundaryEnd; j++) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}

		} else {
			for (int k = 0; k < data.GetRElements(); k++) {
			for (int j = ixBBoundaryBegin; j < ixBBoundaryEnd; j++) {
			for (int i = ixABoundaryBegin; i < ixABoundaryEnd; i++) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}
		}

	// Pack data to send topleftward 
	} else if (m_dir == Direction_TopLeft) {
		ixABoundaryBegin = data.GetHaloElements();
		ixABoundaryEnd   = 2 * data.GetHaloElements();
		ixBBoundaryBegin = data.GetBElements() - 2 * data.GetHaloElements();
		ixBBoundaryEnd   = data.GetBElements() - data.GetHaloElements();
 
		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  data.GetRElements()
			* (ixABoundaryEnd - ixABoundaryBegin)
			* (ixBBoundaryEnd - ixBBoundaryBegin);

		if (m_ixSendBuffer + nTotalValues > m_vecSendBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in SendBuffer for operation.");
		}

		// Pack the SendBuffer
		if (m_fReverseDirection) {
			for (int k = 0; k < data.GetRElements(); k++) {
			for (int i = ixABoundaryEnd-1; i >= ixABoundaryBegin; i--) {
			for (int j = ixBBoundaryBegin; j < ixBBoundaryEnd; j++) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}

		} else {
			for (int k = 0; k < data.GetRElements(); k++) {
			for (int j = ixBBoundaryBegin; j < ixBBoundaryEnd; j++) {
			for (int i = ixABoundaryEnd-1; i >= ixABoundaryBegin; i--) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}
		}

	// Pack data to send bottomleftward 
	} else if (m_dir == Direction_BottomLeft) {
		ixABoundaryBegin = data.GetHaloElements();
		ixABoundaryEnd   = 2 * data.GetHaloElements();
		ixBBoundaryBegin = data.GetHaloElements();
		ixBBoundaryEnd   = 2 * data.GetHaloElements();
 
		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  data.GetRElements()
			* (ixABoundaryEnd - ixABoundaryBegin)
			* (ixBBoundaryEnd - ixBBoundaryBegin);

		if (m_ixSendBuffer + nTotalValues > m_vecSendBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in SendBuffer for operation.");
		}

		// Pack the SendBuffer
		if (m_fReverseDirection) {
			for (int k = 0; k < data.GetRElements(); k++) {
			for (int i = ixABoundaryEnd-1; i >= ixABoundaryBegin; i--) {
			for (int j = ixBBoundaryEnd-1; j >= ixBBoundaryBegin; j--) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}

		} else {
			for (int k = 0; k < data.GetRElements(); k++) {
			for (int j = ixBBoundaryEnd-1; j >= ixBBoundaryBegin; j--) {
			for (int i = ixABoundaryEnd-1; i >= ixABoundaryBegin; i--) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}
		}

	// Pack data to send bottomrightward 
	} else if (m_dir == Direction_BottomRight) {
		ixABoundaryBegin = data.GetAElements() - 2 * data.GetHaloElements();
		ixABoundaryEnd   = data.GetAElements() - data.GetHaloElements();
		ixBBoundaryBegin = data.GetHaloElements();
		ixBBoundaryEnd   = 2 * data.GetHaloElements();
 
		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  data.GetRElements()
			* (ixABoundaryEnd - ixABoundaryBegin)
			* (ixBBoundaryEnd - ixBBoundaryBegin);

		if (m_ixSendBuffer + nTotalValues > m_vecSendBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in SendBuffer for operation.");
		}

		// Pack the SendBuffer
		if (m_fReverseDirection) {
			for (int k = 0; k < data.GetRElements(); k++) {
			for (int i = ixABoundaryBegin; i < ixABoundaryEnd; i++) {
			for (int j = ixBBoundaryEnd-1; j >= ixBBoundaryBegin; j--) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}

		} else {
			for (int k = 0; k < data.GetRElements(); k++) {
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
	const GridData4D & data
) {
	// Model grid
	const Grid & grid = m_pConnect->GetGridPatch().GetGrid();

	// 3D Grid Data
	GridData3D data3D;

	// For state data exclude non-collacted data points
	if (data.GetDataType() == DataType_State) {
		// List of variable indices to send
		// - exclude variables which are not-collocated with this data structure
		for (int c = 0; c < data.GetComponents(); c++) {
			if (grid.GetVarLocation(c) != data.GetDataLocation()) {
				continue;
			}
			data.GetAsGridData3D(c, data3D);
			Pack(data3D);
		}

	// Send everything
	} else {
		for (int c = 0; c < data.GetComponents(); c++) {
			data.GetAsGridData3D(c, data3D);
			Pack(data3D);
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
/*
	if (grid.GetPatch(m_ixNeighbor)->GetPatchIndex() == 0) {
		std::cout << ixPatch << " sending to patch 0" << std::endl;
	}
*/
	// Tag
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
}

///////////////////////////////////////////////////////////////////////////////

void ExteriorNeighbor::Unpack(
	GridData3D & data
) {
	// Check receive status
	int fRecvWaiting;
	MPI_Status status;
	MPI_Request_get_status(m_reqRecv, &fRecvWaiting, &status);
	if (!fRecvWaiting) {
		_EXCEPTIONT("Receive buffer access error");
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

	// Unpack data from right
	if (m_dir == Direction_Right) {
		ixBoundaryBegin = data.GetAElements() - data.GetHaloElements();
		ixBoundaryEnd   = data.GetAElements();

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  data.GetRElements()
			* (ixBoundaryEnd - ixBoundaryBegin)
			* (m_ixNodeEnd - m_ixNodeBegin);

		if (m_ixRecvBuffer + nTotalValues > m_vecRecvBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in RecvBuffer for operation.");
		}

		// Unpack data
		for (int k = 0; k < data.GetRElements(); k++) {
		for (int i = ixBoundaryEnd-1; i >= ixBoundaryBegin; i--) {
		for (int j = m_ixNodeBegin; j < m_ixNodeEnd; j++) {
			data[k][i][j] = m_vecRecvBuffer[m_ixRecvBuffer++];
		}
		}
		}

	// Unpack data from top
	} else if (m_dir == Direction_Top) {
		ixBoundaryBegin = data.GetBElements() - data.GetHaloElements();
		ixBoundaryEnd   = data.GetBElements();

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  data.GetRElements()
			* (ixBoundaryEnd - ixBoundaryBegin)
			* (m_ixNodeEnd - m_ixNodeBegin);

		if (m_ixRecvBuffer + nTotalValues > m_vecRecvBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in RecvBuffer for operation.");
		}

		// Unpack data
		for (int k = 0; k < data.GetRElements(); k++) {
		for (int j = ixBoundaryEnd-1; j >= ixBoundaryBegin; j--) {
		for (int i = m_ixNodeBegin; i < m_ixNodeEnd; i++) {
			data[k][i][j] = m_vecRecvBuffer[m_ixRecvBuffer++];
		}
		}
		}

	// Unpack data from left
	} else if (m_dir == Direction_Left) {
		ixBoundaryBegin = 0;
		ixBoundaryEnd   = data.GetHaloElements();

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  data.GetRElements()
			* (ixBoundaryEnd - ixBoundaryBegin)
			* (m_ixNodeEnd - m_ixNodeBegin);

		if (m_ixRecvBuffer + nTotalValues > m_vecRecvBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in RecvBuffer for operation.");
		}

		// Unpack data
		for (int k = 0; k < data.GetRElements(); k++) {
		for (int i = ixBoundaryBegin; i < ixBoundaryEnd; i++) {
		for (int j = m_ixNodeBegin; j < m_ixNodeEnd; j++) {
			data[k][i][j] = m_vecRecvBuffer[m_ixRecvBuffer++];
		}
		}
		}

	// Unpack data from bottom
	} else if (m_dir == Direction_Bottom) {
		ixBoundaryBegin = 0;
		ixBoundaryEnd   = data.GetHaloElements();

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  data.GetRElements()
			* (ixBoundaryEnd - ixBoundaryBegin)
			* (m_ixNodeEnd - m_ixNodeBegin);

		if (m_ixRecvBuffer + nTotalValues > m_vecRecvBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in RecvBuffer for operation.");
		}

		// Unpack data
		for (int k = 0; k < data.GetRElements(); k++) {
		for (int j = ixBoundaryBegin; j < ixBoundaryEnd; j++) {
		for (int i = m_ixNodeBegin; i < m_ixNodeEnd; i++) {
			data[k][i][j] = m_vecRecvBuffer[m_ixRecvBuffer++];
		}
		}
		}

	// Unpack data from top-right
	} else if (m_dir == Direction_TopRight) {
		ixABoundaryBegin = data.GetAElements() - data.GetHaloElements();
		ixABoundaryEnd   = data.GetAElements();
		ixBBoundaryBegin = data.GetBElements() - data.GetHaloElements();
		ixBBoundaryEnd   = data.GetBElements();

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  data.GetRElements()
			* (ixABoundaryEnd - ixABoundaryBegin)
			* (ixBBoundaryEnd - ixBBoundaryBegin);

		// Unpack data
		for (int k = 0; k < data.GetRElements(); k++) {
		for (int j = ixBBoundaryEnd-1; j >= ixBBoundaryBegin; j--) {
		for (int i = ixABoundaryEnd-1; i >= ixABoundaryBegin; i--) {
			data[k][i][j] = m_vecRecvBuffer[m_ixRecvBuffer++];
		}
		}
		}

	// Unpack data from top-left
	} else if (m_dir == Direction_TopLeft) {
		ixABoundaryBegin = 0;
		ixABoundaryEnd   = data.GetHaloElements();
		ixBBoundaryBegin = data.GetBElements() - data.GetHaloElements();
		ixBBoundaryEnd   = data.GetBElements();

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  data.GetRElements()
			* (ixABoundaryEnd - ixABoundaryBegin)
			* (ixBBoundaryEnd - ixBBoundaryBegin);

		// Unpack data
		for (int k = 0; k < data.GetRElements(); k++) {
		for (int j = ixBBoundaryEnd-1; j >= ixBBoundaryBegin; j--) {
		for (int i = ixABoundaryBegin; i < ixABoundaryEnd; i++) {
			data[k][i][j] = m_vecRecvBuffer[m_ixRecvBuffer++];
		}
		}
		}

	// Unpack data from bottom-left
	} else if (m_dir == Direction_BottomLeft) {
		ixABoundaryBegin = 0;
		ixABoundaryEnd   = data.GetHaloElements();
		ixBBoundaryBegin = 0;
		ixBBoundaryEnd   = data.GetHaloElements();

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  data.GetRElements()
			* (ixABoundaryEnd - ixABoundaryBegin)
			* (ixBBoundaryEnd - ixBBoundaryBegin);

		// Unpack data
		for (int k = 0; k < data.GetRElements(); k++) {
		for (int j = ixBBoundaryBegin; j < ixBBoundaryEnd; j++) {
		for (int i = ixABoundaryBegin; i < ixABoundaryEnd; i++) {
			data[k][i][j] = m_vecRecvBuffer[m_ixRecvBuffer++];
		}
		}
		}

	// Unpack data from bottom-right
	} else if (m_dir == Direction_BottomRight) {
		ixABoundaryBegin = data.GetAElements() - data.GetHaloElements();
		ixABoundaryEnd   = data.GetAElements();
		ixBBoundaryBegin = 0;
		ixBBoundaryEnd   = data.GetHaloElements();

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  data.GetRElements()
			* (ixABoundaryEnd - ixABoundaryBegin)
			* (ixBBoundaryEnd - ixBBoundaryBegin);

		// Unpack data
		for (int k = 0; k < data.GetRElements(); k++) {
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
	GridData4D & data
) {
	// Model grid
	const Grid & grid = m_pConnect->GetGridPatch().GetGrid();

	// 3D data
	GridData3D data3D;

	// List of variable indices to receive
	// - exclude variables which are not-collocated with this data structure
	if (data.GetDataType() == DataType_State) {
		for (int c = 0; c < data.GetComponents(); c++) {
			if (grid.GetVarLocation(c) != data.GetDataLocation()) {
				continue;
			}
			data.GetAsGridData3D(c, data3D);
			Unpack(data3D);
		}

	// Unpack all variables
	} else {
		for (int c = 0; c < data.GetComponents(); c++) {
			data.GetAsGridData3D(c, data3D);
			Unpack(data3D);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
// Connectivity
///////////////////////////////////////////////////////////////////////////////

void Connectivity::ExteriorConnect(
	GridPatch * pPatchFirst,
	Direction dirFirst,
	GridPatch * pPatchSecond,
	Direction dirSecond,
	bool fReverseDirection
) {
	// Get connectivities of individual patches
	Connectivity & connectFirst  = pPatchFirst->GetConnectivity();
	//Connectivity & connectSecond = pPatchSecond->GetConnectivity();

	// First patch coordinate index
	int ixFirstBegin;
	int ixFirstEnd;

	if ((dirFirst == Direction_Right) ||
		(dirFirst == Direction_Left)
	) {
		ixFirstBegin = pPatchFirst->GetPatchBox().GetBInteriorBegin();
		ixFirstEnd   = pPatchFirst->GetPatchBox().GetBInteriorEnd();

	} else if (
		(dirFirst == Direction_Top) ||
		(dirFirst == Direction_Bottom)
	) {
		ixFirstBegin = pPatchFirst->GetPatchBox().GetAInteriorBegin();
		ixFirstEnd   = pPatchFirst->GetPatchBox().GetAInteriorEnd();

	} else {
		ixFirstBegin = pPatchFirst->GetPatchBox().GetHaloElements();
		ixFirstEnd   = pPatchFirst->GetPatchBox().GetHaloElements() + 1;
	}
/*
	// Second patch coordinate index
	int ixSecondBegin;
	int ixSecondEnd;

	if ((dirSecond == Direction_Right) ||
		(dirSecond == Direction_Left)
	) {
		ixSecondBegin = pPatchSecond->GetPatchBox().GetBInteriorBegin();
		ixSecondEnd   = pPatchSecond->GetPatchBox().GetBInteriorEnd();

	} else if (
		(dirSecond == Direction_Top) ||
		(dirSecond == Direction_Bottom)
	) {
		ixSecondBegin = pPatchSecond->GetPatchBox().GetAInteriorBegin();
		ixSecondEnd   = pPatchSecond->GetPatchBox().GetAInteriorEnd();

	} else {
		ixSecondBegin = pPatchSecond->GetPatchBox().GetHaloElements();
		ixSecondEnd   = pPatchSecond->GetPatchBox().GetHaloElements() + 1;
	}
*/
	// Add connectivity to first patch
	connectFirst.m_vecExteriorNeighbors.push_back(
		ExteriorNeighbor(
			&connectFirst,
			dirFirst,
			dirSecond,
			pPatchSecond->GetPatchIndex(),
			fReverseDirection,
			ixFirstBegin,
			ixFirstEnd
		));
/*
	// Add connectivity to second patch
	connectSecond.m_vecExteriorNeighbors.push_back(
		ExteriorNeighbor(
			&connectSecond,
			dirSecond,
			dirFirst,
			pPatchFirst->GetPatchIndex(),
			fReverseDirection,
			ixSecondBegin,
			ixSecondEnd
		));
*/
}

///////////////////////////////////////////////////////////////////////////////

void Connectivity::InitializeBuffers(
	int nRElements,
	int nHaloElements,
	int nVariables
) {
	// Initialize exterior neighbor buffers
	for (int m = 0; m < m_vecExteriorNeighbors.size(); m++) {
		m_vecExteriorNeighbors[m].InitializeBuffers(
			nRElements, nHaloElements, nVariables);
	}

	// Initialize interior neighbor buffers
	for (int m = 0; m < m_vecNestedNeighbors.size(); m++) {
		m_vecNestedNeighbors[m].InitializeBuffers(
			nRElements, nHaloElements, nVariables);
	}
}

///////////////////////////////////////////////////////////////////////////////

void Connectivity::PrepareExchange() {

	// Prepare for asynchronous receives from each neighbor
	for (int m = 0; m < m_vecExteriorNeighbors.size(); m++) {
		m_vecExteriorNeighbors[m].PrepareExchange();
	}
}

///////////////////////////////////////////////////////////////////////////////

void Connectivity::Pack(
	const GridData3D & data
) {
	// Send data to exterior neighbors
	for (int m = 0; m < m_vecExteriorNeighbors.size(); m++) {
		m_vecExteriorNeighbors[m].Pack(data);
	}
}

///////////////////////////////////////////////////////////////////////////////

void Connectivity::Pack(
	const GridData4D & data
) {
	// Send data to exterior neighbors
	for (int m = 0; m < m_vecExteriorNeighbors.size(); m++) {
		m_vecExteriorNeighbors[m].Pack(data);
	}
}

///////////////////////////////////////////////////////////////////////////////

void Connectivity::Send() {
	// Send data to exterior neighbors
	for (int m = 0; m < m_vecExteriorNeighbors.size(); m++) {
		m_vecExteriorNeighbors[m].Send();
	}
}

///////////////////////////////////////////////////////////////////////////////

Neighbor * Connectivity::WaitReceive() {

	// Receive data from exterior neighbors
	int nRecvMessageCount = 0;

	while (nRecvMessageCount != m_vecExteriorNeighbors.size()) {

		nRecvMessageCount = 0;
		for (int m = 0; m < m_vecExteriorNeighbors.size(); m++) {
			if (m_vecExteriorNeighbors[m].IsComplete()) {
				nRecvMessageCount++;
				continue;
			}
			if (m_vecExteriorNeighbors[m].CheckReceive()) {
				m_vecExteriorNeighbors[m].SetComplete();
				return &(m_vecExteriorNeighbors[m]);
			}
		}
	}

	return NULL;
}

///////////////////////////////////////////////////////////////////////////////

int Connectivity::GetExpectedMessageCount() const {
	return (int)(m_vecExteriorNeighbors.size() + m_vecNestedNeighbors.size());
}

///////////////////////////////////////////////////////////////////////////////

