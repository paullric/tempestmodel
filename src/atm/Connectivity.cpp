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
// ExchangeBufferRegistry
///////////////////////////////////////////////////////////////////////////////

ExchangeBufferRegistry::ExchangeBufferRegistry() :
	m_fOwnsData(true)
{ }

///////////////////////////////////////////////////////////////////////////////

ExchangeBufferRegistry::~ExchangeBufferRegistry() {
	if (m_fOwnsData) {
		for (int i = 0; i < m_vecRecvBuffers.size(); i++) {
			if (m_vecRecvBuffers[i] != NULL) {
				delete m_vecRecvBuffers[i];
			}
		}
		for (int i = 0; i < m_vecSendBuffers.size(); i++) {
			if (m_vecSendBuffers[i] != NULL) {
				delete m_vecSendBuffers[i];
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void ExchangeBufferRegistry::SetNoDataOwnership() {
	for (int i = 0; i < m_vecRecvBuffers.size(); i++) {
		if (m_vecRecvBuffers[i] != NULL) {
			_EXCEPTIONT("ExchangeBufferRegistry::RecvBuffers is active");
		}
		if (m_vecSendBuffers[i] != NULL) {
			_EXCEPTIONT("ExchangeBufferRegistry::SendBuffers is active");
		}
	}

	m_fOwnsData = false;
}

///////////////////////////////////////////////////////////////////////////////

void ExchangeBufferRegistry::Register(
	const ExchangeBufferInfo & info
) {
	m_vecRegistry.push_back(info);
	m_vecRecvBuffers.resize(m_vecRegistry.size(), NULL);
	m_vecSendBuffers.resize(m_vecRegistry.size(), NULL);
}

///////////////////////////////////////////////////////////////////////////////

void ExchangeBufferRegistry::GetExchangeBuffersBySourcePatchIx(
	int ixSourcePatch,
	std::vector<int> & vecExchangeBufferIndices
) const {
	vecExchangeBufferIndices.clear();
	for (int m = 0; m < m_vecRegistry.size(); m++) {
		if (m_vecRegistry[m].ixSourcePatch == ixSourcePatch) {
			vecExchangeBufferIndices.push_back(m);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void ExchangeBufferRegistry::Allocate(int ix) {
	if (!m_fOwnsData) {
		_EXCEPTIONT("ExchangeBufferRegistry does not own its data");
	}

	const ExchangeBufferInfo & info = GetExchangeBufferInfo(ix);

	if (m_vecRecvBuffers[ix] != NULL) {
		_EXCEPTIONT("Recv buffer already allocated");
	}
	if (m_vecSendBuffers[ix] != NULL) {
		_EXCEPTIONT("Send buffer already allocated");
	}

	if (info.sByteSize == 0) {
		_EXCEPTIONT("Invalid ByteSize on ExchangeBuffer");
	}

	m_vecRecvBuffers[ix] = new unsigned char[info.sByteSize];
	m_vecSendBuffers[ix] = new unsigned char[info.sByteSize];
}

///////////////////////////////////////////////////////////////////////////////

void ExchangeBufferRegistry::Assign(
	int ix,
	unsigned char * pRecvBuffer,
	unsigned char * pSendBuffer
) {
	if (m_fOwnsData) {
		_EXCEPTIONT("ExchangeBufferRegistry owns its data");
	}
	if ((ix < 0) || (ix > m_vecRegistry.size())) {
		_EXCEPTIONT("ExchangeBuffer index out of range");
	}

	m_vecRecvBuffers[ix] = pRecvBuffer;
	m_vecSendBuffers[ix] = pSendBuffer;
}

///////////////////////////////////////////////////////////////////////////////

void ExchangeBufferRegistry::Unassign(
	int ix
) { 
	if (m_fOwnsData) {
		_EXCEPTIONT("ExchangeBufferRegistry owns its data");
	}
	if ((ix < 0) || (ix > m_vecRegistry.size())) {
		_EXCEPTIONT("ExchangeBuffer index out of range");
	}

	m_vecRecvBuffers[ix] = NULL;
	m_vecSendBuffers[ix] = NULL;
}

///////////////////////////////////////////////////////////////////////////////
// Neighbor
///////////////////////////////////////////////////////////////////////////////

Neighbor::Neighbor(
	const ExchangeBufferInfo & info
) : 
	m_info(info),
	m_fComplete(false),
	m_ixSendBuffer(0),
	m_ixRecvBuffer(0)
{
	m_vecRecvBuffer.SetSize(
		  m_info.sBoundarySize
		* m_info.sMaxRElements
		* m_info.sHaloElements
		* m_info.sComponents);

	m_vecSendBuffer.SetSize(
		  m_info.sBoundarySize
		* m_info.sMaxRElements
		* m_info.sHaloElements
		* m_info.sComponents);
}

///////////////////////////////////////////////////////////////////////////////

void Neighbor::AllocateBuffers() {
	m_vecRecvBuffer.Allocate();
	m_vecSendBuffer.Allocate();
}

///////////////////////////////////////////////////////////////////////////////

void Neighbor::AttachBuffers(
	unsigned char * pRecvBuffer,
	unsigned char * pSendBuffer
) {
	if (pRecvBuffer == NULL) {
		_EXCEPTIONT("Invalid RecvBuffer pointer");
	}
	if (pSendBuffer == NULL) {
		_EXCEPTIONT("Invalid SendBuffer pointer");
	}

	m_vecRecvBuffer.AttachToData(pRecvBuffer);
	m_vecSendBuffer.AttachToData(pSendBuffer);
}

///////////////////////////////////////////////////////////////////////////////

bool Neighbor::CheckReceive() {

	// Check if message already received and processed
	if (m_fComplete) {
		return false;
	}

#ifdef TEMPEST_MPIOMP
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

#ifdef TEMPEST_MPIOMP
	// Information for receive
	int iProcessor = grid.GetPatchProcessor(m_info.ixTargetPatch);

	// Tag of the received message
	int nTag = (m_info.ixTargetPatch << 16) + (m_info.ixSourcePatch << 4) + (int)(m_info.dir);

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
	const DataArray3D<double> & data
) {
	const size_t sRElements = data.GetSize(0);
	const size_t sAElements = data.GetSize(1);
	const size_t sBElements = data.GetSize(2);

	// Check matrix bounds
	if (((m_info.dir == Direction_Right) || (m_info.dir == Direction_Left)) &&
		(m_info.ixSecond > sBElements)
	) {
		_EXCEPTIONT("GridData / ExteriorNeighbor inconsistency.");
	}
	if (((m_info.dir == Direction_Top) || (m_info.dir == Direction_Bottom)) &&
		(m_info.ixSecond > sAElements)
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
	if (m_info.dir == Direction_Right) {
		ixBoundaryBegin = sAElements - 2 * m_info.sHaloElements;
		ixBoundaryEnd   = sAElements - m_info.sHaloElements;

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixBoundaryEnd - ixBoundaryBegin)
			* (m_info.ixSecond - m_info.ixFirst);

		// Pack the SendBuffer
		if (m_info.fReverseDirection) {
			for (int k = 0; k < sRElements; k++) {
			for (int i = ixBoundaryBegin; i < ixBoundaryEnd; i++) {
			for (int j = m_info.ixSecond-1; j >= m_info.ixFirst; j--) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}

		} else {
			for (int k = 0; k < sRElements; k++) {
			for (int i = ixBoundaryBegin; i < ixBoundaryEnd; i++) {
			for (int j = m_info.ixFirst; j < m_info.ixSecond; j++) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];	
			}
			}
			}
		}
	
	// Pack data to send topward
	} else if (m_info.dir == Direction_Top) {
		ixBoundaryBegin = sBElements - 2 * m_info.sHaloElements;
		ixBoundaryEnd   = sBElements - m_info.sHaloElements;

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixBoundaryEnd - ixBoundaryBegin)
			* (m_info.ixSecond - m_info.ixFirst);

		if (m_ixSendBuffer + nTotalValues > m_vecSendBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in SendBuffer for operation.");
		}

		// Pack the SendBuffer
		if (m_info.fReverseDirection) {
			for (int k = 0; k < sRElements; k++) {
			for (int j = ixBoundaryBegin; j < ixBoundaryEnd; j++) {
			for (int i = m_info.ixSecond-1; i >= m_info.ixFirst; i--) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}

		} else {
			for (int k = 0; k < sRElements; k++) {
			for (int j = ixBoundaryBegin; j < ixBoundaryEnd; j++) {
			for (int i = m_info.ixFirst; i < m_info.ixSecond; i++) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}
		}

	// Pack data to send left
	} else if (m_info.dir == Direction_Left) {
		ixBoundaryBegin = m_info.sHaloElements;
		ixBoundaryEnd   = 2 * m_info.sHaloElements;

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixBoundaryEnd - ixBoundaryBegin)
			* (m_info.ixSecond - m_info.ixFirst);

		if (m_ixSendBuffer + nTotalValues > m_vecSendBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in SendBuffer for operation.");
		}

		// Pack the SendBuffer
		if (m_info.fReverseDirection) {
			for (int k = 0; k < sRElements; k++) {
			for (int i = ixBoundaryEnd-1; i >= ixBoundaryBegin; i--) {
			for (int j = m_info.ixSecond-1; j >= m_info.ixFirst; j--) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}

		} else {
			for (int k = 0; k < sRElements; k++) {
			for (int i = ixBoundaryEnd-1; i >= ixBoundaryBegin; i--) {
			for (int j = m_info.ixFirst; j < m_info.ixSecond; j++) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}
		}

	// Pack data to send bottomward
	} else if (m_info.dir == Direction_Bottom) {
		ixBoundaryBegin = m_info.sHaloElements;
		ixBoundaryEnd   = 2 * m_info.sHaloElements;

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixBoundaryEnd - ixBoundaryBegin)
			* (m_info.ixSecond - m_info.ixFirst);

		if (m_ixSendBuffer + nTotalValues > m_vecSendBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in SendBuffer for operation.");
		}

		// Pack the SendBuffer
		if (m_info.fReverseDirection) {
			for (int k = 0; k < sRElements; k++) {
			for (int j = ixBoundaryEnd-1; j >= ixBoundaryBegin; j--) {
			for (int i = m_info.ixSecond-1; i >= m_info.ixFirst; i--) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}

		} else {
			for (int k = 0; k < sRElements; k++) {
			for (int j = ixBoundaryEnd-1; j >= ixBoundaryBegin; j--) {
			for (int i = m_info.ixFirst; i < m_info.ixSecond; i++) {
				m_vecSendBuffer[m_ixSendBuffer++] = data[k][i][j];
			}
			}
			}
		}

	// Pack data to send toprightward 
	} else if (m_info.dir == Direction_TopRight) {
		ixABoundaryBegin = m_info.ixFirst - m_info.sHaloElements + 1;
		ixABoundaryEnd   = m_info.ixFirst + 1;
		ixBBoundaryBegin = m_info.ixSecond - m_info.sHaloElements + 1;
		ixBBoundaryEnd   = m_info.ixSecond + 1;
 
		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixABoundaryEnd - ixABoundaryBegin)
			* (ixBBoundaryEnd - ixBBoundaryBegin);

		if (m_ixSendBuffer + nTotalValues > m_vecSendBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in SendBuffer for operation.");
		}

		// Pack the SendBuffer
		if (m_info.fReverseDirection) {
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
	} else if (m_info.dir == Direction_TopLeft) {
		ixABoundaryBegin = m_info.ixFirst;
		ixABoundaryEnd   = m_info.ixFirst + m_info.sHaloElements;
		ixBBoundaryBegin = m_info.ixSecond - m_info.sHaloElements + 1;
		ixBBoundaryEnd   = m_info.ixSecond + 1;
 
		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixABoundaryEnd - ixABoundaryBegin)
			* (ixBBoundaryEnd - ixBBoundaryBegin);

		if (m_ixSendBuffer + nTotalValues > m_vecSendBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in SendBuffer for operation.");
		}

		// Pack the SendBuffer
		if (m_info.fReverseDirection) {
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
	} else if (m_info.dir == Direction_BottomLeft) {
		ixABoundaryBegin = m_info.ixFirst;
		ixABoundaryEnd   = m_info.ixFirst + m_info.sHaloElements;
		ixBBoundaryBegin = m_info.ixSecond;
		ixBBoundaryEnd   = m_info.ixSecond + m_info.sHaloElements;
 
		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixABoundaryEnd - ixABoundaryBegin)
			* (ixBBoundaryEnd - ixBBoundaryBegin);

		if (m_ixSendBuffer + nTotalValues > m_vecSendBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in SendBuffer for operation.");
		}

		// Pack the SendBuffer
		if (m_info.fReverseDirection) {
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
	} else if (m_info.dir == Direction_BottomRight) {
		ixABoundaryBegin = m_info.ixFirst - m_info.sHaloElements + 1;
		ixABoundaryEnd   = m_info.ixFirst + 1;
		ixBBoundaryBegin = m_info.ixSecond;
		ixBBoundaryEnd   = m_info.ixSecond + m_info.sHaloElements;
 
		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixABoundaryEnd - ixABoundaryBegin)
			* (ixBBoundaryEnd - ixBBoundaryBegin);

		if (m_ixSendBuffer + nTotalValues > m_vecSendBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in SendBuffer for operation.");
		}

		// Pack the SendBuffer
		if (m_info.fReverseDirection) {
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

void ExteriorNeighbor::Send(
	const Grid & grid
) {

#ifdef TEMPEST_MPIOMP
	// Send the data
	int iProcessor = grid.GetPatchProcessor(m_info.ixTargetPatch);

	// Tag
	if (m_info.ixTargetPatch >= (2 << 12)) {
		_EXCEPTIONT("Maximum neighbor index exceeded.");
	}

	int nTag = (m_info.ixSourcePatch << 16) + (m_info.ixTargetPatch << 4) + (int)(m_info.dirOpposing);

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
	if (m_info.dir == Direction_Right) {
		ixBoundaryBegin = sAElements - m_info.sHaloElements;
		ixBoundaryEnd   = sAElements;

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixBoundaryEnd - ixBoundaryBegin)
			* (m_info.ixSecond - m_info.ixFirst);

		if (m_ixRecvBuffer + nTotalValues > m_vecRecvBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in RecvBuffer for operation.");
		}

		// Unpack data
		for (int k = 0; k < sRElements; k++) {
		for (int i = ixBoundaryEnd-1; i >= ixBoundaryBegin; i--) {
		for (int j = m_info.ixFirst; j < m_info.ixSecond; j++) {
			data[k][i][j] = m_vecRecvBuffer[m_ixRecvBuffer++];
		}
		}
		}

	// Unpack data from top
	} else if (m_info.dir == Direction_Top) {
		ixBoundaryBegin = sBElements - m_info.sHaloElements;
		ixBoundaryEnd   = sBElements;

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixBoundaryEnd - ixBoundaryBegin)
			* (m_info.ixSecond - m_info.ixFirst);

		if (m_ixRecvBuffer + nTotalValues > m_vecRecvBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in RecvBuffer for operation.");
		}

		// Unpack data
		for (int k = 0; k < sRElements; k++) {
		for (int j = ixBoundaryEnd-1; j >= ixBoundaryBegin; j--) {
		for (int i = m_info.ixFirst; i < m_info.ixSecond; i++) {
			data[k][i][j] = m_vecRecvBuffer[m_ixRecvBuffer++];
		}
		}
		}

	// Unpack data from left
	} else if (m_info.dir == Direction_Left) {
		ixBoundaryBegin = 0;
		ixBoundaryEnd   = m_info.sHaloElements;

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixBoundaryEnd - ixBoundaryBegin)
			* (m_info.ixSecond - m_info.ixFirst);

		if (m_ixRecvBuffer + nTotalValues > m_vecRecvBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in RecvBuffer for operation.");
		}

		// Unpack data
		for (int k = 0; k < sRElements; k++) {
		for (int i = ixBoundaryBegin; i < ixBoundaryEnd; i++) {
		for (int j = m_info.ixFirst; j < m_info.ixSecond; j++) {
			data[k][i][j] = m_vecRecvBuffer[m_ixRecvBuffer++];
		}
		}
		}

	// Unpack data from bottom
	} else if (m_info.dir == Direction_Bottom) {
		ixBoundaryBegin = 0;
		ixBoundaryEnd   = m_info.sHaloElements;

		// Check that sufficient data remains in send buffer
		int nTotalValues =
			  sRElements
			* (ixBoundaryEnd - ixBoundaryBegin)
			* (m_info.ixSecond - m_info.ixFirst);

		if (m_ixRecvBuffer + nTotalValues > m_vecRecvBuffer.GetRows()) {
			_EXCEPTIONT("Insufficient space in RecvBuffer for operation.");
		}

		// Unpack data
		for (int k = 0; k < sRElements; k++) {
		for (int j = ixBoundaryBegin; j < ixBoundaryEnd; j++) {
		for (int i = m_info.ixFirst; i < m_info.ixSecond; i++) {
			data[k][i][j] = m_vecRecvBuffer[m_ixRecvBuffer++];
		}
		}
		}

	// Unpack data from top-right
	} else if (m_info.dir == Direction_TopRight) {
		ixABoundaryBegin = m_info.ixFirst + 1;
		ixABoundaryEnd   = m_info.ixFirst + m_info.sHaloElements + 1;
		ixBBoundaryBegin = m_info.ixSecond + 1;
		ixBBoundaryEnd   = m_info.ixSecond + m_info.sHaloElements + 1;

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
	} else if (m_info.dir == Direction_TopLeft) {
		ixABoundaryBegin = m_info.ixFirst - m_info.sHaloElements;
		ixABoundaryEnd   = m_info.ixFirst;
		ixBBoundaryBegin = m_info.ixSecond + 1;
		ixBBoundaryEnd   = m_info.ixSecond + m_info.sHaloElements + 1;

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
	} else if (m_info.dir == Direction_BottomLeft) {
		ixABoundaryBegin = m_info.ixFirst - m_info.sHaloElements;
		ixABoundaryEnd   = m_info.ixFirst;
		ixBBoundaryBegin = m_info.ixSecond - m_info.sHaloElements;
		ixBBoundaryEnd   = m_info.ixSecond;

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
	} else if (m_info.dir == Direction_BottomRight) {
		ixABoundaryBegin = m_info.ixFirst + 1;
		ixABoundaryEnd   = m_info.ixFirst + m_info.sHaloElements + 1;
		ixBBoundaryBegin = m_info.ixSecond - m_info.sHaloElements;
		ixBBoundaryEnd   = m_info.ixSecond;

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

void ExteriorNeighbor::WaitSend() {

#ifdef TEMPEST_MPIOMP
	MPI_Status status;
	MPI_Wait(&m_reqSend, &status);
#endif

}

///////////////////////////////////////////////////////////////////////////////
// Connectivity
///////////////////////////////////////////////////////////////////////////////

Connectivity::~Connectivity() {
	ClearNeighbors();
}

///////////////////////////////////////////////////////////////////////////////

void Connectivity::ClearNeighbors() {
	for (int n = 0; n < m_vecExteriorNeighbors.size(); n++) {
		delete m_vecExteriorNeighbors[n];
	}
	m_vecExteriorNeighbors.resize(0);
/*
	for (int n = 0; n < m_vecNestedNeighbors; n++) {
		delete m_vecNestedNeighbors[n];
	}
	m_vecNestedNeighbors.resize(0);
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

		int iDir = static_cast<int>(pNeighbor->m_info.dir);

		// Right or Left side of this PatchBox
		if ((pNeighbor->m_info.dir == Direction_Right) || 
		    (pNeighbor->m_info.dir == Direction_Left)
		) {
			int j = pNeighbor->m_info.ixFirst;
			for (; j < pNeighbor->m_info.ixSecond; j++) {
				if ((j < box.GetBInteriorBegin()) ||
					(j >= box.GetBInteriorEnd())
				) {
					_EXCEPTIONT("Edge index out of range");
				}

				m_vecExteriorEdge[iDir][j] = pNeighbor;
			}

		// Top or Bottom side of this PatchBox
		} else if (
		    (pNeighbor->m_info.dir == Direction_Top) ||
		    (pNeighbor->m_info.dir == Direction_Bottom)
		) {
			int i = pNeighbor->m_info.ixFirst;
			for (; i < pNeighbor->m_info.ixSecond; i++) {
				if ((i < box.GetAInteriorBegin()) ||
					(i >= box.GetAInteriorEnd())
				) {
					_EXCEPTIONT("Edge index out of range");
				}

				m_vecExteriorEdge[iDir][i] = pNeighbor;
			}

		// Corners of this PatchBox
		} else if (
		    (pNeighbor->m_info.dir == Direction_TopRight) ||
		    (pNeighbor->m_info.dir == Direction_TopLeft) ||
		    (pNeighbor->m_info.dir == Direction_BottomRight) ||
		    (pNeighbor->m_info.dir == Direction_BottomLeft)
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
		m_vecExteriorNeighbors[m]->Pack(data);
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

