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

#ifndef _CONNECTIVITY_H_
#define _CONNECTIVITY_H_

///////////////////////////////////////////////////////////////////////////////

#include "Direction.h"

#include "DataContainer.h"
#include "DataArray1D.h"
#include "DataArray3D.h"
#include "DataArray4D.h"

#ifdef TEMPEST_MPIOMP
#include <mpi.h>
#endif

#include <vector>

///////////////////////////////////////////////////////////////////////////////

class Grid;
class GridPatch;
class Neighbor;
class Connectivity;

///////////////////////////////////////////////////////////////////////////////

#ifdef TEMPEST_MPIOMP
typedef int ExchangeBufferId;
#endif
#ifdef TEMPEST_OCR
typedef int ExchangeBufferId;
#endif
#ifdef TEMPEST_HPX
typedef int ExchangeBufferId;
#endif

///////////////////////////////////////////////////////////////////////////////

class ExchangeBufferRegistry;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Class containing information on buffers used for exchange of data
///		between GridPatches.
///	</summary>
class ExchangeBuffer {

	friend class ExchangeBufferRegistry;

public:
	///	<summary>
	///		Message header for MPI message exchange.
	///	</summary>
	struct MessageHeader {
		public:
			///	<summary>
			///		Default constructor.
			///	</summary>
			MessageHeader() :
				m_ixReserved(0x01010101),
				m_ixFirstPatch(-1),
				m_ixSecondPatch(-1),
				m_ixDirection(Direction_Middle)
			{ }

			///	<summary>
			///		Constructor.
			///	</summary>
			MessageHeader(
				int ixFirstPatch,
				int ixSecondPatch,
				int ixDirection
			) :
				m_ixReserved(0x01010101),
				m_ixFirstPatch(ixFirstPatch),
				m_ixSecondPatch(ixSecondPatch),
				m_ixDirection(ixDirection)
			{ }

			///	<summary>
			///		Equality comparator.
			///	</summary>
			bool operator==(const MessageHeader & msghead) const {
				if ((m_ixReserved == msghead.m_ixReserved) &&
				    (m_ixFirstPatch == msghead.m_ixFirstPatch) &&
				    (m_ixSecondPatch == msghead.m_ixSecondPatch) &&
				    (m_ixDirection == msghead.m_ixDirection)
				) {
					return true;
				}
				return false;
			}

			///	<summary>
			///		Inequality comparator.
			///	</summary>
			bool operator!=(const MessageHeader & msghead) const {
				return !((*this) == msghead);
			}

		public:
			int m_ixReserved;
			int m_ixFirstPatch;
			int m_ixSecondPatch;
			int m_ixDirection;
	};

public:
	///	<summary>
	///		Default constructor.
	///	</summary>
	ExchangeBuffer() :
		m_ixSourcePatch(-1),
		m_ixTargetPatch(-1),
		m_dir(Direction_Middle),
		m_dirOpposing(Direction_Middle),
		m_ixSourceProcessor(-1),
		m_ixTargetProcessor(-1),
		m_ixLocalActiveSourcePatch(-1),
		m_sHaloElements(0),
		m_sComponents(0),
		m_sMaxRElements(0),
		m_sBoundarySize(0),
		m_sByteSize(0),
		m_ixFirst(0),
		m_ixSecond(0),
		m_fReverseDirection(false),
		m_fFlippedCoordinate(false)
	{ }

public:
	///	<summary>
	///		Calculate total data size (in bytes).
	///	</summary>
	void CalculateByteSize() {
		m_sByteSize =
			  m_sBoundarySize
			* m_sMaxRElements
			* m_sHaloElements
			* m_sComponents
			* sizeof(double);
	}

	///	<summary>
	///		Get the message size (in bytes).
	///	</summary>
	size_t GetMessageSize() {
		return (m_sByteSize + sizeof(MessageHeader));
	}

	///	<summary>
	///		Get the outgoing message header for this ExchangeBuffer.
	///	</summary>
	void GetSendMessageHeader(
		MessageHeader * pMessageHead
	) {
		(*pMessageHead) =
			MessageHeader(
				m_ixSourcePatch,
				m_ixTargetPatch,
				m_dirOpposing);
	}

	///	<summary>
	///		Get the incoming message header for this ExchangeBuffer.
	///	</summary>
	void GetRecvMessageHeader(
		MessageHeader * pMessageHead
	) {
		(*pMessageHead) =
			MessageHeader(
				m_ixTargetPatch,
				m_ixSourcePatch,
				m_dir);
	}

public:
	///	<summary>
	///		Reset indices and pack the MessageHeader into the given buffer.
	///	</summary>
	void Reset();

	///	<summary>
	///		Pack DataArray3D into the send buffer.
	///	</summary>
	void Pack(
		const DataArray3D<double> & data
	);

	///	<summary>
	///		Pack DataArray4D into the send buffer.
	///	</summary>
	void Pack(
		const Grid & grid,
		const DataArray4D<double> & data
	);

	///	<summary>
	///		Unpack receive buffer into the given DataArray3D.
	///	</summary>
	void Unpack(
		DataArray3D<double> & data
	);

	///	<summary>
	///		Unpack receive buffer into the given DataArray4D.
	///	</summary>
	void Unpack(
		const Grid & grid,
		DataArray4D<double> & data
	);

public:
	///	<summary>
	///		Unique global id of this exchange buffer.
	///	</summary>
	ExchangeBufferId m_id;

	///	<summary>
	///		Source GridPatch index.
	///	</summary>
	int m_ixSourcePatch;

	///	<summary>
	///		Target GridPatch index.
	///	</summary>
	int m_ixTargetPatch;

	///	<summary>
	///		Direction.
	///	</summary>
	Direction m_dir;

	///	<summary>
	///		Opposing direction.
	///	</summary>
	Direction m_dirOpposing;

	///	<summary>
	///		Source processor index.
	///	</summary>
	int m_ixSourceProcessor;

	///	<summary>
	///		Target processor index.
	///	</summary>
	int m_ixTargetProcessor;

	///	<summary>
	///		Local active GridPatch index.
	///	</summary>
	int m_ixLocalActiveSourcePatch;

	///	<summary>
	///		Number of halo elements.
	///	</summary>
	size_t m_sHaloElements;

	///	<summary>
	///		Number of components.
	///	</summary>
	size_t m_sComponents;

	///	<summary>
	///		Number of radial elements.
	///	</summary>
	size_t m_sMaxRElements;

	///	<summary>
	///		Boundary size (lateral length)
	///	</summary>
	size_t m_sBoundarySize;

	///	<summary>
	///		Total data size (in bytes).
	///	</summary>
	size_t m_sByteSize;

	///	<summary>
	///		First local index (begin node for Right, Top, Left, Bottom and
	///		alpha coordinate for TopRight, TopLeft, BottomLeft, BottomRight)
	///	</summary>
	int m_ixFirst;

	///	<summary>
	///		Second local index (end node for Right, Top, Left, Bottom and
	///		beta coordinate for TopRight, TopLeft, BottomLeft, BottomRight)
	///	</summary>
	int m_ixSecond;

	///	<summary>
	///		Indicator of whether or not directions are flipped across the edge.
	///	</summary>
	bool m_fReverseDirection;

	///	<summary>
	///		Indicator of whether or not the direction of increasing coordinate
	///		indices flips across the edge.
	///	</summary>
	bool m_fFlippedCoordinate;

protected:
	///	<summary>
	///		Current RecvBuffer.
	///	</summary>
	DataArray1D<double> m_dRecvBuffer;

	///	<summary>
	///		Current SendBuffer.
	///	</summary>
	DataArray1D<double> m_dSendBuffer;

	///	<summary>
	///		Current RecvBuffer index.
	///	</summary>
	int m_ixRecvBuffer;

	///	<summary>
	///		Current SendBuffer index.
	///	</summary>
	int m_ixSendBuffer;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Class storing all local ExchangeBuffers.
///	</summary>
class ExchangeBufferRegistry {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	ExchangeBufferRegistry();

	///	<summary>
	///		Destructor.
	///	</summary>
	virtual ~ExchangeBufferRegistry();

public:
	///	<summary>
	///		Register a new ExchangeBuffer.
	///	</summary>
	void Register(
		const ExchangeBuffer & exbuf
	);

	///	<summary>
	///		Get the vector of ExchangeBuffers.
	///	</summary>
	std::vector<ExchangeBuffer> & GetExchangeBuffers() {
		return m_vecRegistry;
	}

	///	<summary>
	///		Allocate ExchangeBuffer storage.
	///	</summary>
	void Allocate();

	///	<summary>
	///		Set up asynchronous receives.
	///	</summary>
	void PrepareExchange();

	///	<summary>
	///		Set up asynchronous sends.
	///	</summary>
	void Send();

protected:
	///	<summary>
	///		Attach RecvBuffers based on a received message.
	///	</summary>
	void AttachRecvBuffers(int ixProc);

public:
	///	<summary>
	///		Wait for a message to be received.
	///	</summary>
	///	<returns>
	///		A vector of pointers to ExchangeBuffers that have been filled,
	///		or NULL if all messages have been received.
	///	</returns>
	const std::vector<ExchangeBuffer *> * WaitReceive();

	///	<summary>
	///		Wait for all sends to complete.
	///	</summary>
	void WaitSend();

protected:
	///	<summary>
	///		Flag indicating that ExchangeBuffer RecvBuffers have not yet been
	///		attached to data.
	///	</summary>
	std::vector<bool> m_vecAreRecvBuffersAttached;

	///	<summary>
	///		Vector of ExchangeBuffers.
	///	</summary>
	std::vector<ExchangeBuffer> m_vecRegistry;

	///	<summary>
	///		Buffer size for each processor.
	///	</summary>
	std::vector<int> m_vecBufferSize;

	///	<summary>
	///		Processor associated with each buffer.
	///	</summary>
	std::vector<int> m_vecProcessors;

	///	<summary>
	///		Vector of pointers to receive buffers.
	///	</summary>
	std::vector<char *> m_vecRecvBuffers;

	///	<summary>
	///		Vector of pointers to send buffers.
	///	</summary>
	std::vector<char *> m_vecSendBuffers;

	///	<summary>
	///		Vector of MPI_Requests.
	///	</summary>
	std::vector<MPI_Request> m_vecRecvRequest;

	///	<summary>
	///		Vector of MPI_Requests.
	///	</summary>
	std::vector<MPI_Request> m_vecSendRequest;

	///	<summary>
	///		Vector of booleans indicating which processors have sent messages.
	///	</summary>
	std::vector<bool> m_vecMessageReceived;

protected:
	///	<summary>
	///		A lookup table mapping processor to ExchangeBuffer pointers.
	///	</summary>
	typedef std::vector< std::vector<ExchangeBuffer *> > ProcessorExBufferVector;

	///	<summary>
	///		Vector of ExchangeBuffers organized by processor.
	///	</summary>
	ProcessorExBufferVector m_vecRegistryByProcessor;

};

///////////////////////////////////////////////////////////////////////////////

#endif

