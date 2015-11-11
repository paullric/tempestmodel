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

#ifndef _CONNECTIVITY_H_
#define _CONNECTIVITY_H_

///////////////////////////////////////////////////////////////////////////////

#include "Direction.h"

#include "DataContainer.h"
#include "DataArray1D.h"
#include "DataArray3D.h"
#include "DataArray4D.h"

#ifdef USE_MPI
#include "mpi.h"
#endif

#include <vector>

///////////////////////////////////////////////////////////////////////////////

class Grid;
class GridPatch;
class Neighbor;
class Connectivity;

///////////////////////////////////////////////////////////////////////////////

#ifdef USE_MPI
typedef int ExchangeBufferId;
#endif
#ifdef USE_OCR
typedef int ExchangeBufferId;
#endif
#ifdef USE_HPX
typedef int ExchangeBufferId;
#endif

///////////////////////////////////////////////////////////////////////////////

class ExchangeBufferInfo {

public:
	///	<summary>
	///		Default constructor.
	///	</summary>
	ExchangeBufferInfo() :
		ixSourcePatch(-1),
		ixTargetPatch(-1),
		dir(Direction_Middle),
		dirOpposing(Direction_Middle),
		sHaloElements(0),
		sComponents(0),
		sMaxRElements(0),
		sBoundarySize(0),
		sByteSize(0),
		ixFirst(0),
		ixSecond(0),
		fReverseDirection(false),
		fFlippedCoordinate(false)
	{ }

public:
	///	<summary>
	///		Calculate total data size (in bytes).
	///	</summary>
	void CalculateByteSize() {
		sByteSize =
		  sBoundarySize
		* sMaxRElements
		* sHaloElements
		* sComponents
		* sizeof(double);
	}

public:
	///	<summary>
	///		Unique global id of this exchange buffer.
	///	</summary>
	ExchangeBufferId id;

	///	<summary>
	///		Source GridPatch index.
	///	</summary>
	int ixSourcePatch;

	///	<summary>
	///		Target GridPatch index.
	///	</summary>
	int ixTargetPatch;

	///	<summary>
	///		Direction.
	///	</summary>
	Direction dir;

	///	<summary>
	///		Opposing direction.
	///	</summary>
	Direction dirOpposing;

	///	<summary>
	///		Number of halo elements.
	///	</summary>
	size_t sHaloElements;

	///	<summary>
	///		Number of components.
	///	</summary>
	size_t sComponents;

	///	<summary>
	///		Number of radial elements.
	///	</summary>
	size_t sMaxRElements;

	///	<summary>
	///		Boundary size (lateral length)
	///	</summary>
	size_t sBoundarySize;

	///	<summary>
	///		Total data size (in bytes).
	///	</summary>
	size_t sByteSize;

	///	<summary>
	///		First local index (begin node for Right, Top, Left, Bottom and
	///		alpha coordinate for TopRight, TopLeft, BottomLeft, BottomRight)
	///	</summary>
	int ixFirst;

	///	<summary>
	///		Second local index (end node for Right, Top, Left, Bottom and
	///		beta coordinate for TopRight, TopLeft, BottomLeft, BottomRight)
	///	</summary>
	int ixSecond;

	///	<summary>
	///		Indicator of whether or not directions are flipped across the edge.
	///	</summary>
	bool fReverseDirection;

	///	<summary>
	///		Indicator of whether or not the direction of increasing coordinate
	///		indices flips across the edge.
	///	</summary>
	bool fFlippedCoordinate;

};

///////////////////////////////////////////////////////////////////////////////

class ExchangeBufferRegistry {

public:
	///	<summary>
	///		Vector of ExchangeBufferInfo objects.
	///	</summary>
	typedef std::vector<ExchangeBufferInfo> ExchangeBufferInfoVector;

	///	<summary>
	///		Vector of pointers to data.
	///	</summary>
	typedef std::vector<unsigned char *> DataPointerVector;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	ExchangeBufferRegistry();

	///	<summary>
	///		Destructor.
	///	</summary>
	virtual ~ExchangeBufferRegistry();

	///	<summary>
	///		Reclassify this ExchangeBuffer as one that does not own
	///		its data.
	///	</summary>
	void SetNoDataOwnership();

public:
	///	<summary>
	///		Get a reference to the registry.
	///	</summary>
	const ExchangeBufferInfoVector & GetRegistry() const {
		return m_vecRegistry;
	}

	///	<summary>
	///		Get the ExchangeBufferInfo object stored at index ix in the
	///		registry.
	///	</summary>
	const ExchangeBufferInfo & GetExchangeBufferInfo(int ix) {
		if ((ix < 0) || (ix > m_vecRegistry.size())) {
			_EXCEPTIONT("ExchangeBuffer index out of range");
		}
		return m_vecRegistry[ix];
	}

	///	<summary>
	///		Get the pointer for the specified recv buffer.
	///	</summary>
	unsigned char * GetRecvBufferPtr(int ix) {
		if ((ix < 0) || (ix > m_vecRegistry.size())) {
			_EXCEPTIONT("ExchangeBuffer index out of range");
		}
		return m_vecRecvBuffers[ix];
	}

	///	<summary>
	///		Get the pointer for the specified send buffer.
	///	</summary>
	unsigned char * GetSendBufferPtr(int ix) {
		if ((ix < 0) || (ix > m_vecRegistry.size())) {
			_EXCEPTIONT("ExchangeBuffer index out of range");
		}
		return m_vecSendBuffers[ix];
	}

public:
	///	<summary>
	///		Register a new ExchangeBuffer.
	///	</summary>
	void Register(
		const ExchangeBufferInfo & info
	);

	///	<summary>
	///		Get all ExchangeBuffers with source patch index ixPatch.
	///	</summary>
	void GetExchangeBuffersBySourcePatchIx(
		int ixSourcePatch,
		std::vector<int> & vecExchangeBufferIndices
	) const;

	///	<summary>
	///		Allocate the specified recv/send buffers.
	///	</summary>
	void Allocate(int ix);

	///	<summary>
	///		Assign pointers to specified recv/send buffers.
	///	</summary>
	void Assign(
		int ix,
		unsigned char * pRecvBuffer,
		unsigned char * pSendBuffer
	);

	///	<summary>
	///		Unassign pointers to specified recv/send buffers.
	///	</summary>
	void Unassign(int ix);

protected:
	///	<summary>
	///		A flag indicating this ExchangeBufferRegistry owns its own data.
	///	</summary>
	bool m_fOwnsData;

	///	<summary>
	///		Vector of exchange buffers.
	///	</summary>
	ExchangeBufferInfoVector m_vecRegistry;

	///	<summary>
	///		Vector of pointers to receive buffers.
	///	</summary>
	DataPointerVector m_vecRecvBuffers;

	///	<summary>
	///		Vector of pointers to send buffers.
	///	</summary>
	DataPointerVector m_vecSendBuffers;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Interface for GridPatch neighbors.
///	</summary>
class Neighbor {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	Neighbor(
		const ExchangeBufferInfo & info
	);

	///	<summary>
	///		Virtual destructor.
	///	</summary>
	virtual ~Neighbor() {
	}

public:
	///	<summary>
	///		Allocate the exchange buffers.
	///	</summary>
	void AllocateBuffers();

	///	<summary>
	///		Initialize the exchange buffers by pointer.
	///	</summary>
	void AttachBuffers(
		unsigned char * pRecvBuffer,
		unsigned char * pSendBuffer
	);

public:
	///	<summary>
	///		Check if data has been received by this neighbor.
	///	</summary>
	bool CheckReceive();

	///	<summary>
	///		Set the complete flag.
	///	</summary>
	void SetComplete() {
		if (m_fComplete) {
			_EXCEPTIONT("Logic error");
		}
		m_fComplete = true;
	}

public:
	///	<summary>
	///		Prepare an asynchronous receive.
	///	</summary>
	virtual void PrepareExchange(
		const Grid & grid
	) {
		// Reset the send buffer index
		m_ixSendBuffer = 0;

		// Reset the receive buffer index
		m_ixRecvBuffer = 0;

		// Reset the complete flag
		m_fComplete = false;
	}

	///	<summary>
	///		Pack data into the send buffer.
	///	</summary>
	virtual void Pack(
		const DataArray3D<double> & data
	) = 0;

	///	<summary>
	///		Pack data into the send buffer.
	///	</summary>
	virtual void Pack(
		const Grid & grid,
		const DataArray4D<double> & data
	) = 0;

	///	<summary>
	///		Reset the send buffer size.
	///	</summary>
	void ResetSendBufferSize() {
		m_ixSendBuffer =
			  m_info.sBoundarySize
			* m_info.sMaxRElements
			* m_info.sHaloElements
			* m_info.sComponents;
	}

	///	<summary>
	///		Send data to other processors.
	///	</summary>
	virtual void Send(const Grid & grid) = 0;

	///	<summary>
	///		Unpack data from the receive buffer.
	///	</summary>
	virtual void Unpack(
		DataArray3D<double> & data
	) = 0;

	///	<summary>
	///		Unpack data from the receive buffer.
	///	</summary>
	virtual void Unpack(
		const Grid & grid,
		DataArray4D<double> & data
	) = 0;

	///	<summary>
	///		Wait for the send request to complete.
	///	</summary>
	virtual void WaitSend() = 0;

public:
	///	<summary>
	///		Returns true if this neighbor is complete.
	///	</summary>
	bool IsComplete() const {
		return m_fComplete;
	}

public:
	///	<summary>
	///		Exchange buffer information.
	///	</summary>
	ExchangeBufferInfo m_info;

	///	<summary>
	///		Flag indicating if this neighbor has been received and processed.
	///	</summary>
	bool m_fComplete;

	///	<summary>
	///		Index into the send buffer.
	///	</summary>
	int m_ixSendBuffer;

	///	<summary>
	///		Index into the receive buffer.
	///	</summary>
	int m_ixRecvBuffer;

#ifdef USE_MPI
	///	<summary>
	///		MPI_Request objects used for asynchronous exchange of data.
	///	</summary>
	MPI_Request m_reqSend;
	MPI_Request m_reqRecv;
#endif

	///	<summary>
	///		DataArray1D used for exchange of data and fluxes.
	///	</summary>
	DataArray1D<double> m_vecSendBuffer;
	DataArray1D<double> m_vecRecvBuffer;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Neighbor objects which are along GridPatch boundaries.
///	</summary>
class ExteriorNeighbor : public Neighbor {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	ExteriorNeighbor(
		const ExchangeBufferInfo & info
	) :
		Neighbor(info)
	{ }

public:
	///	<summary>
	///		Set a value in the SendBuffer directly.
	///	</summary>
	inline void SetSendBuffer(
		int iC,
		int iK,
		int iA,
		double dValue
	) {
		int ix = (iC * (m_info.sMaxRElements-1) + iK) * m_info.sBoundarySize;

		if (m_info.fReverseDirection) {
			ix += m_info.ixSecond - iA - 1;
		} else {
			ix += iA - m_info.ixFirst;
		}

		m_vecSendBuffer[ix] = dValue;
	}

	///	<summary>
	///		Get a value in the RecvBuffer directly.
	///	</summary>
	inline double GetRecvBuffer(
		int iC,
		int iK,
		int iA
	) const {
		int ix =
			(iC * (m_info.sMaxRElements-1) + iK) * m_info.sBoundarySize
				+ (iA - m_info.ixFirst);

		return m_vecRecvBuffer[ix];
	}

public:
	///	<summary>
	///		Prepare an asynchronous receive.
	///	</summary>
	virtual void PrepareExchange(
		const Grid & grid
	);

	///	<summary>
	///		Pack data into the send buffer.
	///	</summary>
	virtual void Pack(
		const DataArray3D<double> & data
	);

	///	<summary>
	///		Pack data into the send buffer.
	///	</summary>
	virtual void Pack(
		const Grid & grid,
		const DataArray4D<double> & data
	);

	///	<summary>
	///		Send data to other processors.
	///	</summary>
	virtual void Send(
		const Grid & grid
	);

	///	<summary>
	///		Unpack data from the receive buffer.
	///	</summary>
	void Unpack(
		DataArray3D<double> & data
	);

	///	<summary>
	///		Unpack data from the receive buffer.
	///	</summary>
	void Unpack(
		const Grid & grid,
		DataArray4D<double> & data
	);

	///	<summary>
	///		Wait for the send request to complete.
	///	</summary>
	virtual void WaitSend();
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Nested neighbors of the GridPatch.
///	</summary>
class InteriorNeighbor : public Neighbor {
	///	<summary>
	///		Send whole patch to coarse levels.
	///	</summary>

public:
	///	<summary>
	///		Prepare an asynchronous receive.
	///	</summary>
	virtual void PrepareExchange() {
	}

	///	<summary>
	///		Pack data into the send buffer.
	///	</summary>
	virtual void Pack(
		const Grid & grid,
		const DataArray3D<double> & data
	) {
	}

	///	<summary>
	///		Pack data into the send buffer.
	///	</summary>
	virtual void Pack(
		const Grid & grid,
		const DataArray4D<double> & data
	) {
	}

	///	<summary>
	///		Send data to other processors.
	///	</summary>
	virtual void Send(const Grid & grid) {
	}

	///	<summary>
	///		Unpack data from the receive buffer.
	///	</summary>
	void Unpack(
		const Grid & grid,
		DataArray3D<double> & data
	) {
	}

	///	<summary>
	///		Unpack data from the receive buffer.
	///	</summary>
	virtual void Unpack(
		const Grid & grid,
		const DataArray4D<double> & data
	) {
	}

	///	<summary>
	///		Wait for the send request to complete.
	///	</summary>
	virtual void WaitSend() {
	}
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Connectivity between GridPatch objects.
///	</summary>
class Connectivity {

public:
	///	<summary>
	///		Vector of exterior boundary neighbors.
	///	</summary>
	typedef std::vector<ExteriorNeighbor *> ExteriorNeighborVector;

	///	<summary>
	///		Vector of interior boundary neighbors.
	///	</summary>
	typedef std::vector<InteriorNeighbor *> InteriorNeighborVector;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	Connectivity(GridPatch & patch) :
		m_patch(patch)
	{ }

	///	<summary>
	///		Destructor.
	///	</summary>
	virtual ~Connectivity();

public:
	///	<summary>
	///		Clear all neighbors from this connectivity object.
	///	</summary>
	void ClearNeighbors();

	///	<summary>
	///		Add a new exterior neighbor.
	///	</summary>
	void AddExteriorNeighbor(
		ExteriorNeighbor * pNeighbor
	) {
		m_vecExteriorNeighbors.push_back(pNeighbor);
	}

public:
	///	<summary>
	///		Done adding new exterior neighbors.  Set up flux connectivity.
	///	</summary>
	void BuildFluxConnectivity();

public:
	///	<summary>
	///		Prepare for the exchange of data between processors.
	///	</summary>
	void PrepareExchange();

	///	<summary>
	///		Pack data into the send buffer in preparation for a send.
	///	</summary>
	void Pack(
		const DataArray3D<double> & data
	);

	///	<summary>
	///		Pack data into the send buffer in preparation for a send.
	///	</summary>
	void Pack(
		const DataArray4D<double> & data
	);

	///	<summary>
	///		Send data to other processors.
	///	</summary>
	void Send();

	///	<summary>
	///		Send buffers to other processors.
	///	</summary>
	void SendBuffers();

	///	<summary>
	///		Wait for a neighbor to receive a message and return a pointer
	///		to the neighbor when one has been received.
	///	</summary>
	Neighbor * WaitReceive();

	///	<summary>
	///		Wait for all asynchronous send requests to complete.
	///	</summary>
	void WaitSend();

public:
	///	<summary>
	///		Verify valid connectivity.
	///	</summary>
	void Verify() const {
		_EXCEPTIONT("Not implemented.");
	}

public:
	///	<summary>
	///		Get the reference to the grid patch.
	///	</summary>
	const GridPatch & GetGridPatch() const {
		return m_patch;
	}
/*
	///	<summary>
	///		Get the vector of exterior boundary neighbors.
	///	</summary>
	const ExteriorNeighborVector & GetExteriorNeighbors() const {
		return m_vecExteriorNeighbors;
	}

	///	<summary>
	///		Get the vector of interior boundary neighbors.
	///	</summary>
	const InteriorNeighborVector & GetInteriorNeighbors() const {
		return m_vecNestedNeighbors;
	}
*/
	///	<summary>
	///		Get the expected number of messages received during an exchange.
	///	</summary>
	int GetExpectedMessageCount() const;

	///	<summary>
	///		Set a value directly in send buffer (one halo element).
	///	</summary>
	///	<param name="iA">
	///		Interior index along edge (index 0 corresponds to InteriorBegin)
	///	</param>
	inline void SetSendBuffer(
		Direction dir,
		int iC,
		int iK,
		int iA,
		double dValue
	) {
		_EXCEPTIONT("Disabled");
		//m_vecExteriorEdge[(int)(dir)][iA]->SetSendBuffer(iC, iK, iA, dValue);
	}

	///	<summary>
	///		Get a value directly in recv buffer (one halo element).
	///	</summary>
	///	<param name="iA">
	///		Interior index along edge (index 0 corresponds to InteriorBegin)
	///	</param>
	inline double GetRecvBuffer(
		Direction dir,
		int iC,
		int iK,
		int iA
	) {
		_EXCEPTIONT("Disabled");
		//return m_vecExteriorEdge[(int)(dir)][iA]->GetRecvBuffer(iC, iK, iA);
	}

	///	<summary>
	///		Determine if the direction of increasing coordinate is flipped
	///		across the edge.
	///	</summary>
	inline bool IsCoordinateFlipped(
		Direction dir,
		int iA
	) const {
		_EXCEPTIONT("Disabled");
		//return m_vecExteriorEdge[(int)(dir)][iA]->m_fFlippedCoordinate;
	}

private:
	///	<summary>
	///		Reference to parent object.
	///	</summary>
	GridPatch & m_patch;

	///	<summary>
	///		Vector of exterior boundary neighbors.
	///	</summary>
	ExteriorNeighborVector m_vecExteriorNeighbors;

	///	<summary>
	///		Vector of interior boundary indices.
	///	</summary>
	////InteriorNeighborVector m_vecNestedNeighbors;

	///	<summary>
	///		Array indicating the panel in the specified direction.
	///	</summary>
	std::vector<int> m_ixPanelDirection;
/*
	///	<summary>
	///		Pointer to exterior neighbors along edges, by index.
	///	</summary>
	std::vector<ExteriorNeighbor *> m_vecExteriorEdge[8];
*/
};

///////////////////////////////////////////////////////////////////////////////

#endif

