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

#include "Direction.h"

#include "DataVector.h"
#include "GridData4D.h"

#include "mpi.h"
#include <vector>

///////////////////////////////////////////////////////////////////////////////

class Connectivity;
class GridPatch;

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
		const Connectivity * pConnect,
		int ixNeighbor,
		bool fReverseDirection,
		bool fFlippedCoordinate,
		int nBoundarySize
	) : 
		m_pConnect(pConnect),
		m_ixNeighbor(ixNeighbor),
		m_fReverseDirection(fReverseDirection),
		m_fFlippedCoordinate(fFlippedCoordinate),
		m_nBoundarySize(nBoundarySize),
		m_nMaxRElements(0),
		m_nHaloElements(0),
		m_nComponents(0),
		m_fComplete(false),
		m_ixSendBuffer(0),
		m_ixRecvBuffer(0)
	{
		if (m_nBoundarySize < 0) {
			_EXCEPTIONT("Invalid boundary size");
		}
	}

	///	<summary>
	///		Virtual destructor.
	///	</summary>
	virtual ~Neighbor() {
	}

public:
	///	<summary>
	///		Initialize the exchange buffers.
	///	</summary>
	void InitializeBuffers(
		int nRElements,
		int nHaloElements,
		int nComponents
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
	virtual void PrepareExchange() {

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
		const GridData3D & data
	) = 0;

	///	<summary>
	///		Pack data into the send buffer.
	///	</summary>
	virtual void Pack(
		const GridData4D & data
	) = 0;

	///	<summary>
	///		Reset the send buffer size.
	///	</summary>
	void ResetSendBufferSize() {
		m_ixSendBuffer =
			  m_nBoundarySize
			* m_nMaxRElements
			* m_nHaloElements
			* m_nComponents;
	}

	///	<summary>
	///		Send data to other processors.
	///	</summary>
	virtual void Send() = 0;

	///	<summary>
	///		Receive data from other processors.
	///	</summary>
	virtual void Unpack(
		GridData3D & data
	) = 0;

	///	<summary>
	///		Receive data from other processors.
	///	</summary>
	virtual void Unpack(
		GridData4D & data
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
	///		Parent Connectivity object.
	///	</summary>
	const Connectivity * m_pConnect;

	///	<summary>
	///		Index of the neighboring patch.
	///	</summary>
	int m_ixNeighbor;

	///	<summary>
	///		Indicator of whether or not directions are flipped across the edge.
	///	</summary>
	bool m_fReverseDirection;

	///	<summary>
	///		Indicator of whether or not the direction of increasing coordinate
	///		indices flips across the edge.
	///	</summary>
	bool m_fFlippedCoordinate;

	///	<summary>
	///		Number of independent grid cells along boundary.
	///	</summary>
	int m_nBoundarySize;

	///	<summary>
	///		Number of radial elements.
	///	</summary>
	int m_nMaxRElements;

	///	<summary>
	///		Number of halo elements.
	///	</summary>
	int m_nHaloElements;

	///	<summary>
	///		Number of variables.
	///	</summary>
	int m_nComponents;

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

	///	<summary>
	///		MPI_Request objects used for asynchronous exchange of data.
	///	</summary>
	MPI_Request m_reqSend;
	MPI_Request m_reqRecv;

	///	<summary>
	///		DataVector used for exchange of data and fluxes.
	///	</summary>
	DataVector<double> m_vecSendBuffer;
	DataVector<double> m_vecRecvBuffer;
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
		const Connectivity * pConnect,
		Direction dir,
		Direction dirOpposite,
		int ixNeighbor,
		bool fReverseDirection,
		bool fFlippedCoordinate,
		int nBoundarySize,
		int ixFirst,
		int ixSecond
	) :
		m_dir(dir),
		m_dirOpposite(dirOpposite),
		Neighbor(
			pConnect,
			ixNeighbor,
			fReverseDirection,
			fFlippedCoordinate,
			nBoundarySize),
		m_ixFirst(ixFirst),
		m_ixSecond(ixSecond)
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
		int ix = (iC * m_nMaxRElements + iK) * m_nBoundarySize;

		if (m_fReverseDirection) {
			ix += m_ixSecond - iA - 1;
		} else {
			ix += iA - m_ixFirst;
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
			(iC * m_nMaxRElements + iK) * m_nBoundarySize
				+ (iA - m_ixFirst);

		return m_vecRecvBuffer[ix];
	}

public:
	///	<summary>
	///		Prepare an asynchronous receive.
	///	</summary>
	virtual void PrepareExchange();

	///	<summary>
	///		Pack data into the send buffer.
	///	</summary>
	virtual void Pack(
		const GridData3D & data
	);

	///	<summary>
	///		Pack data into the send buffer.
	///	</summary>
	virtual void Pack(
		const GridData4D & data
	);

	///	<summary>
	///		Send data to other processors.
	///	</summary>
	virtual void Send();

	///	<summary>
	///		Unpack data from the receive buffer.
	///	</summary>
	void Unpack(
		GridData3D & data
	);

	///	<summary>
	///		Unpack data from the receive buffer.
	///	</summary>
	void Unpack(
		GridData4D & data
	);

	///	<summary>
	///		Wait for the send request to complete.
	///	</summary>
	virtual void WaitSend();

public:
	///	<summary>
	///		Direction on this block towards neighbor block.
	///	</summary>
	Direction m_dir;

	///	<summary>
	///		Directon on neighbor block towards this block.
	///	</summary>
	Direction m_dirOpposite;

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
		const GridData3D & data
	) {
	}

	///	<summary>
	///		Pack data into the send buffer.
	///	</summary>
	virtual void Pack(
		const GridData4D & data
	) {
	}

	///	<summary>
	///		Send data to other processors.
	///	</summary>
	virtual void Send() {
	}

	///	<summary>
	///		Unpack data from the receive buffer.
	///	</summary>
	void Unpack(
		GridData3D & data
	) {
	}

	///	<summary>
	///		Unpack data from the receive buffer.
	///	</summary>
	virtual void Unpack(
		const GridData4D & data
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
		const GridData3D & data
	);

	///	<summary>
	///		Pack data into the send buffer in preparation for a send.
	///	</summary>
	void Pack(
		const GridData4D & data
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
		m_vecExteriorEdge[(int)(dir)][iA]->SetSendBuffer(iC, iK, iA, dValue);
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
		return m_vecExteriorEdge[(int)(dir)][iA]->GetRecvBuffer(iC, iK, iA);
	}

	///	<summary>
	///		Determine if the direction of increasing coordinate is flipped
	///		across the edge.
	///	</summary>
	inline bool IsCoordinateFlipped(
		Direction dir,
		int iA
	) const {
		return m_vecExteriorEdge[(int)(dir)][iA]->m_fFlippedCoordinate;
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

	///	<summary>
	///		Pointer to exterior neighbors along edges, by index.
	///	</summary>
	std::vector<ExteriorNeighbor *> m_vecExteriorEdge[8];
};

///////////////////////////////////////////////////////////////////////////////

#endif

