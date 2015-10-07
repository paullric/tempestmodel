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

#define USE_MPI
//#define USE_HPX
//#define USE_OCR

///////////////////////////////////////////////////////////////////////////////

#include "Direction.h"

#include "DataArray1D.h"
#include "DataArray3D.h"
#include "DataArray4D.h"

#ifdef USE_MPI
#include "mpi.h"
#endif

#include <map>
#include <vector>

///////////////////////////////////////////////////////////////////////////////

class Grid;
class GridPatch;
class Neighbor;
class Connectivity;

///////////////////////////////////////////////////////////////////////////////

class ExchangeBufferKey {

public:
	///	<summary>
	///		Default constructor.
	///	</summary>
	ExchangeBufferKey() :
		ixSourcePatch(-1),
		ixTargetPatch(-1),
		dir(Direction_Middle)
	{ }

	///	<summary>
	///		Constructor.
	///	</summary>
	ExchangeBufferKey(
		int a_ixSourcePatch,
		int a_ixTargetPatch,
		Direction a_dir
	) :
		ixSourcePatch(a_ixSourcePatch),
		ixTargetPatch(a_ixTargetPatch),
		dir(a_dir)
	{ }

public:
	///	<summary>
	///		Comparator.
	///	</summary>
	bool operator<(const ExchangeBufferKey & key) const {
		if (ixSourcePatch < key.ixSourcePatch) {
			return true;
		}
		if (ixSourcePatch > key.ixSourcePatch) {
			return false;
		}
		if (ixTargetPatch < key.ixTargetPatch) {
			return true;
		}
		if (ixTargetPatch > key.ixTargetPatch) {
			return false;
		}
		if (dir < key.dir) {
			return true;
		}
		if (dir > key.dir) {
			return false;
		}
		return false;
	}

public:
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
};

///////////////////////////////////////////////////////////////////////////////

class ExchangeBufferMeta {

public:
	///	<summary>
	///		Initial coordinate index.
	///	</summary>
	int ixFirst;

	///	<summary>
	///		Final coordinate index.
	///	</summary>
	int ixSecond;

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
	///		Total data size.
	///	</summary>
	size_t sDataSize;
};

///////////////////////////////////////////////////////////////////////////////

class ExchangeBufferRegistry {

public:
	///	<summary>
	///		Map from ExchangeBufferKey to unique identifier for ExchangeBuffer.
	///	</summary>
	typedef std::map<ExchangeBufferKey, size_t> MapExchangeBuffer;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	ExchangeBufferRegistry() {
	}

public:
	///	<summary>
	///		Register a new ExchangeBuffer.
	///	</summary>
	void Register(
		const ExchangeBufferKey & key,
		const ExchangeBufferMeta & meta
	) {
		size_t sNextEntry = m_vecMeta.size();

		m_mapKeyToIndex.insert(
			MapExchangeBuffer::value_type(key, sNextEntry));

		m_vecMeta.push_back(meta);
	}

public:
	///	<summary>
	///		Map of ExchangeBufferKeys to vector indices.
	///	</summary>
	MapExchangeBuffer m_mapKeyToIndex;

	///	<summary>
	///		Metadata for each exchange buffer.
	///	</summary>
	std::vector<ExchangeBufferMeta> m_vecMeta;

#ifdef USE_MPI
	///	<summary>
	///		Data pointers for exchange buffer.
	///	</summary>
	std::vector<unsigned char *> m_vecData;
#endif
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
		int ixPatch,
		int ixNeighbor,
		bool fReverseDirection,
		bool fFlippedCoordinate,
		size_t sBoundarySize
	) : 
		m_ixPatch(ixPatch),
		m_ixNeighbor(ixNeighbor),
		m_fReverseDirection(fReverseDirection),
		m_fFlippedCoordinate(fFlippedCoordinate),
		m_sBoundarySize(sBoundarySize),
		m_sMaxRElements(0),
		m_sHaloElements(0),
		m_sComponents(0),
		m_fComplete(false),
		m_ixSendBuffer(0),
		m_ixRecvBuffer(0)
	{ }

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
		size_t sMaxRElements,
		size_t sHaloElements,
		size_t sComponents
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
		const Grid & grid,
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
			  m_sBoundarySize
			* m_sMaxRElements
			* m_sHaloElements
			* m_sComponents;
	}

	///	<summary>
	///		Send data to other processors.
	///	</summary>
	virtual void Send(const Grid & grid) = 0;

	///	<summary>
	///		Receive data from other processors.
	///	</summary>
	virtual void Unpack(
		const Grid & grid,
		DataArray3D<double> & data
	) = 0;

	///	<summary>
	///		Receive data from other processors.
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
	///		Index of the source patch.
	///	</summary>
	int m_ixPatch;

	///	<summary>
	///		Index of the target patch.
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
	size_t m_sBoundarySize;

	///	<summary>
	///		Number of radial elements.
	///	</summary>
	size_t m_sMaxRElements;

	///	<summary>
	///		Number of halo elements.
	///	</summary>
	size_t m_sHaloElements;

	///	<summary>
	///		Number of variables.
	///	</summary>
	size_t m_sComponents;

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
		Direction dir,
		Direction dirOpposite,
		int ixPatch,
		int ixNeighbor,
		bool fReverseDirection,
		bool fFlippedCoordinate,
		size_t sBoundarySize,
		int ixFirst,
		int ixSecond
	) :
		m_dir(dir),
		m_dirOpposite(dirOpposite),
		Neighbor(
			ixPatch,
			ixNeighbor,
			fReverseDirection,
			fFlippedCoordinate,
			sBoundarySize),
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
		int ix = (iC * (m_sMaxRElements-1) + iK) * m_sBoundarySize;

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
			(iC * (m_sMaxRElements-1) + iK) * m_sBoundarySize
				+ (iA - m_ixFirst);

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
		const Grid & grid,
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
		const Grid & grid,
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

