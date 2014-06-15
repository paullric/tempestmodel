///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridPatch.h
///	\author  Paul Ullrich
///	\version February 25, 2013
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

#ifndef _GRIDPATCH_H_
#define _GRIDPATCH_H_

#include "DataVector.h"
#include "DataMatrix.h"
#include "DataMatrix3D.h"
#include "GridData3D.h"
#include "GridData4D.h"
#include "FluxBuffer.h"
#include "PatchBox.h"
#include "Connectivity.h"
#include "ChecksumType.h"

///////////////////////////////////////////////////////////////////////////////

class Time;
class Grid;

class TestCase;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class for storing a rectangular patch of the grid and all metric
///		terms relevant for integration on this patch.
///	</summary>
class GridPatch {

friend class Grid;

public:
	///	<summary>
	///		A vector containing panel indices.
	///	</summary>
	typedef std::vector<int> PanelIndexVector;

public:
	///	<summary>
	///		The PatchIndex indicating an invalid patch.
	///	</summary>
	static const int InvalidIndex = (-1);

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	GridPatch(
		Grid & grid,
		int ixPatch,
		const PatchBox & box
	);

	///	<summary>
	///		Virtual destructor.
	///	</summary>
	virtual ~GridPatch()
	{ }

public:
	///	<summary>
	///		Total number of nodes in 2D (ignoring vertical nodes or staggering)
	///	</summary>
	int GetTotalNodeCount2D() const;

	///	<summary>
	///		Total number of nodes in data.
	///	</summary>
	int GetTotalNodeCount(
		DataLocation loc = DataLocation_Node
	) const;

	///	<summary>
	///		Total number of degrees of freedom in data.
	///	</summary>
	int GetTotalDegreesOfFreedom(
		DataType eDataType,
		DataLocation eDataLocation = DataLocation_None
	) const;

public:
	///	<summary>
	///		Initialize patch as remote (data on an exterior processor).
	///	</summary>
	virtual void InitializeDataRemote(int iProcessor);

	///	<summary>
	///		Initialize patch as local (data on this processor).
	///	</summary>
	virtual void InitializeDataLocal();

	///	<summary>
	///		Delete data and reset patch to a stub.
	///	</summary>
	virtual void DeinitializeData();

	///	<summary>
	///		Initialize geometric terms.
	///	</summary>
	virtual void EvaluateGeometricTerms() = 0;

	///	<summary>
	///		Initialize state and tracer data from a TestCase.  Also adjust
	///		geometric quantities that are dependent on the TestCase.
	///	</summary>
	virtual void EvaluateTestCase(
		const TestCase & test,
		const Time & time,
		int iDataIndex = 0
	) = 0;

	///	<summary>
	///		Initialize state and tracer data from a TestCase.
	///	</summary>
	virtual void EvaluateTestCase_StateOnly(
		const TestCase & test,
		const Time & time,
		int iDataIndex = 0
	) {
		_EXCEPTIONT("Unimplemented.");
	}

public:
	///	<summary>
	///		Connect one patch to a second along entire edge.
	///	</summary>
	void ExteriorConnect(
		Direction dirFirst,
		const GridPatch * pPatchSecond
	);

	///	<summary>
	///		Add an exterior connection to another patch.
	///	</summary>
	///	<param name="ixFirst">
	///		First index (begin node for Right, Top, Left, Bottom and alpha
	///		coordinate for TopRight, TopLeft, BottomLeft, BottomRight)
	///	</param>
	///	<param name="ixSecond">
	///		Second index (end node for Right, Top, Left, Bottom and beta
	///		coordinate for TopRight, TopLeft, BottomLeft, BottomRight)
	///	</param>
	void ExteriorConnect(
		Direction dirFirst,
		const GridPatch * pPatchSecond,
		int ixFirst,
		int ixSecond
	);

public:
	///	<summary>
	///		Compute the radial component of the curl on the grid given two
	///		contravariant vector fields.
	///	</summary>
	virtual void ComputeCurlAndDiv(
		const GridData3D & dataUa,
		const GridData3D & dataUb
	) const {
		_EXCEPTIONT("Unimplemented");
	}

	///	<summary>
	///		Compute vorticity on the GridPatch.
	///	</summary>
	virtual void ComputeVorticityDivergence(
		int iDataIndex
	) {
		_EXCEPTIONT("Unimplemented");
	}

	///	<summary>
	///		Compute temperature on the GridPatch.
	///	</summary>
	virtual void ComputeTemperature(
		int iDataIndex,
		DataLocation loc = DataLocation_Node
	);

public:
	///	<summary>
	///		Add local masses to checksum total.
	///	</summary>
	void Checksum(
		DataType eDataType,
		DataVector<double> & dChecksums,
		int iDataIndex,
		ChecksumType eChecksumType
	) const;

	///	<summary>
	///		Compute total energy over this GridPatch.
	///	</summary>
	double ComputeTotalEnergy(
		int iDataIndex
	) const;

	///	<summary>
	///		Compute total potential enstrophy over this GridPatch.
	///	</summary>
	double ComputeTotalPotentialEnstrophy(
		int iDataIndex
	);

public:
	///	<summary>
	///		Prepare for the exchange of halo data between processors.
	///	</summary>
	void PrepareExchange();

	///	<summary>
	///		Send halo data to other processors.
	///	</summary>
	void Send(
		DataType eDataType,
		int iDataIndex
	);

	///	<summary>
	///		Receive halo data from other processors.
	///	</summary>
	void Receive(
		DataType eDataType,
		int iDataIndex
	);

	///	<summary>
	///		Send buffers to other processors.
	///	</summary>
	void SendBuffers();

	///	<summary>
	///		Receive buffers from other processors.
	///	</summary>
	void ReceiveBuffers();

	///	<summary>
	///		Complete the exchange of data between processors.
	///	</summary>
	void CompleteExchange();

public:
	///	<summary>
	///		Copy data from one data index to another.
	///	</summary>
	void CopyData(
		int ixSource,
		int ixDest,
		DataType eDataType
	);

	///	<summary>
	///		Compute a linear combination of data and store at specified index.
	///	</summary>
	void LinearCombineData(
		const DataVector<double> & dCoeff,
		int ixDest,
		DataType eDataType
	);

	///	<summary>
	///		Zero the data at the specified index.
	///	</summary>
	void ZeroData(
		int ixData,
		DataType eDataType
	);

	///	<summary>
	///		Add the reference state to the specified state data index.
	///	</summary>
	void AddReferenceState(
		int ix
	);

public:
	///	<summary>
	///		Interpolate data vertically from Nodes to REdges.
	///	</summary>
	virtual void InterpolateNodeToREdge(
		int iVar,
		int iDataIndex
	);

	///	<summary>
	///		Interpolate data vertically from REdges to Nodes.
	///	</summary>
	virtual void InterpolateREdgeToNode(
		int iVar,
		int iDataIndex
	);

public:
	///	<summary>
	///		Compute velocity across xi surfaces (xi_dot) at nodes.
	///	</summary>
	inline double CalculateXiDotNode(
		int k,
		int iA,
		int iB,
		double dUa,
		double dUb,
		double dUr
	) {
		return
			( m_dataOrthonormNode[k][iA][iB][0] * dUa
			+ m_dataOrthonormNode[k][iA][iB][1] * dUb
			+ m_dataOrthonormNode[k][iA][iB][2] * dUr);
	}

	///	<summary>
	///		Compute velocity across xi surfaces (xi_dot) at interfaces.
	///	</summary>
	inline double CalculateXiDotREdge(
		int k,
		int iA,
		int iB,
		double dUa,
		double dUb,
		double dUr
	) {
		return
			( m_dataOrthonormREdge[k][iA][iB][0] * dUa
			+ m_dataOrthonormREdge[k][iA][iB][1] * dUb
			+ m_dataOrthonormREdge[k][iA][iB][2] * dUr);
	}

	///	<summary>
	///		Compute vertical velocity needed to enforce no flow condition
	///		(xi_dot = 0) at interfaces.
	///	</summary>
	inline double CalculateNoFlowUrREdge(
		int k,
		int iA,
		int iB,
		double dUa,
		double dUb
	) {
		return - 1.0 / m_dataOrthonormREdge[k][iA][iB][2] * (
			  m_dataOrthonormREdge[k][iA][iB][0] * dUa
			+ m_dataOrthonormREdge[k][iA][iB][1] * dUb);
	}

public:
	///	<summary>
	///		Linearly interpolate data horizontally to the specified points.
	///	</summary>
	virtual void InterpolateData(
		const DataVector<double> & dAlpha,
		const DataVector<double> & dBeta,
		const DataVector<int> & iPanel,
		DataType eDataType,
		DataLocation eDataLocation,
		bool fInterpAllVariables,
		DataMatrix3D<double> & dInterpData,
		bool fIncludeReferenceState = true,
		bool fConvertToPrimitive = false
	);

public:
	///	<summary>
	///		Get the reference to the parent grid.
	///	</summary>
	const Grid & GetGrid() const {
		return m_grid;
	}

	///	<summary>
	///		Get the global index of this patch.
	///	</summary>
	int GetPatchIndex() const {
		return m_ixPatch;
	}

	///	<summary>
	///		Get the processor id of this patch.
	///	</summary>
	int GetProcessor() const {
		return m_iProcessor;
	}

	///	<summary>
	///		Get the PatchBox defining the location of this patch on the grid.
	///	</summary>
	const PatchBox & GetPatchBox() const {
		return m_box;
	}

	///	<summary>
	///		Get the specified neighboring panel index.
	///	</summary>
	int GetNeighborPanel(Direction dir) const {
		return m_ixNeighborPanel[(int)(dir)];
	}

	///	<summary>
	///		Get the connectivity of this patch.
	///	</summary>
	const Connectivity & GetConnectivity() const {
		return m_connect;
	}

	///	<summary>
	///		Get the connectivity of this patch.
	///	</summary>
	Connectivity & GetConnectivity() {
		return m_connect;
	}

	///	<summary>
	///		Returns true if this GridPatch contains data, ie. is not a stub.
	///	</summary>
	bool ContainsData() const {
		return m_fContainsData;
	}

	///	<summary>
	///		Get the 2D Jacobian matrix.
	///	</summary>
	const DataMatrix<double> & GetJacobian2D() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataJacobian2D;
	}

	///	<summary>
	///		Get the  components of the contravariant metric (alpha)
	///	</summary>
	const DataMatrix3D<double> & GetContraMetric2DA() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataContraMetric2DA;
	}

	///	<summary>
	///		Get the  components of the contravariant metric (beta)
	///	</summary>
	const DataMatrix3D<double> & GetContraMetric2DB() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataContraMetric2DB;
	}

	///	<summary>
	///		Get the nodal Jacobian matrix.
	///	</summary>
	const DataMatrix3D<double> & GetJacobian() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataJacobian;
	}

	///	<summary>
	///		Get the interface Jacobian matrix.
	///	</summary>
	const DataMatrix3D<double> & GetJacobianREdge() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataJacobianREdge;
	}

	///	<summary>
	///		Get the  components of the contravariant metric (alpha)
	///	</summary>
	const DataMatrix4D<double> & GetContraMetricA() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataContraMetricA;
	}

	///	<summary>
	///		Get the  components of the contravariant metric (beta)
	///	</summary>
	const DataMatrix4D<double> & GetContraMetricB() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataContraMetricB;
	}

	///	<summary>
	///		Get the  components of the contravariant metric (xi)
	///	</summary>
	const DataMatrix4D<double> & GetContraMetricXi() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataContraMetricXi;
	}

	///	<summary>
	///		Get the nodal components of the Christoffel alpha symbols.
	///	</summary>
	const DataMatrix3D<double> & GetChristoffelA() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataChristoffelA;
	}

	///	<summary>
	///		Get the nodal components of the Christoffel beta symbols.
	///	</summary>
	const DataMatrix3D<double> & GetChristoffelB() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataChristoffelB;
	}

	///	<summary>
	///		Get the nodal components of the Christoffel xi symbols.
	///	</summary>
	const DataMatrix4D<double> & GetChristoffelXi() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataChristoffelXi;
	}
 
 	///	<summary>
	///		Get the orthonormalization components at nodes.
	///	</summary>
	const DataMatrix4D<double> & GetOrthonormNode() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataOrthonormNode;
	}

	///	<summary>
	///		Get the orthonormalization components at interfaces.
	///	</summary>
	const DataMatrix4D<double> & GetOrthonormREdge() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataOrthonormREdge;
	}

	///	<summary>
	///		Get the nodal element area matrix.
	///	</summary>
	const DataMatrix3D<double> & GetElementArea() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataElementArea;
	}

	///	<summary>
	///		Get the interface element area matrix.
	///	</summary>
	const DataMatrix3D<double> & GetElementAreaREdge() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataElementAreaREdge;
	}

	///	<summary>
	///		Get the nodal topography matrix.
	///	</summary>
	DataMatrix<double> & GetTopography() {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataTopography;
	}

	///	<summary>
	///		Get the nodal topography matrix.
	///	</summary>
	const DataMatrix<double> & GetTopography() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataTopography;
	}

	///	<summary>
	///		Get the nodal longitude matrix.
	///	</summary>
	const DataMatrix<double> & GetLongitude() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataLon;
	}

	///	<summary>
	///		Get the nodal latitude matrix.
	///	</summary>
	const DataMatrix<double> & GetLatitude() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataLat;
	}

	///	<summary>
	///		Get the nodal Coriolis parameter.
	///	</summary>
	const DataMatrix<double> & GetCoriolisF() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataCoriolisF;
	}

	///	<summary>
	///		Get the radial coordinate matrix on model levels.
	///	</summary>
	const DataMatrix3D<double> & GetZLevels() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataZLevels;
	}

	///	<summary>
	///		Get the radial coordinate matrix on model interfaces.
	///	</summary>
	const DataMatrix3D<double> & GetZInterfaces() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataZInterfaces;
	}

	///	<summary>
	///		Get the reference state.
	///	</summary>
	GridData4D & GetReferenceState(
		DataLocation loc = DataLocation_Node
	) {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}
		if (loc == DataLocation_Node) {
			return m_dataRefStateNode;
		} else if (loc == DataLocation_REdge) {
			return m_dataRefStateREdge;
		} else {
			_EXCEPTIONT("Invalid DataLocation");
		}
	}

	///	<summary>
	///		Get the reference state.
	///	</summary>
	const GridData4D & GetReferenceState(
		DataLocation loc = DataLocation_Node
	) const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}
		if (loc == DataLocation_Node) {
			return m_dataRefStateNode;
		} else if (loc == DataLocation_REdge) {
			return m_dataRefStateREdge;
		} else {
			_EXCEPTIONT("Invalid DataLocation");
		}
	}

	///	<summary>
	///		Get the state data matrix with the specified index.
	///	</summary>
	GridData4D & GetDataState(
		int ix,
		DataLocation loc = DataLocation_Node
	) {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}
		if (loc == DataLocation_Node) {
			if ((ix < 0) || (ix > m_datavecStateNode.size())) {
				_EXCEPTIONT("Invalid index in StateData vector.");
			}
			return m_datavecStateNode[ix];

		} else if (loc == DataLocation_REdge) {
			if ((ix < 0) || (ix > m_datavecStateREdge.size())) {
				_EXCEPTIONT("Invalid index in StateData vector.");
			}
			return m_datavecStateREdge[ix];

		} else {
			_EXCEPTIONT("Invalid DataLocation");
		}
	}

	///	<summary>
	///		Get the state data matrix with the specified index.
	///	</summary>
	const GridData4D & GetDataState(
		int ix,
		DataLocation loc = DataLocation_Node
	) const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}
		if (loc == DataLocation_Node) {
			if ((ix < 0) || (ix > m_datavecStateNode.size())) {
				_EXCEPTIONT("Invalid index in StateData vector.");
			}
			return m_datavecStateNode[ix];

		} else if (loc == DataLocation_REdge) {
			if ((ix < 0) || (ix > m_datavecStateREdge.size())) {
				_EXCEPTIONT("Invalid index in StateData vector.");
			}
			return m_datavecStateREdge[ix];

		} else {
			_EXCEPTIONT("Invalid DataLocation");
		}
	}

	///	<summary>
	///		Get the tracer data matrix with the specified index.
	///	</summary>
	GridData4D & GetDataTracers(int ix) {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}
		if ((ix < 0) || (ix > m_datavecTracers.size())) {
			_EXCEPTIONT("Invalid index in TracersData vector.");
		}
		return m_datavecTracers[ix];
	}

	///	<summary>
	///		Get the tracer data matrix with the specified index.
	///	</summary>
	const GridData4D & GetDataTracers(int ix) const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}
		if ((ix < 0) || (ix > m_datavecTracers.size())) {
			_EXCEPTIONT("Invalid index in TracersData vector.");
		}
		return m_datavecTracers[ix];
	}

	///	<summary>
	///		Get the auxiliary grid data with the specified index.
	///	</summary>
	GridData3D & GetHorizontalDynamicsAuxData(
		int ix,
		DataLocation loc = DataLocation_Node
	) {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}
		if ((ix < 0) || (ix > m_datavecAuxNode[0].size())) {
			_EXCEPTIONT("Invalid index in TracersData vector.");
		}
		if (loc == DataLocation_Node) {
			return m_datavecAuxNode[0][ix];

		} else if (loc == DataLocation_REdge) {
			return m_datavecAuxREdge[0][ix];

		} else {
			_EXCEPTIONT("Invalid DataLocation.");
		}
	}

	///	<summary>
	///		Get the auxiliary grid data with the specified index.
	///	</summary>
	const GridData3D & GetHorizontalDynamicsAuxData(
		int ix,
		DataLocation loc = DataLocation_Node
	) const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}
		if ((ix < 0) || (ix > m_datavecAuxNode[0].size())) {
			_EXCEPTIONT("Invalid index in TracersData vector.");
		}
		if (loc == DataLocation_Node) {
			return m_datavecAuxNode[0][ix];

		} else if (loc == DataLocation_REdge) {
			return m_datavecAuxREdge[0][ix];

		} else {
			_EXCEPTIONT("Invalid DataLocation.");
		}
	}

	///	<summary>
	///		Get the auxiliary grid data with the specified index.
	///	</summary>
	GridData3D & GetVerticalDynamicsAuxData(
		int ix,
		DataLocation loc = DataLocation_Node
	) {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}
		if ((ix < 0) || (ix > m_datavecAuxNode[1].size())) {
			_EXCEPTIONT("Invalid index in TracersData vector.");
		}
		if (loc == DataLocation_Node) {
			return m_datavecAuxNode[1][ix];

		} else if (loc == DataLocation_REdge) {
			return m_datavecAuxREdge[1][ix];

		} else {
			_EXCEPTIONT("Invalid DataLocation.");
		}
	}

	///	<summary>
	///		Get the auxiliary grid data with the specified index.
	///	</summary>
	const GridData3D & GetVerticalDynamicsAuxData(
		int ix,
		DataLocation loc = DataLocation_Node
	) const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}
		if ((ix < 0) || (ix > m_datavecAuxNode[1].size())) {
			_EXCEPTIONT("Invalid index in TracersData vector.");
		}
		if (loc == DataLocation_Node) {
			return m_datavecAuxNode[1][ix];

		} else if (loc == DataLocation_REdge) {
			return m_datavecAuxREdge[1][ix];

		} else {
			_EXCEPTIONT("Invalid DataLocation.");
		}
	}

	///	<summary>
	///		Get the pressure data.
	///	</summary>
	GridData3D & GetDataPressure() {
		return m_dataPressure;
	}

	///	<summary>
	///		Get the pressure data.
	///	</summary>
	const GridData3D & GetDataPressure() const {
		return m_dataPressure;
	}

	///	<summary>
	///		Get the alpha direction pressure derivative data.
	///	</summary>
	GridData3D & GetDataDaPressure() {
		return m_dataDaPressure;
	}

	///	<summary>
	///		Get the alpha direction pressure derivative data.
	///	</summary>
	const GridData3D & GetDataDaPressure() const {
		return m_dataDaPressure;
	}

	///	<summary>
	///		Get the beta direction pressure derivative data.
	///	</summary>
	GridData3D & GetDataDbPressure() {
		return m_dataDbPressure;
	}

	///	<summary>
	///		Get the beta direction pressure derivative data.
	///	</summary>
	const GridData3D & GetDataDbPressure() const {
		return m_dataDbPressure;
	}

	///	<summary>
	///		Get the pressure derivative data.
	///	</summary>
	GridData3D & GetDataDxPressure() {
		return m_dataDxPressure;
	}

	///	<summary>
	///		Get the pressure derivative data.
	///	</summary>
	const GridData3D & GetDataDxPressure() const {
		return m_dataDxPressure;
	}

	///	<summary>
	///		Get the vorticity data.
	///	</summary>
	GridData3D & GetDataVorticity() {
		return m_dataVorticity;
	}

	///	<summary>
	///		Get the vorticity data.
	///	</summary>
	const GridData3D & GetDataVorticity() const {
		return m_dataVorticity;
	}

	///	<summary>
	///		Get the divergence data.
	///	</summary>
	GridData3D & GetDataDivergence() {
		return m_dataDivergence;
	}

	///	<summary>
	///		Get the divergence data.
	///	</summary>
	const GridData3D & GetDataDivergence() const {
		return m_dataDivergence;
	}

	///	<summary>
	///		Get the temperature data.
	///	</summary>
	GridData3D & GetDataTemperature() {
		return m_dataTemperature;
	}

	///	<summary>
	///		Get the temperature data.
	///	</summary>
	const GridData3D & GetDataTemperature() const {
		return m_dataTemperature;
	}

	///	<summary>
	///		Get the Rayleigh friction strength.
	///	</summary>
	GridData3D & GetRayleighStrength(
		DataLocation loc = DataLocation_Node
	) {
		if (loc == DataLocation_Node) {
			return m_dataRayleighStrengthNode;
		} else if (loc == DataLocation_REdge) {
			return m_dataRayleighStrengthREdge;
		} else {
			_EXCEPTIONT("Invalid location");
		}
	}

	///	<summary>
	///		Get the Rayleigh friction strength.
	///	</summary>
	const GridData3D & GetRayleighStrength(
		DataLocation loc = DataLocation_Node
	) const {
		if (loc == DataLocation_Node) {
			return m_dataRayleighStrengthNode;
		} else if (loc == DataLocation_REdge) {
			return m_dataRayleighStrengthREdge;
		} else {
			_EXCEPTIONT("Invalid location");
		}
	}
/*
	///	<summary>
	///		Get the count of FluxBuffers.
	///	</summary>
	inline int GetFluxBufferCount() const {
		return static_cast<int>(m_vecFluxBuffers.size());
	}

	///	<summary>
	///		Get the specified FluxBuffer.
	///	</summary>
	FluxBuffer & GetFluxBuffer(int ix) {
		if ((ix < 0) || (ix >= GetFluxBufferCount())) {
			_EXCEPTIONT("FluxBuffer index out of range");
		}
		return m_vecFluxBuffers[ix];
	}

	///	<summary>
	///		Get the specified FluxBuffer.
	///	</summary>
	const FluxBuffer & GetFluxBuffer(int ix) const {
		if ((ix < 0) || (ix >= GetFluxBufferCount())) {
			_EXCEPTIONT("FluxBuffer index out of range");
		}
		return m_vecFluxBuffers[ix];
	}
*/
protected:
	///	<summary>
	///		Reference to parent grid.
	///	</summary>
	Grid & m_grid;

	///	<summary>
	///		Global index of this patch in the grid array.
	///	</summary>
	int m_ixPatch;

	///	<summary>
	///		Processor which stores this patch.
	///	</summary>
	int m_iProcessor;

	///	<summary>
	///		PatchBox defining the location of this patch on the grid.
	///	</summary>
	PatchBox m_box;

	///	<summary>
	///		Panels in each coordinate direction.
	///	</summary>
	PanelIndexVector m_ixNeighborPanel;

	///	<summary>
	///		Connectivity of this patch to other patches on the grid.
	///	</summary>
	Connectivity m_connect;

	///	<summary>
	///		Flag indicating that this patch contains data.
	///	</summary>
	bool m_fContainsData;

	///	<summary>
	///		2D Jacobian at each node.
	///	</summary>
	DataMatrix<double> m_dataJacobian2D;

	///	<summary>
	///		2D Contravariant metric (alpha) components.
	///	</summary>
	DataMatrix3D<double> m_dataContraMetric2DA;

	///	<summary>
	///		2D Contravariant metric (beta) components.
	///	</summary>
	DataMatrix3D<double> m_dataContraMetric2DB;

	///	<summary>
	///		Jacobian at each node.
	///	</summary>
	DataMatrix3D<double> m_dataJacobian;

	///	<summary>
	///		Jacobian at each edge.
	///	</summary>
	DataMatrix3D<double> m_dataJacobianREdge;

	///	<summary>
	///		Contravariant metric (alpha) components.
	///	</summary>
	DataMatrix4D<double> m_dataContraMetricA;

	///	<summary>
	///		Contravariant metric (beta) components.
	///	</summary>
	DataMatrix4D<double> m_dataContraMetricB;

	///	<summary>
	///		Contravariant metric (xi) components.
	///	</summary>
	DataMatrix4D<double> m_dataContraMetricXi;

	///	<summary>
	///		Christoffel symbol (alpha) components at each node.
	///	</summary>
	DataMatrix3D<double> m_dataChristoffelA;

	///	<summary>
	///		Christoffel symbol (beta) components at each node.
	///	</summary>
	DataMatrix3D<double> m_dataChristoffelB;

	///	<summary>
	///		Christoffel symbol (xi) components at each node.
	///	</summary>
	DataMatrix4D<double> m_dataChristoffelXi;

 	///	<summary>
	///		Orthonormalization coefficients at each node.
	///	</summary>
	DataMatrix4D<double> m_dataOrthonormNode;

	///	<summary>
	///		Orthonormalization coefficients at each interface.
	///	</summary>
	DataMatrix4D<double> m_dataOrthonormREdge;

	///	<summary>
	///		Element area at each node.
	///	</summary>
	DataMatrix3D<double> m_dataElementArea;

	///	<summary>
	///		Element area at each interface.
	///	</summary>
	DataMatrix3D<double> m_dataElementAreaREdge;

	///	<summary>
	///		Topography height at each node.
	///	</summary>
	DataMatrix<double> m_dataTopography;

	///	<summary>
	///		Longitude at each node.
	///	</summary>
	DataMatrix<double> m_dataLon;

	///	<summary>
	///		Latitude at each node.
	///	</summary>
	DataMatrix<double> m_dataLat;

	///	<summary>
	///		Coriolis parameter at each node.
	///	</summary>
	DataMatrix<double> m_dataCoriolisF;

	///	<summary>
	///		Altitude at each level.
	///	</summary>
	DataMatrix3D<double> m_dataZLevels;

	///	<summary>
	///		Altitude at each interface.
	///	</summary>
	DataMatrix3D<double> m_dataZInterfaces;

	///	<summary>
	///		Grid data for the reference state on model levels.
	///	</summary>
	GridData4D m_dataRefStateNode;

	///	<summary>
	///		Grid data for state variables on model levels.
	///	</summary>
	GridData4DVector m_datavecStateNode;

	///	<summary>
	///		Grid data on model interfaces for the reference state.
	///	</summary>
	GridData4D m_dataRefStateREdge;

	///	<summary>
	///		Grid data on model interfaces for state variables.
	///	</summary>
	GridData4DVector m_datavecStateREdge;

	///	<summary>
	///		Grid data for tracer variables.
	///	</summary>
	GridData4DVector m_datavecTracers;

	///	<summary>
	///		Auxiliary grid data on model levels.
	///	</summary>
	GridData3DVectorVector m_datavecAuxNode;

	///	<summary>
	///		Auxiliary grid data on model interfaces.
	///	</summary>
	GridData3DVectorVector m_datavecAuxREdge;

	///	<summary>
	///		Computed pointwise pressures.
	///	</summary>
	GridData3D m_dataPressure;

	///	<summary>
	///		Computed pointwise alpha pressure derivatives.
	///	</summary>
	GridData3D m_dataDaPressure;

	///	<summary>
	///		Computed pointwise beta pressure derivatives.
	///	</summary>
	GridData3D m_dataDbPressure;

	///	<summary>
	///		Computed pointwise vertical pressure derivatives.
	///	</summary>
	GridData3D m_dataDxPressure;

	///	<summary>
	///		Computed vorticity.
	///	</summary>
	GridData3D m_dataVorticity;

	///	<summary>
	///		Computed divergence.
	///	</summary>
	GridData3D m_dataDivergence;

	///	<summary>
	///		Computed temperature.
	///	</summary>
	GridData3D m_dataTemperature;

	///	<summary>
	///		Rayleigh friction strength on nodes.
	///	</summary>
	GridData3D m_dataRayleighStrengthNode;

	///	<summary>
	///		Rayleigh friction strength on interfaces.
	///	</summary>
	GridData3D m_dataRayleighStrengthREdge;

	///	<summary>
	///		Vector of flux buffers.
	///	</summary>
	FluxBufferVector m_vecFluxBuffers;

};

///////////////////////////////////////////////////////////////////////////////

//typedef std::vector<GridPatch*> GridPatchVector;

class GridPatchVector
	: public std::vector<GridPatch*>
{
	public:
		GridPatch * operator[](int i) const {
			if (i == GridPatch::InvalidIndex) {
				return NULL;
			}
			return std::vector<GridPatch*>::operator[](i);
		}
};

///////////////////////////////////////////////////////////////////////////////

#endif

