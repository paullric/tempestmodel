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

#include "DataContainer.h"
#include "DataArray1D.h"
#include "DataArray2D.h"
#include "DataArray3D.h"
#include "DataArray4D.h"

#include "PatchBox.h"
#include "ChecksumType.h"
#include "DataContainer.h"
#include "Connectivity.h"

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

protected:
	///	<summary>
	///		A vector storing DataArray4D<double>.
	///	</summary>
	typedef std::vector< DataArray4D<double> > DataArray4DVector;

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
/*
public:
	///	<summary>
	///		Total number of degrees of freedom in data.
	///	</summary>
	int GetTotalDegreesOfFreedom(
		DataType eDataType,
		DataLocation eDataLocation = DataLocation_None
	) const;
*/
public:
	///	<summary>
	///		Initialize patch as remote (data on an exterior processor).
	///	</summary>
	virtual void InitializeDataRemote(int iProcessor);

	///	<summary>
	///		Initialize patch as local (data on this processor).
	///	</summary>
	virtual void InitializeDataLocal(
		bool fAllocateGeometric = true,
		bool fAllocateActiveState = true,
		bool fAllocateBufferState = true,
		bool fAllocateAuxiliary = true
	);

	///	<summary>
	///		Delete data and reset patch to a stub.
	///	</summary>
	virtual void DeinitializeData();

	///	<summary>
	///		Initialize coordinate data on patch.
	///	</summary>
	virtual void InitializeCoordinateData() = 0;

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

	///	<summary>
	///		Apply boundary conditions to the state.
	///	</summary>
	virtual void ApplyBoundaryConditions(
		int iDataIndex = 0,
		DataType eDataType = DataType_State
	) {
	}

public:
	///	<summary>
	///		Compute the radial component of the curl on the grid given two
	///		contravariant vector fields.
	///	</summary>
	virtual void ComputeCurlAndDiv(
		const DataArray3D<double> & dataUa,
		const DataArray3D<double> & dataUb
	) {
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

	///	<summary>
	///		Compute temperature on the GridPatch.
	///	</summary>
	virtual void ComputeSurfacePressure(
		int iDataIndex
	);
/*
	///	<summary>
	///		Compute Richardson number on the GridPatch.
	///	</summary>
	virtual void ComputeRichardson(
		int iDataIndex,
		DataLocation loc = DataLocation_Node
	);
*/
public:
	///	<summary>
	///		Add local masses to checksum total.
	///	</summary>
	void Checksum(
		DataType eDataType,
		DataArray1D<double> & dChecksums,
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

	///	<summary>
	///		Compute total vertical momentum over this GridPatch.
	///	</summary>
	double ComputeTotalVerticalMomentum(
		int iDataIndex
	);

public:
	///	<summary>
	///		Pack the ExchangeBuffer.
	///	</summary>
	void PackExchangeBuffer(
		DataType eDataType,
		int iDataIndex,
		ExchangeBuffer & exbuf
	);

	///	<summary>
	///		Unpack the ExchangeBuffer.
	///	</summary>
	void UnpackExchangeBuffer(
		DataType eDataType,
		int iDataIndex,
		ExchangeBuffer & exbuf
	);

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
		const DataArray1D<double> & dCoeff,
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
	///		Linearly interpolate data horizontally to the specified points.
	///	</summary>
	virtual void InterpolateData(
		DataType eDataType,
		const DataArray1D<double> & dREta,
		const DataArray1D<double> & dAlpha,
		const DataArray1D<double> & dBeta,
		const DataArray1D<int> & iPatch,
		DataArray3D<double> & dInterpData,
		DataLocation eOnlyVariablesAt = DataLocation_None,
		bool fIncludeReferenceState = true,
		bool fConvertToPrimitive = true
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

public:
	///	<summary>
	///		Get the DataContainer for storing geometric data.
	///	</summary>
	DataContainer & GetDataContainerGeometric() {
		return m_dcGeometric;
	}

	///	<summary>
	///		Get the DataContainer for storing active state data.
	///	</summary>
	DataContainer & GetDataContainerActiveState() {
		return m_dcActiveState;
	}

	///	<summary>
	///		Get the DataContainer for storing buffer state data.
	///	</summary>
	DataContainer & GetDataContainerBufferState() {
		return m_dcBufferState;
	}

	///	<summary>
	///		Get the DataContainer for storing auxiliary data.
	///	</summary>
	DataContainer & GetDataContainerAuxiliary() {
		return m_dcAuxiliary;
	}

public:
	///	<summary>
	///		Get the DataContainer for storing geometric data.
	///	</summary>
	const DataContainer & GetDataContainerGeometric() const {
		return m_dcGeometric;
	}

	///	<summary>
	///		Get the DataContainer for storing active state data.
	///	</summary>
	const DataContainer & GetDataContainerActiveState() const {
		return m_dcActiveState;
	}

	///	<summary>
	///		Get the DataContainer for storing buffer state data.
	///	</summary>
	const DataContainer & GetDataContainerBufferState() const {
		return m_dcBufferState;
	}

	///	<summary>
	///		Get the DataContainer for storing auxiliary data.
	///	</summary>
	const DataContainer & GetDataContainerAuxiliary() const {
		return m_dcAuxiliary;
	}

public:
	///	<summary>
	///		Returns true if this GridPatch contains data, ie. is not a stub.
	///	</summary>
	bool ContainsData() const {
		return m_fContainsData;
	}

public:
	///	<summary>
	///		Get the alpha node with specified local index.
	///	</summary>
	inline double GetANode(int ix) const {
		return m_dANode[ix];
	}

	///	<summary>
	///		Get the alpha edge with specified local alpha edge index.
	///	</summary>
	inline double GetAEdge(int ix) const {
		return m_dAEdge[ix];
	}

	///	<summary>
	///		Get the array of alpha nodes.
	///	</summary>
	inline const DataArray1D<double> & GetANodes() const {
		return m_dANode;
	}

	///	<summary>
	///		Get the array of alpha edges.
	///	</summary>
	inline const DataArray1D<double> & GetAEdges() const {
		return m_dAEdge;
	}

	///	<summary>
	///		Get the beta node with specified local index.
	///	</summary>
	inline double GetBNode(int ix) const {
		return m_dBNode[ix];
	}

	///	<summary>
	///		Get the beta edge with specified local beta edge index.
	///	</summary>
	inline double GetBEdge(int ix) const {
		return m_dBEdge[ix];
	}

	///	<summary>
	///		Get the array of beta nodes.
	///	</summary>
	inline const DataArray1D<double> & GetBNodes() const {
		return m_dBNode;
	}

	///	<summary>
	///		Get the array of beta edges.
	///	</summary>
	inline const DataArray1D<double> & GetBEdges() const {
		return m_dBEdge;
	}

	///	<summary>
	///		Get the 2D Jacobian matrix.
	///	</summary>
	const DataArray2D<double> & GetJacobian2D() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataJacobian2D;
	}

	///	<summary>
	///		Get the components of the contravariant metric (alpha)
	///	</summary>
	const DataArray3D<double> & GetContraMetric2DA() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataContraMetric2DA;
	}

	///	<summary>
	///		Get the components of the contravariant metric (beta)
	///	</summary>
	const DataArray3D<double> & GetContraMetric2DB() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataContraMetric2DB;
	}

	///	<summary>
	///		Get the components of the covariant metric (alpha)
	///	</summary>
	const DataArray3D<double> & GetCovMetric2DA() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataCovMetric2DA;
	}

	///	<summary>
	///		Get the components of the covariant metric (beta)
	///	</summary>
	const DataArray3D<double> & GetCovMetric2DB() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataCovMetric2DB;
	}

	///	<summary>
	///		Get the nodal Jacobian matrix.
	///	</summary>
	const DataArray3D<double> & GetJacobian() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataJacobian;
	}

	///	<summary>
	///		Get the interface Jacobian matrix.
	///	</summary>
	const DataArray3D<double> & GetJacobianREdge() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataJacobianREdge;
	}

	///	<summary>
	///		Get the components of the contravariant metric (alpha)
	///	</summary>
	const DataArray4D<double> & GetContraMetricA() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataContraMetricA;
	}

	///	<summary>
	///		Get the components of the contravariant metric (beta)
	///	</summary>
	const DataArray4D<double> & GetContraMetricB() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataContraMetricB;
	}

	///	<summary>
	///		Get the components of the contravariant metric (xi)
	///	</summary>
	const DataArray4D<double> & GetContraMetricXi() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataContraMetricXi;
	}

	///	<summary>
	///		Get the components of the contravariant metric (alpha)
	///		on interfaces.
	///	</summary>
	const DataArray4D<double> & GetContraMetricAREdge() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataContraMetricAREdge;
	}

	///	<summary>
	///		Get the components of the contravariant metric (beta)
	///		on interfaces.
	///	</summary>
	const DataArray4D<double> & GetContraMetricBREdge() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataContraMetricBREdge;
	}

	///	<summary>
	///		Get the components of the contravariant metric (xi)
	///		on interfaces.
	///	</summary>
	const DataArray4D<double> & GetContraMetricXiREdge() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataContraMetricXiREdge;
	}
/*
	///	<summary>
	///		Get the components of the covariant metric (alpha)
	///	</summary>
	const DataArray4D<double> & GetCovMetricA() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataCovMetricA;
	}

	///	<summary>
	///		Get the  components of the covariant metric (beta)
	///	</summary>
	const DataArray4D<double> & GetCovMetricB() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataCovMetricB;
	}

	///	<summary>
	///		Get the  components of the covariant metric (xi)
	///	</summary>
	const DataArray4D<double> & GetCovMetricXi() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataCovMetricXi;
	}
*/
 	///	<summary>
	///		Get the vertical coordinate transform (derivatives of the
	///		radius) at nodes.
	///	</summary>
	const DataArray4D<double> & GetDerivRNode() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataDerivRNode;
	}

 	///	<summary>
	///		Get the vertical coordinate transform (derivatives of the
	///		radius) at edges.
	///	</summary>
	const DataArray4D<double> & GetDerivRREdge() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataDerivRREdge;
	}

	///	<summary>
	///		Get the nodal element area matrix.
	///	</summary>
	const DataArray3D<double> & GetElementArea() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataElementArea;
	}

	///	<summary>
	///		Get the interface element area matrix.
	///	</summary>
	const DataArray3D<double> & GetElementAreaREdge() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataElementAreaREdge;
	}

	///	<summary>
	///		Get the nodal topography.
	///	</summary>
	DataArray2D<double> & GetTopography() {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataTopography;
	}

	///	<summary>
	///		Get the nodal topography.
	///	</summary>
	const DataArray2D<double> & GetTopography() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataTopography;
	}

	///	<summary>
	///		Get the derivatives of the topography.
	///	</summary>
	DataArray3D<double> & GetTopographyDeriv() {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataTopographyDeriv;
	}

	///	<summary>
	///		Get the derivatives of the topography.
	///	</summary>
	const DataArray3D<double> & GetTopographyDeriv() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataTopographyDeriv;
	}

	///	<summary>
	///		Get the nodal longitude matrix.
	///	</summary>
	const DataArray2D<double> & GetLongitude() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataLon;
	}

	///	<summary>
	///		Get the nodal latitude matrix.
	///	</summary>
	const DataArray2D<double> & GetLatitude() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataLat;
	}

	///	<summary>
	///		Get the nodal Coriolis parameter.
	///	</summary>
	const DataArray2D<double> & GetCoriolisF() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataCoriolisF;
	}

	///	<summary>
	///		Get the radial coordinate matrix on model levels.
	///	</summary>
	DataArray3D<double> & GetZLevels() {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataZLevels;
	}

	///	<summary>
	///		Get the radial coordinate matrix on model levels.
	///	</summary>
	const DataArray3D<double> & GetZLevels() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataZLevels;
	}

	///	<summary>
	///		Get the radial coordinate matrix on model interfaces.
	///	</summary>
	DataArray3D<double> & GetZInterfaces() {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataZInterfaces;
	}

	///	<summary>
	///		Get the radial coordinate matrix on model interfaces.
	///	</summary>
	const DataArray3D<double> & GetZInterfaces() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}

		return m_dataZInterfaces;
	}

	///	<summary>
	///		Get the reference state.
	///	</summary>
	DataArray4D<double> & GetReferenceState(
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
	const DataArray4D<double> & GetReferenceState(
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
	DataArray4D<double> & GetDataState(
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
	const DataArray4D<double> & GetDataState(
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
	DataArray4D<double> & GetDataTracers(int ix) {
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
	const DataArray4D<double> & GetDataTracers(int ix) const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}
		if ((ix < 0) || (ix > m_datavecTracers.size())) {
			_EXCEPTIONT("Invalid index in TracersData vector.");
		}
		return m_datavecTracers[ix];
	}

	///	<summary>
	///		Get the tracer data reference state.
	///	</summary>
	DataArray4D<double> & GetReferenceTracers() {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}
		return m_dataRefTracers;
	}

	///	<summary>
	///		Get the tracer data reference state.
	///	</summary>
	const DataArray4D<double> & GetReferenceTracers() const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}
		return m_dataRefTracers;
	}

	///	<summary>
	///		Get the tracer data matrix with the specified index.
	///	</summary>
	DataArray4D<double> & GetDataTracersReference(int ix) {
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
	const DataArray4D<double> & GetDataTracersReference(int ix) const {
		if (!m_fContainsData) {
			_EXCEPTIONT("Stub patch does not store data.");
		}
		if ((ix < 0) || (ix > m_datavecTracers.size())) {
			_EXCEPTIONT("Invalid index in TracersData vector.");
		}
		return m_datavecTracers[ix];
	}

	///	<summary>
	///		Get the 2D user data.
	///	</summary>
	DataArray3D<double> & GetUserData2D() {
		return m_dataUserData2D;
	}

	///	<summary>
	///		Get the 2D user data.
	///	</summary>
	const DataArray3D<double> & GetUserData2D() const {
		return m_dataUserData2D;
	}

	///	<summary>
	///		Get the pressure data.
	///	</summary>
	DataArray3D<double> & GetDataPressure() {
		return m_dataPressure;
	}

	///	<summary>
	///		Get the pressure data.
	///	</summary>
	const DataArray3D<double> & GetDataPressure() const {
		return m_dataPressure;
	}

	///	<summary>
	///		Get the pressure derivative data.
	///	</summary>
	DataArray3D<double> & GetDataDxPressure() {
		return m_dataDxPressure;
	}

	///	<summary>
	///		Get the pressure derivative data.
	///	</summary>
	const DataArray3D<double> & GetDataDxPressure() const {
		return m_dataDxPressure;
	}

	///	<summary>
	///		Get the vorticity data.
	///	</summary>
	DataArray3D<double> & GetDataVorticity() {
		return m_dataVorticity;
	}

	///	<summary>
	///		Get the vorticity data.
	///	</summary>
	const DataArray3D<double> & GetDataVorticity() const {
		return m_dataVorticity;
	}

	///	<summary>
	///		Get the divergence data.
	///	</summary>
	DataArray3D<double> & GetDataDivergence() {
		return m_dataDivergence;
	}

	///	<summary>
	///		Get the divergence data.
	///	</summary>
	const DataArray3D<double> & GetDataDivergence() const {
		return m_dataDivergence;
	}

	///	<summary>
	///		Get the temperature data.
	///	</summary>
	DataArray3D<double> & GetDataTemperature() {
		return m_dataTemperature;
	}

	///	<summary>
	///		Get the temperature data.
	///	</summary>
	const DataArray3D<double> & GetDataTemperature() const {
		return m_dataTemperature;
	}

	///	<summary>
	///		Get the surface pressure data.
	///	</summary>
	DataArray2D<double> & GetDataSurfacePressure() {
		return m_dataSurfacePressure;
	}

	///	<summary>
	///		Get the surface pressure data.
	///	</summary>
	const DataArray2D<double> & GetDataSurfacePressure() const {
		return m_dataSurfacePressure;
	}

	///	<summary>
	///		Get the Rayleigh friction strength.
	///	</summary>
	DataArray3D<double> & GetRayleighStrength(
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
	const DataArray3D<double> & GetRayleighStrength(
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

protected:
	///	<summary>
	///		DataContainer for storing geometric data.
	///	</summary>
	DataContainer m_dcGeometric;

	///	<summary>
	///		DataContainer for storing active state data.
	///	</summary>
	DataContainer m_dcActiveState;

	///	<summary>
	///		DataContainer for storing buffer state data.
	///	</summary>
	DataContainer m_dcBufferState;

	///	<summary>
	///		DataContainer for storing auxiliary data.
	///	</summary>
	DataContainer m_dcAuxiliary;

	///	<summary>
	///		Flag indicating that this patch contains data.
	///	</summary>
	bool m_fContainsData;

protected:
	///	<summary>
	///		Geometric patch index.
	///	</summary>
	DataArray1D<size_t> m_iGeometricPatchIx;

	///	<summary>
	///		Array of nodal values in the alpha direction (Geometric).
	///	</summary>
	DataArray1D<double> m_dANode;

	///	<summary>
	///		Array of edge values in the alpha direction (Geometric).
	///	</summary>
	DataArray1D<double> m_dAEdge;

	///	<summary>
	///		Array of nodal values in the beta direction (Geometric).
	///	</summary>
	DataArray1D<double> m_dBNode;

	///	<summary>
	///		Array of edge values in the beta direction (Geometric).
	///	</summary>
	DataArray1D<double> m_dBEdge;

	///	<summary>
	///		2D Jacobian at each node (Geometric).
	///	</summary>
	DataArray2D<double> m_dataJacobian2D;

	///	<summary>
	///		2D Contravariant metric (alpha) components (Geometric).
	///	</summary>
	DataArray3D<double> m_dataContraMetric2DA;

	///	<summary>
	///		2D Contravariant metric (beta) components (Geometric).
	///	</summary>
	DataArray3D<double> m_dataContraMetric2DB;

	///	<summary>
	///		2D Covariant metric (alpha) components (Geometric).
	///	</summary>
	DataArray3D<double> m_dataCovMetric2DA;

	///	<summary>
	///		2D Covariant metric (beta) components (Geometric).
	///	</summary>
	DataArray3D<double> m_dataCovMetric2DB;

	///	<summary>
	///		Jacobian at each node (Geometric).
	///	</summary>
	DataArray3D<double> m_dataJacobian;

	///	<summary>
	///		Jacobian at each edge (Geometric).
	///	</summary>
	DataArray3D<double> m_dataJacobianREdge;

	///	<summary>
	///		Contravariant metric (alpha) components (Geometric).
	///	</summary>
	DataArray4D<double> m_dataContraMetricA;

	///	<summary>
	///		Contravariant metric (beta) components (Geometric).
	///	</summary>
	DataArray4D<double> m_dataContraMetricB;

	///	<summary>
	///		Contravariant metric (xi) components (Geometric).
	///	</summary>
	DataArray4D<double> m_dataContraMetricXi;
/*
	///	<summary>
	///		Covariant metric (alpha) components (Geometric).
	///	</summary>
	DataArray4D<double> m_dataCovMetricA;

	///	<summary>
	///		Covariant metric (beta) components (Geometric).
	///	</summary>
	DataArray4D<double> m_dataCovMetricB;

	///	<summary>
	///		Covariant metric (xi) components (Geometric).
	///	</summary>
	DataArray4D<double> m_dataCovMetricXi;
*/
	///	<summary>
	///		Contravariant metric (alpha) components on interfaces (Geometric).
	///	</summary>
	DataArray4D<double> m_dataContraMetricAREdge;

	///	<summary>
	///		Contravariant metric (beta) components on interfaces (Geometric).
	///	</summary>
	DataArray4D<double> m_dataContraMetricBREdge;

	///	<summary>
	///		Contravariant metric (xi) components on interfaces (Geometric).
	///	</summary>
	DataArray4D<double> m_dataContraMetricXiREdge;

 	///	<summary>
	///		Vertical coordinate transform (derivatives of the radius)
	///		at each node (Geometric).
	///	</summary>
	DataArray4D<double> m_dataDerivRNode;

 	///	<summary>
	///		Vertical coordinate transform (derivatives of the radius)
	///		at each interface (Geometric).
	///	</summary>
	DataArray4D<double> m_dataDerivRREdge;

	///	<summary>
	///		Element area at each node (Geometric).
	///	</summary>
	DataArray3D<double> m_dataElementArea;

	///	<summary>
	///		Element area at each interface (Geometric).
	///	</summary>
	DataArray3D<double> m_dataElementAreaREdge;

	///	<summary>
	///		Topography height at each node (Geometric).
	///	</summary>
	DataArray2D<double> m_dataTopography;

	///	<summary>
	///		Derivatives of topography at each node (Geometric).
	///	</summary>
	DataArray3D<double> m_dataTopographyDeriv;

	///	<summary>
	///		Longitude at each node (Geometric).
	///	</summary>
	DataArray2D<double> m_dataLon;

	///	<summary>
	///		Latitude at each node (Geometric).
	///	</summary>
	DataArray2D<double> m_dataLat;

	///	<summary>
	///		Coriolis parameter at each node (Geometric).
	///	</summary>
	DataArray2D<double> m_dataCoriolisF;

	///	<summary>
	///		Altitude at each level (Geometric).
	///	</summary>
	DataArray3D<double> m_dataZLevels;

	///	<summary>
	///		Altitude at each interface (Geometric).
	///	</summary>
	DataArray3D<double> m_dataZInterfaces;

	///	<summary>
	///		Rayleigh friction strength on nodes (Geometric).
	///	</summary>
	DataArray3D<double> m_dataRayleighStrengthNode;

	///	<summary>
	///		Rayleigh friction strength on interfaces (Geometric).
	///	</summary>
	DataArray3D<double> m_dataRayleighStrengthREdge;

public:
	///	<summary>
	///		Active State patch index.
	///	</summary>
	DataArray1D<size_t> m_iActiveStatePatchIx;

	///	<summary>
	///		Grid data for the reference state on model levels (State).
	///	</summary>
	DataArray4D<double> m_dataRefStateNode;

	///	<summary>
	///		Grid data on model interfaces for the reference state (State).
	///	</summary>
	DataArray4D<double> m_dataRefStateREdge;

	///	<summary>
	///		Grid data for state variables on model levels (State).
	///	</summary>
	DataArray4DVector m_datavecStateNode;

	///	<summary>
	///		Grid data on model interfaces for state variables (State).
	///	</summary>
	DataArray4DVector m_datavecStateREdge;

	///	<summary>
	///		Grid data for tracer variables (State).
	///	</summary>
	DataArray4DVector m_datavecTracers;

	///	<summary>
	///		Grid data for the reference state on model levels (State).
	///	</summary>
	DataArray4D<double> m_dataRefTracers;

public:
	///	<summary>
	///		Computed pointwise pressures (Auxiliary).
	///	</summary>
	DataArray3D<double> m_dataPressure;

	///	<summary>
	///		Computed pointwise vertical pressure derivatives (Auxiliary).
	///	</summary>
	DataArray3D<double> m_dataDxPressure;

	///	<summary>
	///		Computed vorticity (Auxiliary).
	///	</summary>
	DataArray3D<double> m_dataVorticity;

	///	<summary>
	///		Computed divergence (Auxiliary).
	///	</summary>
	DataArray3D<double> m_dataDivergence;

	///	<summary>
	///		Computed temperature (Auxiliary).
	///	</summary>
	DataArray3D<double> m_dataTemperature;

	///	<summary>
	///		Computed Richardson number (Auxiliary).
	///	</summary>
	DataArray3D<double> m_dataRichardson;

	///	<summary>
	///		Computed surface pressure (Auxiliary).
	///	</summary>
	DataArray2D<double> m_dataSurfacePressure;

public:
	///	<summary>
	///		2D user data (Auxiliary).
	///	</summary>
	DataArray3D<double> m_dataUserData2D;
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

