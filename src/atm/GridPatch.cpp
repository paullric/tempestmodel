///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridPatch.cpp
///	\author  Paul Ullrich
///	\version February 26, 2013
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

#include "GridPatch.h"
#include "Grid.h"
#include "Model.h"
#include "EquationSet.h"
#include "Defines.h"

#include <mpi.h>

///////////////////////////////////////////////////////////////////////////////

GridPatch::GridPatch(
	Grid & grid,
	int ixPatch,
	const PatchBox & box
) :
	m_grid(grid),
	m_ixPatch(ixPatch),
	m_iProcessor(0),
	m_box(box),
	m_connect(*this),
	m_fContainsData(false)
{
}

///////////////////////////////////////////////////////////////////////////////

void GridPatch::InitializeDataRemote(
	int iProcessor
) {
	// Remove existing data
	if (m_fContainsData) {
		DeinitializeData();
	}

	// Set the processor
	m_iProcessor = iProcessor;
}

///////////////////////////////////////////////////////////////////////////////

void GridPatch::InitializeDataLocal(
	bool fAllocateGeometric,
	bool fAllocateActiveState,
	bool fAllocateBufferState,
	bool fAllocateAuxiliary
) {
	if (m_fContainsData) {
		_EXCEPTIONT(
			"Attempting to initialize a previously initialized GridPatch.");
	}

	// This patch contains data
	m_fContainsData = true;

	// Set the processor
	MPI_Comm_rank(MPI_COMM_WORLD, &m_iProcessor);

	// Get the model
	const Model & model = m_grid.GetModel();

	// Get the equation set
	const EquationSet & eqn = model.GetEquationSet();

	// 2D user data metadata
	const UserDataMeta & metaUserData = model.GetUserDataMeta();

	// Geometric patch index
	m_iGeometricPatchIx.SetSize(1);

	// Nodal coordinates in alpha direction
	m_dANode.SetSize(m_box.GetATotalWidth());

	// Edge coordinates in alpha direction
	m_dAEdge.SetSize(m_box.GetATotalWidth()+1);

	// Nodal coordinates in beta direction
	m_dBNode.SetSize(m_box.GetBTotalWidth());

	// Edge coordinates in beta direction
	m_dBEdge.SetSize(m_box.GetBTotalWidth()+1);

	// Jacobian at each node (2D)
	m_dataJacobian2D.SetSize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth());

	// Contravariant metric (2D) components at each node
	m_dataContraMetric2DA.SetSize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		2);

	m_dataContraMetric2DB.SetSize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		2);

	// Covariant metric (2D) components at each node
	m_dataCovMetric2DA.SetSize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		2);

	m_dataCovMetric2DB.SetSize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		2);

	// Jacobian at each node
	m_dataJacobian.SetSize(
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth());

	// Jacobian at each interface
	m_dataJacobianREdge.SetSize(
		m_grid.GetRElements()+1,
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth());

	// Contravariant metric components at each node
	m_dataContraMetricA.SetSize(
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		3);

	m_dataContraMetricB.SetSize(
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		3);

	m_dataContraMetricXi.SetSize(
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		3);

	// Covariant metric components at each node
	m_dataCovMetricA.SetSize(
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		3);

	m_dataCovMetricB.SetSize(
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		3);

	m_dataCovMetricXi.SetSize(
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		3);

	// Xi contravariant metric on interfaces
	m_dataContraMetricAREdge.SetSize(
		m_grid.GetRElements()+1,
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		3);

	m_dataContraMetricBREdge.SetSize(
		m_grid.GetRElements()+1,
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		3);

	m_dataContraMetricXiREdge.SetSize(
		m_grid.GetRElements()+1,
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		3);

	// Vertical coordinate transform (derivatives of the radius)
	m_dataDerivRNode.SetSize(
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		3);

	m_dataDerivRREdge.SetSize(
		m_grid.GetRElements()+1,
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		3);

	// Element area at each node
	m_dataElementArea.SetSize(
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth());

	// Element area at each interface
	m_dataElementAreaREdge.SetSize(
		m_grid.GetRElements()+1,
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth());

	// Topography height at each node
	m_dataTopography.SetSize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth());

	// Topography derivatives at each node
	m_dataTopographyDeriv.SetDataType(DataType_TopographyDeriv);
	m_dataTopographyDeriv.SetSize(
		2,
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth());

	// Longitude at each node
	m_dataLon.SetSize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth());

	// Latitude at each node
	m_dataLat.SetSize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth());

	// Coriolis parameter at each node
	m_dataCoriolisF.SetSize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth());

	// Radial coordinate at each level
	m_dataZLevels.SetSize(
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth());

	// Radial coordinate at each interface
	m_dataZInterfaces.SetSize(
		m_grid.GetRElements()+1,
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth());

	// Active State Patch Index
	m_iActiveStatePatchIx.SetSize(1);

	// Initialize reference state
	m_dataRefStateNode.SetDataType(DataType_State);
	m_dataRefStateNode.SetDataLocation(DataLocation_Node);
	m_dataRefStateNode.SetSize(
		eqn.GetComponents(),
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth());

	m_dataRefStateREdge.SetDataType(DataType_State);
	m_dataRefStateREdge.SetDataLocation(DataLocation_REdge);
	m_dataRefStateREdge.SetSize(
		eqn.GetComponents(),
		m_grid.GetRElements()+1,
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth());

	m_dataRefTracers.SetDataType(DataType_Tracers);
	m_dataRefTracers.SetDataLocation(DataLocation_Node);
	m_dataRefTracers.SetSize(
		eqn.GetTracers(),
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth());

	// Rayleigh friction strength
	m_dataRayleighStrengthNode.SetDataType(DataType_None);
	m_dataRayleighStrengthNode.SetDataLocation(DataLocation_Node);
	m_dataRayleighStrengthNode.SetSize(
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth());

	// Rayleigh friction strength
	m_dataRayleighStrengthREdge.SetDataType(DataType_None);
	m_dataRayleighStrengthREdge.SetDataLocation(DataLocation_Node);
	m_dataRayleighStrengthREdge.SetSize(
		m_grid.GetRElements()+1,
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth());

	// Put all read-only data objects into DataContainer
	m_dcGeometric.PushDataChunk(&m_iGeometricPatchIx);
	m_dcGeometric.PushDataChunk(&m_dANode);
	m_dcGeometric.PushDataChunk(&m_dAEdge);
	m_dcGeometric.PushDataChunk(&m_dBNode);
	m_dcGeometric.PushDataChunk(&m_dBEdge);
	m_dcGeometric.PushDataChunk(&m_dataJacobian2D);
	m_dcGeometric.PushDataChunk(&m_dataContraMetric2DA);
	m_dcGeometric.PushDataChunk(&m_dataContraMetric2DB);
	m_dcGeometric.PushDataChunk(&m_dataCovMetric2DA);
	m_dcGeometric.PushDataChunk(&m_dataCovMetric2DB);
	m_dcGeometric.PushDataChunk(&m_dataJacobian);
	m_dcGeometric.PushDataChunk(&m_dataJacobianREdge);
	m_dcGeometric.PushDataChunk(&m_dataContraMetricA);
	m_dcGeometric.PushDataChunk(&m_dataContraMetricB);
	m_dcGeometric.PushDataChunk(&m_dataContraMetricXi);
	m_dcGeometric.PushDataChunk(&m_dataCovMetricA);
	m_dcGeometric.PushDataChunk(&m_dataCovMetricB);
	m_dcGeometric.PushDataChunk(&m_dataCovMetricXi);
	m_dcGeometric.PushDataChunk(&m_dataContraMetricAREdge);
	m_dcGeometric.PushDataChunk(&m_dataContraMetricBREdge);
	m_dcGeometric.PushDataChunk(&m_dataContraMetricXiREdge);
	m_dcGeometric.PushDataChunk(&m_dataDerivRNode);
	m_dcGeometric.PushDataChunk(&m_dataDerivRREdge);
	m_dcGeometric.PushDataChunk(&m_dataElementArea);
	m_dcGeometric.PushDataChunk(&m_dataElementAreaREdge);
	m_dcGeometric.PushDataChunk(&m_dataTopography);
	m_dcGeometric.PushDataChunk(&m_dataTopographyDeriv);
	m_dcGeometric.PushDataChunk(&m_dataLon);
	m_dcGeometric.PushDataChunk(&m_dataLat);
	m_dcGeometric.PushDataChunk(&m_dataCoriolisF);
	m_dcGeometric.PushDataChunk(&m_dataZLevels);
	m_dcGeometric.PushDataChunk(&m_dataZInterfaces);
	m_dcGeometric.PushDataChunk(&m_dataRefStateNode);
	m_dcGeometric.PushDataChunk(&m_dataRefStateREdge);
	m_dcGeometric.PushDataChunk(&m_dataRefTracers);

	m_dcGeometric.PushDataChunk(&m_dataRayleighStrengthNode);
	m_dcGeometric.PushDataChunk(&m_dataRayleighStrengthREdge);

	if (fAllocateGeometric) {
		m_dcGeometric.Allocate();
	}

	// Initialize component data
	m_datavecStateNode .resize(model.GetComponentDataInstances());
	m_datavecStateREdge.resize(model.GetComponentDataInstances());

	if (model.GetComponentDataInstances() < 1) {
		_EXCEPTIONT("At least one ComponentDataInstance required");
	}

	for (int m = 0; m < model.GetComponentDataInstances(); m++) {
		m_datavecStateNode[m].SetDataType(DataType_State);
		m_datavecStateNode[m].SetDataLocation(DataLocation_Node);
		m_datavecStateNode[m].SetSize(
			eqn.GetComponents(),
			m_grid.GetRElements(),
			m_box.GetATotalWidth(),
			m_box.GetBTotalWidth());

		m_datavecStateREdge[m].SetDataType(DataType_State);
		m_datavecStateREdge[m].SetDataLocation(DataLocation_REdge);
		m_datavecStateREdge[m].SetSize(
			eqn.GetComponents(),
			m_grid.GetRElements()+1,
			m_box.GetATotalWidth(),
			m_box.GetBTotalWidth());
	}

	m_dcActiveState.PushDataChunk(&m_iActiveStatePatchIx);
	m_dcActiveState.PushDataChunk(&m_datavecStateNode[0]);
	m_dcActiveState.PushDataChunk(&m_datavecStateREdge[0]);

	for (int m = 1; m < model.GetComponentDataInstances(); m++) {
		m_dcBufferState.PushDataChunk(&m_datavecStateNode[m]);
		m_dcBufferState.PushDataChunk(&m_datavecStateREdge[m]);
	}

	// Initialize tracer data
	m_datavecTracers.resize(model.GetTracerDataInstances());

	if (model.GetTracerDataInstances() < 1) {
		_EXCEPTIONT("At least one TracerDataInstance required");
	}

	if (eqn.GetTracers() != 0) {
		for (int m = 0; m < model.GetTracerDataInstances(); m++) {
			m_datavecTracers[m].SetDataType(DataType_Tracers);
			m_datavecTracers[m].SetDataLocation(DataLocation_Node);
			m_datavecTracers[m].SetSize(
				eqn.GetTracers(),
				m_grid.GetRElements(),
				m_box.GetATotalWidth(),
				m_box.GetBTotalWidth());
		}

		// Store active state data
		m_dcActiveState.PushDataChunk(&m_datavecTracers[0]);

		for (int m = 1; m < model.GetTracerDataInstances(); m++) {
			m_dcBufferState.PushDataChunk(&m_datavecTracers[m]);
		}
	}

	// Allocate active and buffer state data
	if (fAllocateActiveState) {
		m_dcActiveState.Allocate();
	}
	if (fAllocateBufferState) {
		m_dcBufferState.Allocate();
	}

	// Pressure data
	m_dataPressure.SetDataType(DataType_Pressure);
	m_dataPressure.SetDataLocation(DataLocation_Node);
	m_dataPressure.SetSize(
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth());

	m_dataDxPressure.SetDataType(DataType_Pressure);
	m_dataDxPressure.SetDataLocation(DataLocation_Node);
	m_dataDxPressure.SetSize(
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth());

	// Vorticity data
	m_dataVorticity.SetDataType(DataType_Vorticity);
	m_dataVorticity.SetDataLocation(DataLocation_Node);
	m_dataVorticity.SetSize(
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth());

	// Divergence data
	m_dataDivergence.SetDataType(DataType_Divergence);
	m_dataDivergence.SetDataLocation(DataLocation_Node);
	m_dataDivergence.SetSize(
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth());

	// Temperature data
	m_dataTemperature.SetDataType(DataType_Temperature);
	m_dataTemperature.SetDataLocation(DataLocation_Node);
	m_dataTemperature.SetSize(
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth());

	// Surface pressure data
	m_dataSurfacePressure.SetDataType(DataType_SurfacePressure);
	m_dataSurfacePressure.SetDataLocation(DataLocation_None);
	m_dataSurfacePressure.SetSize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth());

	// Put auxiliary date into DataContainer
	m_dcAuxiliary.PushDataChunk(&m_dataPressure);
	m_dcAuxiliary.PushDataChunk(&m_dataDxPressure);
	m_dcAuxiliary.PushDataChunk(&m_dataVorticity);
	m_dcAuxiliary.PushDataChunk(&m_dataDivergence);
	m_dcAuxiliary.PushDataChunk(&m_dataTemperature);
	m_dcAuxiliary.PushDataChunk(&m_dataSurfacePressure);

	// 2D User data
	if (metaUserData.GetUserData2DItemCount() != 0) {
		m_dataUserData2D.SetDataType(DataType_Auxiliary2D);
		m_dataUserData2D.SetDataLocation(DataLocation_None);
		m_dataUserData2D.SetSize(
			metaUserData.GetUserData2DItemCount(),
			m_box.GetATotalWidth(),
			m_box.GetBTotalWidth());

		m_dcAuxiliary.PushDataChunk(&m_dataUserData2D);
	}

	// Allocate auxiliary data
	if (fAllocateAuxiliary) {
		m_dcAuxiliary.Allocate();
	}

#pragma message "Remove these two lines?"
	// Mark Patch index
	m_iGeometricPatchIx[0] = m_ixPatch;
	m_iActiveStatePatchIx[0] = m_ixPatch;
}

///////////////////////////////////////////////////////////////////////////////

void GridPatch::DeinitializeData() {
	if (!m_fContainsData) {
		_EXCEPTIONT("Attempting to deinitialize a stub GridPatch.");
	}

	m_fContainsData = false;

	_EXCEPTIONT("Not implemented");
/*
	m_dataJacobian2D.Deallocate();
	m_dataContraMetric2DA.Deallocate();
	m_dataContraMetric2DB.Deallocate();
	m_dataCovMetric2DA.Deallocate();
	m_dataCovMetric2DB.Deallocate();

	m_dataJacobian.Deallocate();
	m_dataJacobianREdge.Deallocate();
	m_dataContraMetricA.Deallocate();
	m_dataContraMetricB.Deallocate();
	m_dataContraMetricXi.Deallocate();
	m_dataCovMetricA.Deallocate();
	m_dataCovMetricB.Deallocate();
	m_dataCovMetricXi.Deallocate();
	m_dataContraMetricXiREdge.Deallocate();
	m_dataDerivRNode.Deallocate();
	m_dataElementArea.Deallocate();
	m_dataElementAreaREdge.Deallocate();
	m_dataTopography.Deallocate();
	m_dataTopographyDeriv.Deallocate();

	m_dataLon.Deallocate();
	m_dataLat.Deallocate();
	m_dataZLevels.Deallocate();
	m_dataZInterfaces.Deallocate();

	m_datavecStateNode.clear();
	m_datavecStateREdge.clear();
	m_datavecTracers.clear();

	m_dataPressure.Deallocate();
	m_dataDxPressure.Deallocate();

	m_dataVorticity.Deallocate();
	m_dataDivergence.Deallocate();
	m_dataTemperature.Deallocate();
	m_dataRayleighStrengthNode.Deallocate();
	m_dataRayleighStrengthREdge.Deallocate();
*/
}

///////////////////////////////////////////////////////////////////////////////

void GridPatch::ComputeSurfacePressure(
	int iDataIndex
) {
	const PhysicalConstants & phys = m_grid.GetModel().GetPhysicalConstants();

	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int PIx = 2;
	const int WIx = 3;
	const int RIx = 4;

	if (m_grid.GetModel().GetEquationSet().GetComponents() < 5) {
		_EXCEPTIONT("Invalid EquationSet.");
	}

	int k;
	int i;
	int j;

	// Calculate hydrostatic surface pressure on nodes
	const DataArray4D<double> & dataNode = m_datavecStateNode[iDataIndex];

	for (i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
	for (j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {

		m_dataSurfacePressure[i][j] = 0.0;

		for (k = 0; k < m_grid.GetRElements(); k++) {
			m_dataSurfacePressure[i][j] +=
				phys.GetG()
				* dataNode[RIx][k][i][j]
				* (m_dataZInterfaces[k+1][i][j]
					- m_dataZInterfaces[k][i][j]);
		}
	}
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatch::ComputeTemperature(
	int iDataIndex,
	DataLocation loc
) {
	const PhysicalConstants & phys = m_grid.GetModel().GetPhysicalConstants();

	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int PIx = 2;
	const int WIx = 3;
	const int RIx = 4;

	if (m_grid.GetModel().GetEquationSet().GetComponents() < 5) {
		_EXCEPTIONT("Invalid EquationSet.");
	}

	int k;
	int i;
	int j;

	// Calculate temperature on nodes
	if (loc == DataLocation_Node) {
		if (m_grid.GetVarLocation(PIx) == DataLocation_REdge) {
			InterpolateREdgeToNode(PIx, iDataIndex);
		}
		if (m_grid.GetVarLocation(RIx) == DataLocation_REdge) {
			InterpolateREdgeToNode(RIx, iDataIndex);
		}

		const DataArray4D<double> & dataNode = m_datavecStateNode[iDataIndex];

		for (k = 0; k < m_grid.GetRElements(); k++) {
		for (i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {

#ifdef FORMULATION_PRESSURE
			double dPressure = dataNode[PIx][k][i][j];
#endif
#if defined(FORMULATION_RHOTHETA_PI) || defined(FORMULATION_RHOTHETA_P)
			double dPressure = phys.PressureFromRhoTheta(dataNode[PIx][k][i][j]);
#endif
#if defined(FORMULATION_THETA) || defined(FORMULATION_THETA_FLUX)
			double dPressure =
				phys.PressureFromRhoTheta(
					dataNode[RIx][k][i][j] * dataNode[PIx][k][i][j]);
#endif
			m_dataTemperature[k][i][j] =
				dPressure / (dataNode[RIx][k][i][j] * phys.GetR());
		}
		}
		}
	}

	// Calculate temperature on interfaces
	if (loc == DataLocation_REdge) {
		_EXCEPTIONT("Temperature not implemented on interfaces");

		if (m_grid.GetVarLocation(PIx) == DataLocation_Node) {
			InterpolateNodeToREdge(PIx, iDataIndex);
		}
		if (m_grid.GetVarLocation(RIx) == DataLocation_Node) {
			InterpolateNodeToREdge(RIx, iDataIndex);
		}

		const DataArray4D<double> & dataNode = m_datavecStateREdge[iDataIndex];

		for (k = 0; k <= m_grid.GetRElements(); k++) {
		for (i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {

#ifdef FORMULATION_PRESSURE
			double dPressure = dataNode[PIx][k][i][j];
#endif
#if defined(FORMULATION_RHOTHETA_PI) || defined(FORMULATION_RHOTHETA_P)
			double dPressure = phys.PressureFromRhoTheta(dataNode[PIx][k][i][j]);
#endif
#if defined(FORMULATION_THETA) || defined(FORMULATION_THETA_FLUX)
			double dPressure =
				phys.PressureFromRhoTheta(
					dataNode[RIx][k][i][j] * dataNode[PIx][k][i][j]);
#endif

			m_dataTemperature[k][i][j] =
				dPressure / (dataNode[RIx][k][i][j] * phys.GetR());
		}
		}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatch::Checksum(
	DataType eDataType,
	DataArray1D<double> & dChecksums,
	int iDataIndex,
	ChecksumType eChecksumType
) const {
	int c;
	int i;
	int j;
	int k;

	// Verify consistency in number of components
	if (!m_fContainsData) {
		_EXCEPTIONT("Checksum called on uninitialized GridPatch");
	}

	// State variables
	DataArray4D<double> const * pDataNode;
	DataArray4D<double> const * pDataREdge;

	std::vector<int> nodevars;
	std::vector<int> redgevars;

	// State data
	if (eDataType == DataType_State) {
		pDataNode  = &(m_datavecStateNode[iDataIndex]);
		pDataREdge = &(m_datavecStateREdge[iDataIndex]);

		// Variables on nodes
		int nComponents = m_grid.GetModel().GetEquationSet().GetComponents();
		for (int c = 0; c < nComponents; c++) {
			if (m_grid.GetVarLocation(c) == DataLocation_Node) {
				nodevars.push_back(c);
			} else if (m_grid.GetVarLocation(c) == DataLocation_REdge) {
				redgevars.push_back(c);
			} else {
				_EXCEPTIONT("Not implemented.");
			}
		}

		if (dChecksums.GetRows() < nComponents) {
			_EXCEPTIONT("Invalid Checksum count");
		}

	// Tracer data
	} else if (eDataType == DataType_Tracers) {
		pDataNode  = &(m_datavecTracers[iDataIndex]);
		pDataREdge = NULL;

		int nTracers = m_grid.GetModel().GetEquationSet().GetTracers();
		for (int c = 0; c < nTracers; c++) {
			nodevars.push_back(c);
		}

		if (dChecksums.GetRows() < nTracers) {
			_EXCEPTIONT("Invalid Checksum count");
		}

	} else {
		_EXCEPTIONT("Invalid DataType.");
	}

	// ChecksumType_Sum
	if (eChecksumType == ChecksumType_Sum) {
		for (c = 0; c < nodevars.size(); c++) {
		for (k = 0; k < m_grid.GetRElements(); k++) {
		for (i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {
			dChecksums[nodevars[c]] +=
				  (*pDataNode)[nodevars[c]][k][i][j]
				* m_dataElementArea[k][i][j];
		}
		}
		}
		}

		for (c = 0; c < redgevars.size(); c++) {
		for (k = 0; k <= m_grid.GetRElements(); k++) {
		for (i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {
			dChecksums[redgevars[c]] +=
				  (*pDataREdge)[redgevars[c]][k][i][j]
				* m_dataElementAreaREdge[k][i][j];
		}
		}
		}
		}

	// ChecksumType_L1
	} else if (eChecksumType == ChecksumType_L1) {
		for (c = 0; c < nodevars.size(); c++) {
		for (k = 0; k < m_grid.GetRElements(); k++) {
		for (i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {
			double dValue = fabs((*pDataNode)[nodevars[c]][k][i][j]);
			dChecksums[nodevars[c]] +=
				dValue * m_dataElementArea[k][i][j];
		}
		}
		}
		}

		for (c = 0; c < redgevars.size(); c++) {
		for (k = 0; k <= m_grid.GetRElements(); k++) {
		for (i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {
			double dValue = fabs((*pDataREdge)[redgevars[c]][k][i][j]);
			dChecksums[redgevars[c]] +=
				dValue * m_dataElementAreaREdge[k][i][j];
		}
		}
		}
		}

	// ChecksumType_L2
	} else if (eChecksumType == ChecksumType_L2) {
		for (c = 0; c < nodevars.size(); c++) {
		for (k = 0; k < m_grid.GetRElements(); k++) {
		for (i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {
			double dValue = (*pDataNode)[nodevars[c]][k][i][j];
			dChecksums[nodevars[c]] +=
				dValue * dValue * m_dataElementArea[k][i][j];
		}
		}
		}
		}

		for (c = 0; c < redgevars.size(); c++) {
		for (k = 0; k <= m_grid.GetRElements(); k++) {
		for (i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {
			double dValue = (*pDataREdge)[redgevars[c]][k][i][j];
			dChecksums[redgevars[c]] +=
				dValue * dValue * m_dataElementAreaREdge[k][i][j];
		}
		}
		}
		}

	// ChecksumType_Linf
	} else if (eChecksumType == ChecksumType_Linf) {
		for (c = 0; c < nodevars.size(); c++) {
		for (k = 0; k < m_grid.GetRElements(); k++) {
		for (i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {
			double dValue = (*pDataNode)[nodevars[c]][k][i][j];
			if (fabs(dValue) > dChecksums[nodevars[c]]) {
				dChecksums[nodevars[c]] = fabs(dValue);
			}
		}
		}
		}
		}

		for (c = 0; c < redgevars.size(); c++) {
		for (k = 0; k <= m_grid.GetRElements(); k++) {
		for (i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {
			double dValue = (*pDataREdge)[redgevars[c]][k][i][j];
			if (fabs(dValue) > dChecksums[redgevars[c]]) {
				dChecksums[redgevars[c]] = fabs(dValue);
			}
		}
		}
		}
		}

	// Invalid data type
	} else {
		_EXCEPTIONT("Invalid DataType in Checksum: Expected State or Tracers");
	}
}

///////////////////////////////////////////////////////////////////////////////

double GridPatch::ComputeTotalEnergy(
	int iDataIndex
) const {
	// Physical constants
	const PhysicalConstants & phys = m_grid.GetModel().GetPhysicalConstants();

	// Accumulated local energy
	double dLocalEnergy = 0.0;

	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int HIx = 2;
	const int PIx = 2;
	const int WIx = 3;
	const int RIx = 4;

	// Determine type of energy to compute from EquationSet
	EquationSet::Type eEquationSetType =
		m_grid.GetModel().GetEquationSet().GetType();

	// Grid data
	if ((iDataIndex < 0) || (iDataIndex >= m_datavecStateNode.size())) {
		_EXCEPTION1("iDataIndex out of range: %i", iDataIndex);
	}
	const DataArray4D<double> & dataNode = m_datavecStateNode[iDataIndex];

	// Shallow Water Energy
	if (eEquationSetType == EquationSet::ShallowWaterEquations) {

		// Loop over all elements
		int k;
		int i;
		int j;

		for (k = 0; k < m_grid.GetRElements(); k++) {
		for (i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {
			double dUdotU =
				+ m_dataContraMetric2DB[i][j][1]
					* dataNode[UIx][k][i][j] * dataNode[UIx][k][i][j]
				- 2.0 * m_dataContraMetric2DA[i][j][1]
					* dataNode[UIx][k][i][j] * dataNode[VIx][k][i][j]
				+ m_dataContraMetric2DA[i][j][0]
					* dataNode[VIx][k][i][j] * dataNode[VIx][k][i][j];

			dUdotU *= m_dataJacobian2D[i][j] * m_dataJacobian2D[i][j];

			double dKineticEnergy =
				0.5 * (dataNode[HIx][k][i][j] - m_dataTopography[i][j])
					* dUdotU;

			double dPotentialEnergy =
				0.5 * phys.GetG()
					* (dataNode[HIx][k][i][j] * dataNode[HIx][k][i][j]
						- m_dataTopography[i][j] * m_dataTopography[i][j]);

			dLocalEnergy += m_dataElementArea[k][i][j]
				* (dKineticEnergy + dPotentialEnergy);
		}
		}
		}

	} else {

		// Loop over all elements
		int k;
		int i;
		int j;
/*
		double dTotalKineticEnergy = 0.0;
		double dTotalInternalEnergy = 0.0;
		double dTotalPotentialEnergy = 0.0;
*/
		for (k = 0; k < m_grid.GetRElements(); k++) {
		//for (k = 0; k < 1; k++) {
		for (i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {
/*
			double dUa = dataNode[UIx][k][i][j];
			double dUb = dataNode[VIx][k][i][j];
			double dUx = dataNode[WIx][k][i][j];

			double dCovUa =
				  m_dataCovMetricA[k][i][j][0] * dUa
				+ m_dataCovMetricA[k][i][j][1] * dUb
				+ m_dataCovMetricA[k][i][j][2] * dUx;

			double dCovUb =
				  m_dataCovMetricB[k][i][j][0] * dUa
				+ m_dataCovMetricB[k][i][j][1] * dUb
				+ m_dataCovMetricB[k][i][j][2] * dUx;

			double dCovUx =
				  m_dataCovMetricXi[k][i][j][0] * dUa
				+ m_dataCovMetricXi[k][i][j][1] * dUb
				+ m_dataCovMetricXi[k][i][j][2] * dUx;

			double dUdotU =
				dCovUa * dUa + dCovUb * dUb + dCovUx * dUx;
*/

			double dCovUa = dataNode[UIx][k][i][j];
			double dCovUb = dataNode[VIx][k][i][j];
			double dCovUx = dataNode[WIx][k][i][j] * m_dataDerivRNode[k][i][j][2];

			double dConUa =
				  m_dataContraMetricA[k][i][j][0] * dCovUa
				+ m_dataContraMetricA[k][i][j][1] * dCovUb
				+ m_dataContraMetricA[k][i][j][2] * dCovUx;

			double dConUb =
				  m_dataContraMetricB[k][i][j][0] * dCovUa
				+ m_dataContraMetricB[k][i][j][1] * dCovUb
				+ m_dataContraMetricB[k][i][j][2] * dCovUx;

			double dConUx =
				  m_dataContraMetricXi[k][i][j][0] * dCovUa
				+ m_dataContraMetricXi[k][i][j][1] * dCovUb
				+ m_dataContraMetricXi[k][i][j][2] * dCovUx;

			double dUdotU = dConUa * dCovUa + dConUb * dCovUb + dConUx * dCovUx;

			double dKineticEnergy =
				0.5 * dataNode[RIx][k][i][j] * dUdotU;

#ifdef FORMULATION_PRESSURE
			double dPressure = dataNode[PIx][k][i][j];
#endif
#if defined(FORMULATION_RHOTHETA_PI) || defined(FORMULATION_RHOTHETA_P)
			double dPressure = phys.PressureFromRhoTheta(dataNode[PIx][k][i][j]);
#endif
#if defined(FORMULATION_THETA) || defined(FORMULATION_THETA_FLUX)
			double dPressure =
				phys.PressureFromRhoTheta(
					dataNode[RIx][k][i][j] * dataNode[PIx][k][i][j]);
#endif

			double dInternalEnergy =
				dPressure / (phys.GetGamma() - 1.0);

			double dPotentialEnergy =
				phys.GetG() * dataNode[RIx][k][i][j] * m_dataZLevels[k][i][j];

			dLocalEnergy += m_dataElementArea[k][i][j]
				* (dKineticEnergy + dInternalEnergy + dPotentialEnergy);
/*
			dTotalKineticEnergy += m_dataElementArea[k][i][j] * dKineticEnergy;
			dTotalInternalEnergy += m_dataElementArea[k][i][j] * dInternalEnergy;
			dTotalPotentialEnergy += m_dataElementArea[k][i][j] * dPotentialEnergy;
*/
		}
		}
		}
/*
		printf("%1.15e %1.15e %1.15e\n",
			dTotalKineticEnergy,
			dTotalInternalEnergy,
			dTotalPotentialEnergy);
*/
	}

	return dLocalEnergy;
}

///////////////////////////////////////////////////////////////////////////////

double GridPatch::ComputeTotalPotentialEnstrophy(
	int iDataIndex
) {
	// Physical constants
	const PhysicalConstants & phys = m_grid.GetModel().GetPhysicalConstants();

	// Accumulated local energy
	double dLocalPotentialEnstrophy = 0.0;

	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int HIx = 2;
	const int PIx = 2;
	const int WIx = 3;
	const int RIx = 4;

	// Determine type of energy to compute from EquationSet
	EquationSet::Type eEquationSetType =
		m_grid.GetModel().GetEquationSet().GetType();

	// Grid data
	if ((iDataIndex < 0) || (iDataIndex >= m_datavecStateNode.size())) {
		_EXCEPTION1("iDataIndex out of range: %i", iDataIndex);
	}
	const DataArray4D<double> & dataNode = m_datavecStateNode[iDataIndex];

	// Shallow Water PotentialEnstrophy
	if (eEquationSetType == EquationSet::ShallowWaterEquations) {

		// Loop over all elements
		int k;
		int i;
		int j;

		for (k = 0; k < m_grid.GetRElements(); k++) {
		for (i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {
			double dPlanetaryVorticity =
				2.0 * phys.GetOmega() * sin(m_dataLat[i][j]);

			double dAbsoluteVorticity =
				m_dataVorticity[k][i][j] + dPlanetaryVorticity;

			dLocalPotentialEnstrophy +=
				m_dataElementArea[k][i][j]
					* 0.5 * dAbsoluteVorticity * dAbsoluteVorticity
					/ (dataNode[HIx][k][i][j] - m_dataTopography[i][j]);
		}
		}
		}

	} else {

#pragma message "Calculate total potential enstrophy here"

		// Set to zero
		dLocalPotentialEnstrophy = 0.0;

		// Loop over all elements
		int k;
		int i;
		int j;

		for (k = 0; k < m_grid.GetRElements(); k++) {
		for (i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {

			// Zonal momentum
			dLocalPotentialEnstrophy +=
				m_dataElementArea[k][i][j]
				* dataNode[RIx][k][i][j]
				* dataNode[UIx][k][i][j];
		}
		}
		}

	}

	return dLocalPotentialEnstrophy;
}

///////////////////////////////////////////////////////////////////////////////

double GridPatch::ComputeTotalVerticalMomentum(
	int iDataIndex
) {
	// Physical constants
	const PhysicalConstants & phys = m_grid.GetModel().GetPhysicalConstants();

	// Accumulated local energy
	double dLocalVerticalMomentum = 0.0;

	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int PIx = 2;
	const int WIx = 3;
	const int RIx = 4;

	// Determine type of energy to compute from EquationSet
	EquationSet::Type eEquationSetType =
		m_grid.GetModel().GetEquationSet().GetType();

	// Grid data
	if ((iDataIndex < 0) || (iDataIndex >= m_datavecStateNode.size())) {
		_EXCEPTION1("iDataIndex out of range: %i", iDataIndex);
	}
	const DataArray4D<double> & dataNode = m_datavecStateNode[iDataIndex];

	// Shallow Water ZonalMomentum
	if (eEquationSetType == EquationSet::ShallowWaterEquations) {
		_EXCEPTIONT("ComputeTotalVerticalMomentum() Not implemented "
			"for ShallowWaterEquations");

	} else {

		// Set to zero
		dLocalVerticalMomentum = 0.0;

		// Loop over all elements
		int k;
		int i;
		int j;

		for (k = 0; k < m_grid.GetRElements(); k++) {
		for (i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {

			// Vertical momentum
			dLocalVerticalMomentum +=
				m_dataElementArea[k][i][j]
				* dataNode[RIx][k][i][j]
				* dataNode[WIx][k][i][j];
		}
		}
		}

	}

	return dLocalVerticalMomentum;
}

///////////////////////////////////////////////////////////////////////////////

void GridPatch::PrepareExchange() {
	m_connect.PrepareExchange();
}

///////////////////////////////////////////////////////////////////////////////

void GridPatch::Send(
	DataType eDataType,
	int iDataIndex
) {
	// State data
	if (eDataType == DataType_State) {
		if ((iDataIndex < 0) || (iDataIndex > m_datavecStateNode.size())) {
			_EXCEPTIONT("Invalid state data instance.");
		}

		m_connect.Pack(m_datavecStateNode[iDataIndex]);
		m_connect.Pack(m_datavecStateREdge[iDataIndex]);
		m_connect.Send();

	// Tracer data
	} else if (eDataType == DataType_Tracers) {
		if ((iDataIndex < 0) || (iDataIndex > m_datavecTracers.size())) {
			_EXCEPTIONT("Invalid tracers data instance.");
		}

		m_connect.Pack(m_datavecTracers[iDataIndex]);
		m_connect.Send();

	// Vorticity data
	} else if (eDataType == DataType_Vorticity) {
		m_connect.Pack(m_dataVorticity);
		m_connect.Send();

	// Divergence data
	} else if (eDataType == DataType_Divergence) {
		m_connect.Pack(m_dataDivergence);
		m_connect.Send();

	// Temperature data
	} else if (eDataType == DataType_Temperature) {
		m_connect.Pack(m_dataTemperature);
		m_connect.Send();

	// Topography derivative data
	} else if (eDataType == DataType_TopographyDeriv) {
		m_connect.Pack(m_dataTopographyDeriv);
		m_connect.Send();

	// Invalid data
	} else {
		_EXCEPTIONT("Invalid DataType");
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatch::Receive(
	DataType eDataType,
	int iDataIndex
) {
	// State data
	if (eDataType == DataType_State) {
		if ((iDataIndex < 0) || (iDataIndex > m_datavecStateNode.size())) {
			_EXCEPTIONT("Invalid state data instance.");
		}

		Neighbor * pNeighbor;
		while ((pNeighbor = m_connect.WaitReceive()) != NULL) {
			pNeighbor->Unpack(m_grid, m_datavecStateNode[iDataIndex]);
			pNeighbor->Unpack(m_grid, m_datavecStateREdge[iDataIndex]);
		}

	// Tracer data
	} else if (eDataType == DataType_Tracers) {
		if ((iDataIndex < 0) || (iDataIndex > m_datavecTracers.size())) {
			_EXCEPTIONT("Invalid tracers data instance.");
		}

		Neighbor * pNeighbor;
		while ((pNeighbor = m_connect.WaitReceive()) != NULL) {
			pNeighbor->Unpack(m_grid, m_datavecTracers[iDataIndex]);
		}

	// Vorticity data
	} else if (eDataType == DataType_Vorticity) {
		Neighbor * pNeighbor;
		while ((pNeighbor = m_connect.WaitReceive()) != NULL) {
			pNeighbor->Unpack(m_dataVorticity);
		}

	// Divergence data
	} else if (eDataType == DataType_Divergence) {
		Neighbor * pNeighbor;
		while ((pNeighbor = m_connect.WaitReceive()) != NULL) {
			pNeighbor->Unpack(m_dataDivergence);
		}

	// Temperature data
	} else if (eDataType == DataType_Temperature) {
		Neighbor * pNeighbor;
		while ((pNeighbor = m_connect.WaitReceive()) != NULL) {
			pNeighbor->Unpack(m_dataTemperature);
		}

	// Topographic derivatives
	} else if (eDataType == DataType_TopographyDeriv) {
		Neighbor * pNeighbor;
		while ((pNeighbor = m_connect.WaitReceive()) != NULL) {
			pNeighbor->Unpack(m_dataTopographyDeriv);
		}

	// Invalid data
	} else {
		_EXCEPTIONT("Invalid DataType");
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatch::SendBuffers() {
	m_connect.SendBuffers();
}

///////////////////////////////////////////////////////////////////////////////

void GridPatch::ReceiveBuffers() {
	Neighbor * pNeighbor;
	while ((pNeighbor = m_connect.WaitReceive()) != NULL);
}

///////////////////////////////////////////////////////////////////////////////

void GridPatch::CompleteExchange() {
	m_connect.WaitSend();
}

///////////////////////////////////////////////////////////////////////////////

void GridPatch::CopyData(
	int ixSource,
	int ixDest,
	DataType eDataType
) {
	// Copy over State data
	if (eDataType == DataType_State) {
		if ((ixSource < 0) || (ixSource >= m_datavecStateNode.size())) {
			_EXCEPTIONT("Invalid ixSource index in CopyData.");
		}
		if ((ixDest < 0) || (ixDest >= m_datavecStateNode.size())) {
			_EXCEPTIONT("Invalid ixDest index in CopyData.");
		}

		m_datavecStateNode[ixDest]  = m_datavecStateNode[ixSource];
		m_datavecStateREdge[ixDest] = m_datavecStateREdge[ixSource];

	// Copy over Tracers data
	} else if (eDataType == DataType_Tracers) {
		if ((ixSource < 0) || (ixSource >= m_datavecTracers.size())) {
			_EXCEPTIONT("Invalid ixSource index in CopyData.");
		}
		if ((ixDest < 0) || (ixDest >= m_datavecTracers.size())) {
			_EXCEPTIONT("Invalid ixDest index in CopyData.");
		}

		m_datavecTracers[ixDest] = m_datavecTracers[ixSource];

	// Invalid datatype; only State or Tracers expected
	} else {
		_EXCEPTIONT("Invalid DataType specified for CopyData.");
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatch::LinearCombineData(
	const DataArray1D<double> & dCoeff,
	int ixDest,
	DataType eDataType
) {
	// Check bounds on Coeff array
	if (ixDest >= dCoeff.GetRows()) {
		_EXCEPTION2("Destination index [%i] out of coefficient bounds [0,%i]",
			ixDest, dCoeff.GetRows()-1);
	}

	// Check bounds on ixDest for State data
	if (eDataType == DataType_State) {
		if ((ixDest < 0) || (ixDest >= m_datavecStateNode.size())) {
			_EXCEPTIONT("Invalid ixDest index in LinearCombineData.");
		}
		if (dCoeff.GetRows() > m_datavecStateNode.size()) {
			_EXCEPTIONT("Too many elements in coefficient vector.");
		}

		// Premultiply
		if (dCoeff[ixDest] == 0.0) {
			m_datavecStateNode [ixDest].Zero();
			m_datavecStateREdge[ixDest].Zero();
		} else {
			m_datavecStateNode [ixDest].Scale(dCoeff[ixDest]);
			m_datavecStateREdge[ixDest].Scale(dCoeff[ixDest]);
		}

		// Consider all other terms
		for (int m = 0; m < dCoeff.GetRows(); m++) {
			if (m == ixDest) {
				continue;
			}
			if (dCoeff[m] == 0.0) {
				continue;
			}

			m_datavecStateNode[ixDest].AddProduct(
				m_datavecStateNode[m], dCoeff[m]);
			m_datavecStateREdge[ixDest].AddProduct(
				m_datavecStateREdge[m], dCoeff[m]);
		}

	// Check bounds on ixDest for Tracers data
	} else if (eDataType == DataType_Tracers) {
		if ((ixDest < 0) || (ixDest >= m_datavecTracers.size())) {
			_EXCEPTIONT("Invalid ixDest index in LinearCombineData.");
		}
		if (dCoeff.GetRows() > m_datavecTracers.size()) {
			_EXCEPTIONT("Too many elements in coefficient vector.");
		}

		// If tracers are not initialized, do nothing
		if (!m_datavecTracers[ixDest].IsAttached()) {
			return;
		}

		// Premultiply
		if (dCoeff[ixDest] == 0.0) {
			m_datavecTracers[ixDest].Zero();
		} else {
			m_datavecTracers[ixDest].Scale(dCoeff[ixDest]);
		}

		// Consider all other terms
		for (int m = 0; m < dCoeff.GetRows(); m++) {
			if (m == ixDest) {
				continue;
			}
			if (dCoeff[m] == 0.0) {
				continue;
			}

			m_datavecTracers[ixDest].AddProduct(
				m_datavecTracers[m], dCoeff[m]);
		}

	// Invalid datatype; only State or Tracers expected
	} else {
		_EXCEPTIONT("Invalid DataType specified for LinearCombineData.");
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatch::ZeroData(
	int ixData,
	DataType eDataType
) {
	// Check bounds on ixDest for State data
	if (eDataType == DataType_State) {
		if ((ixData < 0) || (ixData >= m_datavecStateNode.size())) {
			_EXCEPTIONT("Invalid ixData index in LinearCombineData.");
		}

		m_datavecStateNode [ixData].Zero();
		m_datavecStateREdge[ixData].Zero();

	// Check bounds on ixDest for Tracers data
	} else if (eDataType == DataType_Tracers) {
		if ((ixData < 0) || (ixData >= m_datavecTracers.size())) {
			_EXCEPTIONT("Invalid ixData index in LinearCombineData.");
		}

		if (!m_datavecTracers[ixData].IsAttached()) {
			return;
		}

		m_datavecTracers[ixData].Zero();

	// Invalid datatype; only State or Tracers expected
	} else {
		_EXCEPTIONT("Invalid DataType specified for LinearCombineData.");
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatch::InterpolateNodeToREdge(
	int iVar,
	int iDataIndex
) {
	_EXCEPTIONT("Not implemented.");
}

///////////////////////////////////////////////////////////////////////////////

void GridPatch::InterpolateREdgeToNode(
	int iVar,
	int iDataIndex
) {
	_EXCEPTIONT("Not implemented.");
}

///////////////////////////////////////////////////////////////////////////////

void GridPatch::AddReferenceState(
	int ix
) {
	if ((ix < 0) || (ix >= m_datavecStateNode.size())) {
		_EXCEPTIONT("Invalid ix in AddReferenceState.");
	}

	m_datavecStateNode[ix].AddProduct(m_dataRefStateNode, 1.0);
}

///////////////////////////////////////////////////////////////////////////////

void GridPatch::InterpolateData(
	DataType eDataType,
	const DataArray1D<double> & dREta,
	const DataArray1D<double> & dAlpha,
	const DataArray1D<double> & dBeta,
	const DataArray1D<int> & iPatch,
	DataArray3D<double> & dInterpData,
	DataLocation eOnlyVariablesAt,
	bool fIncludeReferenceState,
	bool fConvertToPrimitive
) {
	_EXCEPTIONT("Unimplemented.");
}

///////////////////////////////////////////////////////////////////////////////

