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

#ifdef TEMPEST_MPIOMP
#include <mpi.h>
#endif

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
	m_iProcessor = 0;
#ifdef TEMPEST_MPIOMP
	MPI_Comm_rank(MPI_COMM_WORLD, &m_iProcessor);
#endif

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
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_grid.GetRElements());

	// Jacobian at each interface
	m_dataJacobianREdge.SetSize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_grid.GetRElements()+1);

	// Contravariant metric components at each node
	m_dataContraMetricA.SetSize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_grid.GetRElements(),
		3);

	m_dataContraMetricB.SetSize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_grid.GetRElements(),
		3);

	m_dataContraMetricXi.SetSize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_grid.GetRElements(),
		3);

	// Xi contravariant metric on interfaces
	m_dataContraMetricAREdge.SetSize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_grid.GetRElements()+1,
		3);

	m_dataContraMetricBREdge.SetSize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_grid.GetRElements()+1,
		3);

	m_dataContraMetricXiREdge.SetSize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_grid.GetRElements()+1,
		3);

	// Vertical coordinate transform (derivatives of the radius)
	// include transformation from xi to z/r dxi/dr
	m_dataDerivRNode.SetSize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_grid.GetRElements(),
		4);

	m_dataDerivRREdge.SetSize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_grid.GetRElements()+1,
		4);

	// Element area at each node
	m_dataElementAreaNode.SetSize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_grid.GetRElements());

	// Element area at each interface
	m_dataElementAreaREdge.SetSize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_grid.GetRElements()+1);

	// Topography height at each node
	m_dataTopography.SetSize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth());

	// Topography derivatives at each node
	m_dataTopographyDeriv.SetDataType(DataType_TopographyDeriv);
	m_dataTopographyDeriv.SetSize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		2);

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
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_grid.GetRElements());

	// Radial coordinate at each interface
	m_dataZInterfaces.SetSize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_grid.GetRElements()+1);

	// Active State Patch Index
	m_iActiveStatePatchIx.SetSize(1);

	// Initialize reference state
	m_dataRefStateNode.SetDataType(DataType_State);
	m_dataRefStateNode.SetDataLocation(DataLocation_Node);
	m_dataRefStateNode.SetSize(
		eqn.GetComponents(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_grid.GetRElements());

	m_dataRefStateREdge.SetDataType(DataType_State);
	m_dataRefStateREdge.SetDataLocation(DataLocation_REdge);
	m_dataRefStateREdge.SetSize(
		eqn.GetComponents(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_grid.GetRElements()+1);

	m_dataRefTracers.SetDataType(DataType_Tracers);
	m_dataRefTracers.SetDataLocation(DataLocation_Node);
	m_dataRefTracers.SetSize(
		eqn.GetTracers(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_grid.GetRElements());

	// Rayleigh friction strength
	m_dataRayleighStrengthNode.SetDataType(DataType_None);
	m_dataRayleighStrengthNode.SetDataLocation(DataLocation_Node);
	m_dataRayleighStrengthNode.SetSize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_grid.GetRElements());

	// Rayleigh friction strength
	m_dataRayleighStrengthREdge.SetDataType(DataType_None);
	m_dataRayleighStrengthREdge.SetDataLocation(DataLocation_REdge);
	m_dataRayleighStrengthREdge.SetSize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_grid.GetRElements()+1);

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
/*
	m_dcGeometric.PushDataChunk(&m_dataCovMetricA);
	m_dcGeometric.PushDataChunk(&m_dataCovMetricB);
	m_dcGeometric.PushDataChunk(&m_dataCovMetricXi);
*/
	m_dcGeometric.PushDataChunk(&m_dataContraMetricAREdge);
	m_dcGeometric.PushDataChunk(&m_dataContraMetricBREdge);
	m_dcGeometric.PushDataChunk(&m_dataContraMetricXiREdge);
	m_dcGeometric.PushDataChunk(&m_dataDerivRNode);
	m_dcGeometric.PushDataChunk(&m_dataDerivRREdge);
	m_dcGeometric.PushDataChunk(&m_dataElementAreaNode);
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
			m_box.GetATotalWidth(),
			m_box.GetBTotalWidth(),
			m_grid.GetRElements());

		m_datavecStateREdge[m].SetDataType(DataType_State);
		m_datavecStateREdge[m].SetDataLocation(DataLocation_REdge);
		m_datavecStateREdge[m].SetSize(
			eqn.GetComponents(),
			m_box.GetATotalWidth(),
			m_box.GetBTotalWidth(),
			m_grid.GetRElements()+1);
	}

	m_dcActiveState.PushDataChunk(&m_iActiveStatePatchIx);
	m_dcActiveState.PushDataChunk(&m_datavecStateNode[0]);
	m_dcActiveState.PushDataChunk(&m_datavecStateREdge[0]);

	for (int m = 1; m < model.GetComponentDataInstances(); m++) {
  		m_dcBufferState.PushDataChunk(&m_datavecStateNode[m]);
    		m_dcBufferState.PushDataChunk(&m_datavecStateREdge[m]);
  	}

	m_dataXiDiffNode.SetDataType(DataType_None);
	m_dataXiDiffNode.SetDataLocation(DataLocation_Node);
	m_dataXiDiffNode.SetSize(
		2,
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_grid.GetRElements());
	m_dcActiveState.PushDataChunk(&m_dataXiDiffNode);

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
				m_box.GetATotalWidth(),
				m_box.GetBTotalWidth(),
				m_grid.GetRElements());
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
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_grid.GetRElements());

	m_dataDxPressure.SetDataType(DataType_Pressure);
	m_dataDxPressure.SetDataLocation(DataLocation_Node);
	m_dataDxPressure.SetSize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_grid.GetRElements());

	// Vorticity data
	m_dataVorticity.SetDataType(DataType_Vorticity);
	m_dataVorticity.SetDataLocation(DataLocation_Node);
	m_dataVorticity.SetSize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_grid.GetRElements());

	// Divergence data
	m_dataDivergence.SetDataType(DataType_Divergence);
	m_dataDivergence.SetDataLocation(DataLocation_Node);
	m_dataDivergence.SetSize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_grid.GetRElements());

	// Temperature data
	m_dataTemperature.SetDataType(DataType_Temperature);
	m_dataTemperature.SetDataLocation(DataLocation_Node);
	m_dataTemperature.SetSize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_grid.GetRElements());

	// Richardson number data
	m_dataRichardson.SetDataType(DataType_Richardson);
	m_dataRichardson.SetDataLocation(DataLocation_Node);
	m_dataRichardson.SetSize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_grid.GetRElements());

	// Convective stability number data
	m_dataConvective.SetDataType(DataType_Convective);
	m_dataConvective.SetDataLocation(DataLocation_Node);
	m_dataConvective.SetSize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_grid.GetRElements());

	// Zonal contravariant momentum tendencty data
        m_dataZonalForce.SetDataType(DataType_ZonalForce);
        m_dataZonalForce.SetDataLocation(DataLocation_Node);
        m_dataZonalForce.SetSize(
                m_box.GetATotalWidth(),
                m_box.GetBTotalWidth(),
                m_grid.GetRElements());
	
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
	m_dcAuxiliary.PushDataChunk(&m_dataRichardson);
	m_dcAuxiliary.PushDataChunk(&m_dataConvective);
	m_dcAuxiliary.PushDataChunk(&m_dataZonalForce);
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
	m_dataElementAreaNode.Deallocate();
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

#if defined(FIXED_RELEMENTS)
	const int nRElements = FIXED_RELEMENTS;
	if (nRElements != m_grid.GetRElements()) {
		_EXCEPTIONT("Command line levels must match FIXED_RELEMENTS");
	}
#else
	const int nRElements = m_grid.GetRElements();
#endif

	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int PIx = 2;
	const int WIx = 3;
	const int RIx = 4;

	if (m_grid.GetModel().GetEquationSet().GetComponents() < 5) {
		_EXCEPTIONT("Invalid EquationSet.");
	}

	// Calculate hydrostatic surface pressure on nodes
	const DataArray4D<double> & dataNode = m_datavecStateNode[iDataIndex];

	for (int i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
	for (int j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {

		m_dataSurfacePressure(i,j) = 0.0;

		for (int k = 0; k < nRElements; k++) {
			m_dataSurfacePressure(i,j) +=
				phys.GetG()
				* dataNode(RIx,i,j,k)
				* (m_dataZInterfaces(i,j,k+1)
					- m_dataZInterfaces(i,j,k));
		}
	}
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatch::ComputePressure(
	int iDataIndex
) {
	const PhysicalConstants & phys = m_grid.GetModel().GetPhysicalConstants();

#if defined(FIXED_RELEMENTS)
	const int nRElements = FIXED_RELEMENTS;
	if (nRElements != m_grid.GetRElements()) {
		_EXCEPTIONT("Command line levels must match FIXED_RELEMENTS");
	}
#else
	const int nRElements = m_grid.GetRElements();
#endif

	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int PIx = 2;
	const int WIx = 3;
	const int RIx = 4;

	if (m_grid.GetModel().GetEquationSet().GetComponents() < 5) {
		_EXCEPTIONT("Invalid EquationSet.");
	}

	const DataArray4D<double> & dataNode = m_datavecStateNode[iDataIndex];

	for (int i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
	for (int j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {
#pragma simd
	for (int k = 0; k < nRElements; k++) {
#ifdef FORMULATION_PRESSURE
		m_dataPressure(i,j,k) =
			dataNode(PIx,i,j,k);
#endif
#if defined(FORMULATION_RHOTHETA_PI) || defined(FORMULATION_RHOTHETA_P)
		m_dataPressure(i,j,k) =
			phys.PressureFromRhoTheta(
				dataNode(PIx,i,j,k));
#endif
#if defined(FORMULATION_THETA) || defined(FORMULATION_THETA_FLUX)
		m_dataPressure(i,j,k) =
			phys.PressureFromRhoTheta(
				dataNode(RIx,i,j,k)
				* dataNode(PIx,i,j,k));
#endif
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

#if defined(FIXED_RELEMENTS)
	const int nRElements = FIXED_RELEMENTS;
	if (nRElements != m_grid.GetRElements()) {
		_EXCEPTIONT("Command line levels must match FIXED_RELEMENTS");
	}
#else
	const int nRElements = m_grid.GetRElements();
#endif

	// Indices of EquationSet variables
	const int UIx = 0;
	const int VIx = 1;
	const int PIx = 2;
	const int WIx = 3;
	const int RIx = 4;

	if (m_grid.GetModel().GetEquationSet().GetComponents() < 5) {
		_EXCEPTIONT("Invalid EquationSet.");
	}

	// Calculate temperature on nodes
	if (loc == DataLocation_Node) {
		if (m_grid.GetVarLocation(PIx) == DataLocation_REdge) {
			InterpolateREdgeToNode(PIx, iDataIndex);
		}
		if (m_grid.GetVarLocation(RIx) == DataLocation_REdge) {
			InterpolateREdgeToNode(RIx, iDataIndex);
		}

		const DataArray4D<double> & dataNode = m_datavecStateNode[iDataIndex];

		for (int i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (int j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {
#pragma simd
		for (int k = 0; k < nRElements; k++) {

#ifdef FORMULATION_PRESSURE
			const double dPressure = dataNode(PIx,i,j,k);
#endif
#if defined(FORMULATION_RHOTHETA_PI) || defined(FORMULATION_RHOTHETA_P)
			const double dPressure = phys.PressureFromRhoTheta(dataNode(PIx,i,j,k));
#endif
#if defined(FORMULATION_THETA) || defined(FORMULATION_THETA_FLUX)
			const double dPressure =
				phys.PressureFromRhoTheta(
					dataNode(RIx,i,j,k)
					* dataNode(PIx,i,j,k));
#endif
			m_dataTemperature(i,j,k) =
				dPressure / (dataNode(RIx,i,j,k) * phys.GetR());
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

		for (int i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (int j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {
#pragma simd
		for (int k = 0; k <= nRElements; k++) {

#ifdef FORMULATION_PRESSURE
			const double dPressure = dataNode(PIx,i,j,k);
#endif
#if defined(FORMULATION_RHOTHETA_PI) || defined(FORMULATION_RHOTHETA_P)
			const double dPressure =
				phys.PressureFromRhoTheta(dataNode(PIx,i,j,k));
#endif
#if defined(FORMULATION_THETA) || defined(FORMULATION_THETA_FLUX)
			const double dPressure =
				phys.PressureFromRhoTheta(
					dataNode(RIx,i,j,k)
					* dataNode(PIx,i,j,k));
#endif

			m_dataTemperature(i,j,k) =
				dPressure / (dataNode(RIx,i,j,k) * phys.GetR());
		}
		}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatch::ComputeZonalForce(
        int iDataIndex,
        DataLocation loc
) {
	const PhysicalConstants & phys = m_grid.GetModel().GetPhysicalConstants();

#if defined(FIXED_RELEMENTS)
        const int nRElements = FIXED_RELEMENTS;
        if (nRElements != m_grid.GetRElements()) {
                _EXCEPTIONT("Command line levels must match FIXED_RELEMENTS");
        }
#else
        const int nRElements = m_grid.GetRElements();
#endif

        // Indices of EquationSet variables
        const int UIx = 0;
        const int VIx = 1;
        const int PIx = 2;
        const int WIx = 3;
        const int RIx = 4;

        if (m_grid.GetModel().GetEquationSet().GetComponents() < 5) {
                _EXCEPTIONT("Invalid EquationSet.");
        }

	// Calculate zonal force on nodes
        if (loc == DataLocation_Node) {
                if (m_grid.GetVarLocation(UIx) == DataLocation_REdge) {
                        InterpolateREdgeToNode(UIx, iDataIndex);
                }
		if (m_grid.GetVarLocation(WIx) == DataLocation_REdge) {
                        InterpolateREdgeToNode(WIx, iDataIndex);
                }
                if (m_grid.GetVarLocation(RIx) == DataLocation_REdge) {
                        InterpolateREdgeToNode(RIx, iDataIndex);
                }

		// LHS tendency data is at the back instance of state vector
		const int iLHSdex = m_grid.GetModel().GetComponentDataInstances() - 1; 
                const DataArray4D<double> & dataNode = m_datavecStateNode[iDataIndex];
		const DataArray4D<double> & dataTendency = m_datavecStateNode[iLHSdex];

                double dUX = 0.0;
		double dUXdt = 0.0;
		for (int i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
                for (int j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {
                for (int k = 0; k < nRElements; k++) {
			const double dZda = m_dataDerivRNode(i,j,k,0);
                        const double dZdxi = m_dataDerivRNode(i,j,k,2);

			// Contravariant zonal velocity
			const double dCovUa = dataNode(UIx,i,j,k);
			const double dCovUx = dataNode(WIx,i,j,k);

			dUX - dCovUa - (dZda / dZdxi) * dCovUx;

			// Contravariant zonal velocity tendency
			const double dCovUadt = dataTendency(UIx,i,j,k);
                        const double dCovUxdt = dataTendency(WIx,i,j,k);

                        dUXdt = dCovUadt - (dZda / dZdxi) * dCovUxdt;

			// Compute the zonal momentum tendency by product rule
			m_dataZonalForce(i,j,k) = dataNode(RIx,i,j,k) * 
				dUXdt + dUX * dataTendency(RIx,i,j,k);

                }
                }
                }
        }

	// Calculate zonal force on interfaces
        if (loc == DataLocation_REdge) {
                _EXCEPTIONT("Zonal force not implemented on interfaces");
        }
}

///////////////////////////////////////////////////////////////////////////////

void GridPatch::Checksum(
	DataType eDataType,
	DataArray1D<double> & dChecksums,
	int iDataIndex,
	ChecksumType eChecksumType
) const {

#if defined(FIXED_RELEMENTS)
	const int nRElements = FIXED_RELEMENTS;
	if (nRElements != m_grid.GetRElements()) {
		_EXCEPTIONT("Command line levels must match FIXED_RELEMENTS");
	}
#else
	const int nRElements = m_grid.GetRElements();
#endif

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
		for (int c = 0; c < nodevars.size(); c++) {
		for (int i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (int j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {
		for (int k = 0; k < nRElements; k++) {
			dChecksums[nodevars[c]] +=
				  (*pDataNode)(nodevars[c],i,j,k)
				* m_dataElementAreaNode(i,j,k);
		}
		}
		}
		}

		for (int c = 0; c < redgevars.size(); c++) {
		for (int i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (int j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {
		for (int k = 0; k <= nRElements; k++) {
			dChecksums[redgevars[c]] +=
				  (*pDataREdge)(redgevars[c],i,j,k)
				* m_dataElementAreaREdge(i,j,k);
		}
		}
		}
		}

	// ChecksumType_L1
	} else if (eChecksumType == ChecksumType_L1) {
		for (int c = 0; c < nodevars.size(); c++) {
		for (int i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (int j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {
		for (int k = 0; k < nRElements; k++) {
			const double dValue = fabs((*pDataNode)(nodevars[c],i,j,k));
			dChecksums[nodevars[c]] +=
				dValue * m_dataElementAreaNode(i,j,k);
		}
		}
		}
		}

		for (int c = 0; c < redgevars.size(); c++) {
		for (int i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (int j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {
		for (int k = 0; k <= nRElements; k++) {
			const double dValue = fabs((*pDataREdge)(redgevars[c],i,j,k));
			dChecksums[redgevars[c]] +=
				dValue * m_dataElementAreaREdge(i,j,k);
		}
		}
		}
		}

	// ChecksumType_L2
	} else if (eChecksumType == ChecksumType_L2) {
		for (int c = 0; c < nodevars.size(); c++) {
		for (int i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (int j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {
		for (int k = 0; k < nRElements; k++) {
			const double dValue = (*pDataNode)(nodevars[c],i,j,k);
			dChecksums[nodevars[c]] +=
				dValue * dValue * m_dataElementAreaNode(i,j,k);
		}
		}
		}
		}

		for (int c = 0; c < redgevars.size(); c++) {
		for (int i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (int j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {
		for (int k = 0; k <= nRElements; k++) {
			const double dValue = (*pDataREdge)(redgevars[c],i,j,k);
			dChecksums[redgevars[c]] +=
				dValue * dValue * m_dataElementAreaREdge(i,j,k);
		}
		}
		}
		}

	// ChecksumType_Linf
	} else if (eChecksumType == ChecksumType_Linf) {
		for (int c = 0; c < nodevars.size(); c++) {
		for (int i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (int j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {
		for (int k = 0; k < nRElements; k++) {
			const double dValue = (*pDataNode)(nodevars[c],i,j,k);
			if (fabs(dValue) > dChecksums[nodevars[c]]) {
				dChecksums[nodevars[c]] = fabs(dValue);
			}
		}
		}
		}
		}

		for (int c = 0; c < redgevars.size(); c++) {
		for (int i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (int j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {
		for (int k = 0; k <= nRElements; k++) {
			const double dValue = (*pDataREdge)(redgevars[c],i,j,k);
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

#if defined(FIXED_RELEMENTS)
	const int nRElements = FIXED_RELEMENTS;
	if (nRElements != m_grid.GetRElements()) {
		_EXCEPTIONT("Command line levels must match FIXED_RELEMENTS");
	}
#else
	const int nRElements = m_grid.GetRElements();
#endif

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
		for (int i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (int j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {
		for (int k = 0; k < nRElements; k++) {
			double dUdotU =
				+ m_dataContraMetric2DB(i,j,1)
					* dataNode(UIx,i,j,k)
					* dataNode(UIx,i,j,k)
				- 2.0 * m_dataContraMetric2DA(i,j,1)
					* dataNode(UIx,i,j,k)
					* dataNode(VIx,i,j,k)
				+ m_dataContraMetric2DA(i,j,0)
					* dataNode(VIx,i,j,k)
					* dataNode(VIx,i,j,k);

			dUdotU *= m_dataJacobian2D(i,j)
				* m_dataJacobian2D(i,j);

			const double dKineticEnergy =
				0.5 * dUdotU * (dataNode(HIx,i,j,k)
						- m_dataTopography(i,j));

			const double dPotentialEnergy =
				0.5 * phys.GetG()
					* (dataNode(HIx,i,j,k)
						* dataNode(HIx,i,j,k)
					- m_dataTopography(i,j)
						* m_dataTopography(i,j));

			dLocalEnergy += m_dataElementAreaNode(i,j,k)
				* (dKineticEnergy + dPotentialEnergy);
		}
		}
		}

	// Non-hydrostatic model with vertical velocity on nodes
	} else if (m_grid.GetVarLocation(WIx) == DataLocation_Node) {

		// Loop over all elements
		for (int i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (int j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {
		for (int k = 0; k < nRElements; k++) {

			const double dCovUa = dataNode(UIx,i,j,k);
			const double dCovUb = dataNode(VIx,i,j,k);
			const double dCovUx = dataNode(WIx,i,j,k);

			const double dConUa =
				  m_dataContraMetricA(i,j,k,0) * dCovUa
				+ m_dataContraMetricA(i,j,k,1) * dCovUb
				+ m_dataContraMetricA(i,j,k,2) * dCovUx;

			const double dConUb =
				  m_dataContraMetricB(i,j,k,0) * dCovUa
				+ m_dataContraMetricB(i,j,k,1) * dCovUb
				+ m_dataContraMetricB(i,j,k,2) * dCovUx;

			const double dConUx =
				  m_dataContraMetricXi(i,j,k,0) * dCovUa
				+ m_dataContraMetricXi(i,j,k,1) * dCovUb
				+ m_dataContraMetricXi(i,j,k,2) * dCovUx;

			const double dUdotU = dConUa * dCovUa + dConUb * dCovUb + dConUx * dCovUx;

			const double dKineticEnergy =
				0.5 * dataNode(RIx,i,j,k) * dUdotU;

#ifdef FORMULATION_PRESSURE
			const double dPressure = dataNode(PIx,i,j,k);
#endif
#if defined(FORMULATION_RHOTHETA_PI) || defined(FORMULATION_RHOTHETA_P)
			const double dPressure = phys.PressureFromRhoTheta(dataNode(PIx,i,j,k));
#endif
#if defined(FORMULATION_THETA) || defined(FORMULATION_THETA_FLUX)
			const double dPressure =
				phys.PressureFromRhoTheta(
					dataNode(RIx,i,j,k)
					* dataNode(PIx,i,j,k));
#endif

			const double dInternalEnergy =
				dPressure / (phys.GetGamma() - 1.0);

			const double dPotentialEnergy =
				phys.GetG()
				* dataNode(RIx,i,j,k)
				* m_dataZLevels(i,j,k);

			dLocalEnergy += m_dataElementAreaNode(i,j,k)
				* (dKineticEnergy + dInternalEnergy + dPotentialEnergy);
		}
		}
		}

	// Non-hydrostatic model with vertical velocity on interfaces
	} else {

		// Data on interfaces
		const DataArray4D<double> & dataREdge =
			m_datavecStateREdge[iDataIndex];

		// Loop over all elements
		for (int i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (int j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {
		for (int k = 0; k < nRElements; k++) {

			const double dCovUa = dataNode(UIx,i,j,k);
			const double dCovUb = dataNode(VIx,i,j,k);
			const double dCovUx = dataNode(WIx,i,j,k);

			const double dConUa =
				  m_dataContraMetricA(i,j,k,0) * dCovUa
				+ m_dataContraMetricA(i,j,k,1) * dCovUb
				+ m_dataContraMetricA(i,j,k,2) * dCovUx;

			const double dConUb =
				  m_dataContraMetricB(i,j,k,0) * dCovUa
				+ m_dataContraMetricB(i,j,k,1) * dCovUb
				+ m_dataContraMetricB(i,j,k,2) * dCovUx;

			double dUdotU = dConUa * dCovUa + dConUb * dCovUb;

			dUdotU += m_dataContraMetricXi(i,j,k,0) * dCovUa * dCovUx;
			dUdotU += m_dataContraMetricXi(i,j,k,1) * dCovUb * dCovUx;

			const double dKineticEnergy =
				0.5 * dataNode(RIx,i,j,k) * dUdotU;

#ifdef FORMULATION_PRESSURE
			const double dPressure = dataNode(PIx,i,j,k);
#endif
#if defined(FORMULATION_RHOTHETA_PI) || defined(FORMULATION_RHOTHETA_P)
			const double dPressure = phys.PressureFromRhoTheta(dataNode(PIx,i,j,k));
#endif
#if defined(FORMULATION_THETA) || defined(FORMULATION_THETA_FLUX)
			const double dPressure =
				phys.PressureFromRhoTheta(
					dataNode(RIx,i,j,k)
					* dataNode(PIx,i,j,k));
#endif

			const double dInternalEnergy =
				dPressure / (phys.GetGamma() - 1.0);

			const double dPotentialEnergy =
				phys.GetG() * dataNode(RIx,i,j,k)
				* m_dataZLevels(i,j,k);

			dLocalEnergy += m_dataElementAreaNode(i,j,k)
				* (dKineticEnergy + dInternalEnergy + dPotentialEnergy);
		}
		}
		}

		for (int i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (int j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {
		for (int k = 0; k <= nRElements; k++) {

			const double dCovUx = dataREdge(WIx,i,j,k);

			const double dKineticEnergy =
				0.5 * dataREdge(RIx,i,j,k)
				* m_dataContraMetricXiREdge(i,j,k,2)
				* dCovUx * dCovUx;

			dLocalEnergy +=
				m_dataElementAreaREdge(i,j,k)
				* dKineticEnergy;
		}
		}
		}
	}

	return dLocalEnergy;
}

///////////////////////////////////////////////////////////////////////////////

double GridPatch::ComputeTotalPotentialEnstrophy(
	int iDataIndex
) {
	// Physical constants
	const PhysicalConstants & phys = m_grid.GetModel().GetPhysicalConstants();

#if defined(FIXED_RELEMENTS)
	const int nRElements = FIXED_RELEMENTS;
	if (nRElements != m_grid.GetRElements()) {
		_EXCEPTIONT("Command line levels must match FIXED_RELEMENTS");
	}
#else
	const int nRElements = m_grid.GetRElements();
#endif

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
		for (int i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (int j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {
		for (int k = 0; k < nRElements; k++) {
			const double dPlanetaryVorticity =
				2.0 * phys.GetOmega() * sin(m_dataLat(i,j));

			const double dAbsoluteVorticity =
				m_dataVorticity(i,j,k) + dPlanetaryVorticity;

			dLocalPotentialEnstrophy +=
				m_dataElementAreaNode(i,j,k)
					* 0.5 * dAbsoluteVorticity * dAbsoluteVorticity
					/ (dataNode(HIx,i,j,k)
						- m_dataTopography(i,j));
		}
		}
		}

	} else {

		// Set to zero
		dLocalPotentialEnstrophy = 0.0;

		// Loop over all elements
		for (int i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (int j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {
		for (int k = 0; k < nRElements; k++) {

			// Zonal momentum
			dLocalPotentialEnstrophy +=
				m_dataElementAreaNode(i,j,k)
				* dataNode(RIx,i,j,k)
				* dataNode(UIx,i,j,k);
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

#if defined(FIXED_RELEMENTS)
	const int nRElements = FIXED_RELEMENTS;
	if (nRElements != m_grid.GetRElements()) {
		_EXCEPTIONT("Command line levels must match FIXED_RELEMENTS");
	}
#else
	const int nRElements = m_grid.GetRElements();
#endif

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
		for (int i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (int j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {
		for (int k = 0; k < nRElements; k++) {

			// Vertical momentum
			dLocalVerticalMomentum +=
				m_dataElementAreaNode(i,j,k)
				* dataNode(RIx,i,j,k)
				* dataNode(WIx,i,j,k);
		}
		}
		}

	}

	return dLocalVerticalMomentum;
}

///////////////////////////////////////////////////////////////////////////////

void GridPatch::PackExchangeBuffer(
	DataType eDataType,
	int iDataIndex,
	ExchangeBuffer & exbuf
) {
	// Check exchange buffer target
	if (exbuf.m_ixSourcePatch != m_ixPatch) {
		_EXCEPTIONT("ExchangeBuffer patch index mismatch");
	}

	// State data
	if (eDataType == DataType_State) {
		if ((iDataIndex < 0) || (iDataIndex > m_datavecStateNode.size())) {
			_EXCEPTIONT("Invalid state data instance.");
		}

		exbuf.Pack(m_grid, m_datavecStateNode[iDataIndex]);
		exbuf.Pack(m_grid, m_datavecStateREdge[iDataIndex]);

	// Tracer data
	} else if (eDataType == DataType_Tracers) {
		if ((iDataIndex < 0) || (iDataIndex > m_datavecTracers.size())) {
			_EXCEPTIONT("Invalid tracers data instance.");
		}

		exbuf.Pack(m_grid, m_datavecTracers[iDataIndex]);
	
	// Vorticity data
	} else if (eDataType == DataType_Vorticity) {
		exbuf.Pack(m_dataVorticity);

	// Divergence data
	} else if (eDataType == DataType_Divergence) {
		exbuf.Pack(m_dataDivergence);

	// Temperature data
	} else if (eDataType == DataType_Temperature) {
		exbuf.Pack(m_dataTemperature);

	// Zonal drag data
	} else if (eDataType == DataType_ZonalForce) {
		exbuf.Pack(m_dataZonalForce);

	// Richardson data
	} else if (eDataType == DataType_Richardson) {
		exbuf.Pack(m_dataRichardson);

	// Convective stability data
	} else if (eDataType == DataType_Convective) {
		exbuf.Pack(m_dataConvective);

	// Topography derivative data
	} else if (eDataType == DataType_TopographyDeriv) {
		exbuf.Pack(m_dataTopographyDeriv);

	// Invalid data
	} else {
		_EXCEPTIONT("Invalid DataType");
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatch::UnpackExchangeBuffer(
	DataType eDataType,
	int iDataIndex,
	ExchangeBuffer & exbuf
) {
	// Check exchange buffer target
	if (exbuf.m_ixSourcePatch != m_ixPatch) {
		_EXCEPTIONT("ExchangeBuffer patch index mismatch");
	}

	// State data
	if (eDataType == DataType_State) {
		if ((iDataIndex < 0) || (iDataIndex > m_datavecStateNode.size())) {
			_EXCEPTIONT("Invalid state data instance.");
		}

		exbuf.Unpack(m_grid, m_datavecStateNode[iDataIndex]);
		exbuf.Unpack(m_grid, m_datavecStateREdge[iDataIndex]);

	// Tracer data
	} else if (eDataType == DataType_Tracers) {
		if ((iDataIndex < 0) || (iDataIndex > m_datavecTracers.size())) {
			_EXCEPTIONT("Invalid tracers data instance.");
		}

		exbuf.Unpack(m_grid, m_datavecTracers[iDataIndex]);

	// Vorticity data
	} else if (eDataType == DataType_Vorticity) {
		exbuf.Unpack(m_dataVorticity);

	// Divergence data
	} else if (eDataType == DataType_Divergence) {
		exbuf.Unpack(m_dataDivergence);

	// Temperature data
	} else if (eDataType == DataType_Temperature) {
		exbuf.Unpack(m_dataTemperature);

	// Zonal drag data
        } else if (eDataType == DataType_ZonalForce) {
                exbuf.Unpack(m_dataZonalForce);

	// Richardson data
	} else if (eDataType == DataType_Richardson) {
		exbuf.Unpack(m_dataRichardson);

	// Convective stability data
	} else if (eDataType == DataType_Convective) {
		exbuf.Unpack(m_dataConvective);

	// Topographic derivatives
	} else if (eDataType == DataType_TopographyDeriv) {
		exbuf.Unpack(m_dataTopographyDeriv);

	// Invalid data
	} else {
		_EXCEPTIONT("Invalid DataType");
	}
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
			_EXCEPTIONT("Invalid ixData index in ZeroData.");
		}

		m_datavecStateNode [ixData].Zero();
		m_datavecStateREdge[ixData].Zero();

	// Check bounds on ixDest for Tracers data
	} else if (eDataType == DataType_Tracers) {
		if ((ixData < 0) || (ixData >= m_datavecTracers.size())) {
			_EXCEPTIONT("Invalid ixData index in ZeroData.");
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
