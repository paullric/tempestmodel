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

#include "mpi.h"

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

int GridPatch::GetTotalNodeCount2D() const {
	return (m_box.GetATotalWidth() * m_box.GetBTotalWidth());
}

///////////////////////////////////////////////////////////////////////////////

int GridPatch::GetTotalNodeCount(
	DataLocation loc
) const {
	if (loc == DataLocation_Node) {
		return (m_box.GetTotalNodes() * m_grid.GetRElements());

	} else if (loc == DataLocation_REdge) {
		return (m_box.GetTotalNodes() * (m_grid.GetRElements()+1));

	} else {
		_EXCEPTIONT("Invalid DataLocation");
	}
}

///////////////////////////////////////////////////////////////////////////////

int GridPatch::GetTotalDegreesOfFreedom(
	DataType eDataType,
	DataLocation eDataLocation
) const {
	
	// Take into account staggering of State data
	if ((eDataType == DataType_State) ||
		(eDataType == DataType_RefState)
	) {
		int nComponents = m_grid.GetModel().GetEquationSet().GetComponents();

		if (eDataLocation == DataLocation_None) {
			return (GetTotalNodeCount2D()
				* m_grid.GetDegreesOfFreedomPerColumn());

		} else if (eDataLocation == DataLocation_Node) {
			return (GetTotalNodeCount2D()
				* m_grid.GetRElements()
				* nComponents);

		} else if (eDataLocation == DataLocation_REdge) {
			return (GetTotalNodeCount2D()
				* (m_grid.GetRElements()+1)
				* nComponents);

		} else {
			_EXCEPTIONT("Invalid DataLocation");
		}

	// All tracers on model levels
	} else if (eDataType == DataType_Tracers) {
		int nTracers = m_grid.GetModel().GetEquationSet().GetTracers();

		return (GetTotalNodeCount2D()
			* m_grid.GetRElements()
			* nTracers);

	// Topography only at surface
	} else if (eDataType == DataType_Topography) {
		return GetTotalNodeCount2D();

	// Rayleigh strength
	} else if (eDataType == DataType_RayleighStrength) {
		if (eDataLocation == DataLocation_Node) {
			return GetTotalNodeCount2D() * m_grid.GetRElements();

		} else if (eDataLocation == DataLocation_REdge) {
			return GetTotalNodeCount2D() * (m_grid.GetRElements() + 1);

		} else {
			_EXCEPTIONT("Invalid DataLocation");
		}

	// Invalid DataType
	} else {
		_EXCEPTIONT("(UNIMPLEMENTED) Invalid DataType");
	}
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

void GridPatch::InitializeDataLocal() {
	if (m_fContainsData) {
		_EXCEPTIONT(
			"Attempting to initialize a previously initialized GridPatch.");
	}

	// This patch contains data
	m_fContainsData = true;

	// Set the processor
	MPI_Comm_rank(MPI_COMM_WORLD, &m_iProcessor);

	// Jacobian at each node (2D)
	m_dataJacobian2D.Initialize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth());

	// Contravariant metric (2D) components at each node
	m_dataContraMetric2DA.Initialize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		2);

	m_dataContraMetric2DB.Initialize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		2);

	// Covariant metric (2D) components at each node
	m_dataCovMetric2DA.Initialize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		2);

	m_dataCovMetric2DB.Initialize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		2);

	// Jacobian at each node
	m_dataJacobian.Initialize(
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth());

	// Jacobian at each interface
	m_dataJacobianREdge.Initialize(
		m_grid.GetRElements()+1,
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth());

	// Contravariant metric components at each node
	m_dataContraMetricA.Initialize(
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		3);

	m_dataContraMetricB.Initialize(
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		3);

	m_dataContraMetricXi.Initialize(
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		3);

	// Covariant metric components at each node
	m_dataCovMetricA.Initialize(
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		3);

	m_dataCovMetricB.Initialize(
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		3);

	m_dataCovMetricXi.Initialize(
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		3);

	// Xi contravariant metric on interfaces
	m_dataContraMetricAREdge.Initialize(
		m_grid.GetRElements()+1,
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		3);

	m_dataContraMetricBREdge.Initialize(
		m_grid.GetRElements()+1,
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		3);

	m_dataContraMetricXiREdge.Initialize(
		m_grid.GetRElements()+1,
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		3);

	// Vertical coordinate transform (derivatives of the radius)
	m_dataDerivRNode.Initialize(
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		3);

	m_dataDerivRREdge.Initialize(
		m_grid.GetRElements()+1,
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		3);

	// Element area at each node
	m_dataElementArea.Initialize(
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth());

	// Element area at each interface
	m_dataElementAreaREdge.Initialize(
		m_grid.GetRElements()+1,
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth());

	// Topography height at each node
	m_dataTopography.Initialize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth());

	// Topography derivatives at each node
	m_dataTopographyDeriv.Initialize(
		DataType_TopographyDeriv,
		DataLocation_Node,
		2,
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_box.GetHaloElements());

	// Longitude at each node
	m_dataLon.Initialize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth());

	// Latitude at each node
	m_dataLat.Initialize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth());

	// Coriolis parameter at each node
	m_dataCoriolisF.Initialize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth());

	// Radial coordinate at each level
	m_dataZLevels.Initialize(
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth());

	// Radial coordinate at each interface
	m_dataZInterfaces.Initialize(
		m_grid.GetRElements()+1,
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth());

	// Get the model
	const Model & model = m_grid.GetModel();

	// Get the equation set
	const EquationSet & eqn = model.GetEquationSet();

	// Initialize reference state
	m_dataRefStateNode.Initialize(
		DataType_State,
		DataLocation_Node,
		eqn.GetComponents(),
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_box.GetHaloElements());

	m_dataRefStateREdge.Initialize(
		DataType_State,
		DataLocation_REdge,
		eqn.GetComponents(),
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_box.GetHaloElements());

	// Initialize component data
	m_datavecStateNode .resize(model.GetComponentDataInstances());
	m_datavecStateREdge.resize(model.GetComponentDataInstances());

	for (int m = 0; m < model.GetComponentDataInstances(); m++) {
		m_datavecStateNode[m].Initialize(
			DataType_State,
			DataLocation_Node,
			eqn.GetComponents(),
			m_grid.GetRElements(),
			m_box.GetATotalWidth(),
			m_box.GetBTotalWidth(),
			m_box.GetHaloElements());

		m_datavecStateREdge[m].Initialize(
			DataType_State,
			DataLocation_REdge,
			eqn.GetComponents(),
			m_grid.GetRElements(),
			m_box.GetATotalWidth(),
			m_box.GetBTotalWidth(),
			m_box.GetHaloElements());
	}

	// Initialize tracer data
	m_datavecTracers.resize(model.GetTracerDataInstances());

	for (int m = 0; m < model.GetTracerDataInstances(); m++) {
		m_datavecTracers[m].Initialize(
			DataType_Tracers,
			DataLocation_Node,
			eqn.GetTracers(),
			m_grid.GetRElements(),
			m_box.GetATotalWidth(),
			m_box.GetBTotalWidth(),
			m_box.GetHaloElements());
	}

#pragma message "Make these processes more generic"
	// Initialize auxiliary data
	m_datavecAuxNode.resize(2);
	m_datavecAuxREdge.resize(2);

	m_datavecAuxNode[0].resize(model.GetHorizontalDynamicsAuxDataCount());
	m_datavecAuxREdge[0].resize(model.GetHorizontalDynamicsAuxDataCount());

	m_datavecAuxNode[1].resize(model.GetVerticalDynamicsAuxDataCount());
	m_datavecAuxREdge[1].resize(model.GetVerticalDynamicsAuxDataCount());

	for (int m = 0; m < model.GetHorizontalDynamicsAuxDataCount(); m++) {
		m_datavecAuxNode[0][m].Initialize(
			DataType_Auxiliary,
			DataLocation_Node,
			m_grid.GetRElements(),
			m_box.GetATotalWidth(),
			m_box.GetBTotalWidth(),
			m_box.GetHaloElements());

		m_datavecAuxREdge[0][m].Initialize(
			DataType_Auxiliary,
			DataLocation_REdge,
			m_grid.GetRElements(),
			m_box.GetATotalWidth(),
			m_box.GetBTotalWidth(),
			m_box.GetHaloElements());
	}

	for (int m = 0; m < model.GetVerticalDynamicsAuxDataCount(); m++) {
		m_datavecAuxNode[1][m].Initialize(
			DataType_Auxiliary,
			DataLocation_Node,
			m_grid.GetRElements(),
			m_box.GetATotalWidth(),
			m_box.GetBTotalWidth(),
			m_box.GetHaloElements());

		m_datavecAuxREdge[1][m].Initialize(
			DataType_Auxiliary,
			DataLocation_REdge,
			m_grid.GetRElements(),
			m_box.GetATotalWidth(),
			m_box.GetBTotalWidth(),
			m_box.GetHaloElements());
	}

	// Pressure data
	m_dataPressure.Initialize(
		DataType_Pressure,
		DataLocation_Node,
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_box.GetHaloElements());
/*
	m_dataDaPressure.Initialize(
		DataType_Pressure,
		DataLocation_Node,
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_box.GetHaloElements());

	m_dataDbPressure.Initialize(
		DataType_Pressure,
		DataLocation_Node,
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_box.GetHaloElements());
*/
	m_dataDxPressure.Initialize(
		DataType_Pressure,
		DataLocation_Node,
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_box.GetHaloElements());

	// Vorticity data
	m_dataVorticity.Initialize(
		DataType_Vorticity,
		DataLocation_Node,
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_box.GetHaloElements());

	// Divergence data
	m_dataDivergence.Initialize(
		DataType_Divergence,
		DataLocation_Node,
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_box.GetHaloElements());

	// Temperature data
	m_dataTemperature.Initialize(
		DataType_Temperature,
		DataLocation_Node,
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_box.GetHaloElements());

	// Rayleigh friction strength
	m_dataRayleighStrengthNode.Initialize(
		DataType_None,
		DataLocation_Node,
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_box.GetHaloElements());

	// Rayleigh friction strength
	m_dataRayleighStrengthREdge.Initialize(
		DataType_None,
		DataLocation_REdge,
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		m_box.GetHaloElements());

}

///////////////////////////////////////////////////////////////////////////////

void GridPatch::DeinitializeData() {
	if (!m_fContainsData) {
		_EXCEPTIONT("Attempting to deinitialize a stub GridPatch.");
	}

	m_fContainsData = false;

	m_dataJacobian2D.Deinitialize();
	m_dataContraMetric2DA.Deinitialize();
	m_dataContraMetric2DB.Deinitialize();
	m_dataCovMetric2DA.Deinitialize();
	m_dataCovMetric2DB.Deinitialize();

	m_dataJacobian.Deinitialize();
	m_dataJacobianREdge.Deinitialize();
	m_dataContraMetricA.Deinitialize();
	m_dataContraMetricB.Deinitialize();
	m_dataContraMetricXi.Deinitialize();
	m_dataCovMetricA.Deinitialize();
	m_dataCovMetricB.Deinitialize();
	m_dataCovMetricXi.Deinitialize();
	m_dataContraMetricXiREdge.Deinitialize();
	m_dataDerivRNode.Deinitialize();
	m_dataElementArea.Deinitialize();
	m_dataElementAreaREdge.Deinitialize();
	m_dataTopography.Deinitialize();
	m_dataTopographyDeriv.Deinitialize();
	m_dataLon.Deinitialize();
	m_dataLat.Deinitialize();
	m_dataZLevels.Deinitialize();
	m_dataZInterfaces.Deinitialize();
	m_datavecStateNode.Deinitialize();
	m_datavecStateREdge.Deinitialize();
	m_datavecTracers.Deinitialize();

	for (int n = 0; n < m_datavecAuxNode.size(); n++) {
	for (int m = 0; m < m_datavecAuxNode[n].size(); m++) {
		m_datavecAuxNode[n][m].Deinitialize();
	}
	}

	for (int n = 0; n < m_datavecAuxREdge.size(); n++) {
	for (int m = 0; m < m_datavecAuxREdge[n].size(); m++) {
		m_datavecAuxREdge[n][m].Deinitialize();
	}
	}

	m_dataPressure.Deinitialize();
	m_dataDxPressure.Deinitialize();

	m_dataVorticity.Deinitialize();
	m_dataDivergence.Deinitialize();
	m_dataTemperature.Deinitialize();
	m_dataRayleighStrengthNode.Deinitialize();
	m_dataRayleighStrengthREdge.Deinitialize();
}

///////////////////////////////////////////////////////////////////////////////

void GridPatch::ExteriorConnect(
	Direction dirFirst,
	const GridPatch * pPatchSecond
) {
	// First patch coordinate index
	int ixFirst;
	int ixSecond;

	if ((dirFirst == Direction_Right) ||
		(dirFirst == Direction_Left)
	) {
		ixFirst  = m_box.GetBInteriorBegin();
		ixSecond = m_box.GetBInteriorEnd();

	} else if (
		(dirFirst == Direction_Top) ||
		(dirFirst == Direction_Bottom)
	) {
		ixFirst  = m_box.GetAInteriorBegin();
		ixSecond = m_box.GetAInteriorEnd();

	} else if (dirFirst == Direction_TopRight) {
		ixFirst  = m_box.GetAInteriorEnd()-1;
		ixSecond = m_box.GetBInteriorEnd()-1;

	} else if (dirFirst == Direction_TopLeft) {
		ixFirst  = m_box.GetAInteriorBegin();
		ixSecond = m_box.GetBInteriorEnd()-1;

	} else if (dirFirst == Direction_BottomLeft) {
		ixFirst  = m_box.GetAInteriorBegin();
		ixSecond = m_box.GetBInteriorBegin();

	} else if (dirFirst == Direction_BottomRight) {
		ixFirst  = m_box.GetAInteriorEnd()-1;
		ixSecond = m_box.GetBInteriorBegin();

	} else {
		_EXCEPTIONT("Invalid direction");
	}

	// Exterior connect
	ExteriorConnect(
		dirFirst,
		pPatchSecond,
		ixFirst,
		ixSecond);
}

///////////////////////////////////////////////////////////////////////////////

void GridPatch::ExteriorConnect(
	Direction dirFirst,
	const GridPatch * pPatchSecond,
	int ixFirst,
	int ixSecond
) {
	// Check for NULL patches (do nothing)
	if (pPatchSecond == NULL) {
		return;
	}

	// Get number of components
	const Model & model = m_grid.GetModel();

	const EquationSet & eqn = model.GetEquationSet();

	int nStateTracerMaxVariables;
	if (eqn.GetComponents() > eqn.GetTracers()) {
		nStateTracerMaxVariables = eqn.GetComponents();
	} else {
		nStateTracerMaxVariables = eqn.GetTracers();
	}

	int nRElements = m_grid.GetRElements();

	int nHaloElements = model.GetHaloElements();

	// Get the opposing direction
	Direction dirOpposing = Direction_Unreachable;
	bool fReverseDirection = false;
	bool fFlippedCoordinate = false;

	m_grid.GetOpposingDirection(
		m_box.GetPanel(),
		pPatchSecond->GetPatchBox().GetPanel(),
		dirFirst,
		dirOpposing,
		fReverseDirection,
		fFlippedCoordinate);

	// Determine the size of the boundary (number of elements along exterior
	// edge).  Used in computing the size of the send/recv buffers.
	int nBoundarySize;
	if ((dirFirst == Direction_Right) ||
		(dirFirst == Direction_Top) ||
		(dirFirst == Direction_Left) ||
		(dirFirst == Direction_Bottom)
	) {
		nBoundarySize = ixSecond - ixFirst;
	} else {
		nBoundarySize = nHaloElements;
	}

	if ((dirFirst == Direction_TopRight) && (
		(ixFirst < m_box.GetAInteriorBegin() + nBoundarySize - 1) ||
		(ixSecond < m_box.GetBInteriorBegin() + nBoundarySize - 1)
	)) {
		_EXCEPTIONT("Insufficient interior elements to build "
			"diagonal connection.");
	}

	if ((dirFirst == Direction_TopLeft) && (
		(ixFirst > m_box.GetAInteriorEnd() - nBoundarySize) ||
		(ixSecond < m_box.GetBInteriorBegin() + nBoundarySize - 1)
	)) {
		_EXCEPTIONT("Insufficient interior elements to build "
			"diagonal connection.");
	}

	if ((dirFirst == Direction_BottomLeft) && (
		(ixFirst > m_box.GetAInteriorEnd() - nBoundarySize) ||
		(ixSecond > m_box.GetBInteriorEnd() - nBoundarySize)
	)) {
		_EXCEPTIONT("Insufficient interior elements to build "
			"diagonal connection.");
	}

	if ((dirFirst == Direction_BottomRight) && (
		(ixFirst < m_box.GetAInteriorBegin() + nBoundarySize - 1) ||
		(ixSecond > m_box.GetBInteriorEnd() - nBoundarySize)
	)) {
		_EXCEPTIONT("Insufficient interior elements to build "
			"diagonal connection.");
	}

	// Add an external neighbor to the patch
	ExteriorNeighbor * pNeighbor =
		new ExteriorNeighbor(
			&m_connect,
			dirFirst,
			dirOpposing,
			pPatchSecond->GetPatchIndex(),
			fReverseDirection,
			fFlippedCoordinate,
			nBoundarySize,
			ixFirst,
			ixSecond);

	pNeighbor->InitializeBuffers(
		nRElements,
		nHaloElements,
		nStateTracerMaxVariables);

	m_connect.AddExteriorNeighbor(pNeighbor);
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

		const GridData4D & dataState = m_datavecStateNode[iDataIndex];

		for (k = 0; k < m_grid.GetRElements(); k++) {
		for (i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {
			m_dataTemperature[k][i][j] =
				dataState[PIx][k][i][j]
				/ (dataState[RIx][k][i][j] * phys.GetR());
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

		const GridData4D & dataState = m_datavecStateREdge[iDataIndex];

		for (k = 0; k <= m_grid.GetRElements(); k++) {
		for (i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {
			m_dataTemperature[k][i][j] =
				dataState[PIx][k][i][j]
				/ (dataState[RIx][k][i][j] * phys.GetR());
		}
		}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void GridPatch::Checksum(
	DataType eDataType,
	DataVector<double> & dChecksums,
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
	GridData4D const * pDataNode;
	GridData4D const * pDataREdge;

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
	const GridData4D & dataNode = m_datavecStateNode[iDataIndex];

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

#ifdef USE_COVARIANT_VELOCITIES
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
#else
			double dUdotU =
				+ m_dataCovMetric2DA[i][j][0]
					* dataNode[UIx][k][i][j] * dataNode[UIx][k][i][j]
				+ (m_dataCovMetric2DA[i][j][1] + m_dataCovMetric2DB[i][j][0])
					* dataNode[UIx][k][i][j] * dataNode[VIx][k][i][j]
				+ m_dataCovMetric2DB[i][j][1]
					* dataNode[VIx][k][i][j] * dataNode[VIx][k][i][j];

			dUdotU += dataNode[WIx][k][i][j] * dataNode[WIx][k][i][j];
#endif

			double dKineticEnergy =
				0.5 * dataNode[RIx][k][i][j] * dUdotU;

#ifdef FORMULATION_PRESSURE
			double dPressure = dataNode[PIx][k][i][j];
#endif
#if defined(FORMULATION_RHOTHETA_PI) || defined(FORMULATION_RHOTHETA_P)
			double dPressure = phys.PressureFromRhoTheta(dataNode[PIx][k][i][j]);
#endif
#ifdef FORMULATION_THETA
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
	const GridData4D & dataNode = m_datavecStateNode[iDataIndex];

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

#pragma message "Total potential enstrophy"
		// Unimplemented
		dLocalPotentialEnstrophy = 0.0;
	}

	return dLocalPotentialEnstrophy;
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
			pNeighbor->Unpack(m_datavecStateNode[iDataIndex]);
			pNeighbor->Unpack(m_datavecStateREdge[iDataIndex]);
		}

	// Tracer data
	} else if (eDataType == DataType_Tracers) {
		if ((iDataIndex < 0) || (iDataIndex > m_datavecTracers.size())) {
			_EXCEPTIONT("Invalid tracers data instance.");
		}

		Neighbor * pNeighbor;
		while ((pNeighbor = m_connect.WaitReceive()) != NULL) {
			pNeighbor->Unpack(m_datavecTracers[iDataIndex]);
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
	const DataVector<double> & dCoeff,
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
	const DataVector<double> & dAlpha,
	const DataVector<double> & dBeta,
	const DataVector<int> & iPanel,
	DataType eDataType,
	DataLocation eDataLocation,
	bool fInterpAllVariables,
	DataMatrix3D<double> & dInterpData,
	bool fIncludeReferenceState,
	bool fConvertToPrimitive
) {
	_EXCEPTIONT("Unimplemented.");
}

///////////////////////////////////////////////////////////////////////////////

