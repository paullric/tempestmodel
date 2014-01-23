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

	// Jacobian at each node
	m_dataJacobian.Initialize(
		m_grid.GetRElements(),
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

	// Christoffel symbol components at each node
	m_dataChristoffelA.Initialize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		3);

	m_dataChristoffelB.Initialize(
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		3);

	// Vertical metric components on nodes or interfaces
#pragma message "This should pull from EquationSet by name"
	int nWLevels = m_grid.GetRElements();
	if (m_grid.GetModel().GetEquationSet().GetDimensionality() == 3) {
		if (m_grid.GetVarLocation(3) == DataLocation_REdge) {
			nWLevels = m_grid.GetRElements()+1;
		}
	}

	m_dataContraMetricXi.Initialize(
		nWLevels,
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		3);

	m_dataChristoffelXi.Initialize(
		nWLevels,
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		6);

	// Orthonormalization components
	m_dataOrthonormNode.Initialize(
		m_grid.GetRElements(),
		m_box.GetATotalWidth(),
		m_box.GetBTotalWidth(),
		3);

	m_dataOrthonormREdge.Initialize(
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

	m_dataJacobian.Deinitialize();
	m_dataContraMetricA.Deinitialize();
	m_dataContraMetricB.Deinitialize();
	m_dataContraMetricXi.Deinitialize();
	m_dataChristoffelA.Deinitialize();
	m_dataChristoffelB.Deinitialize();
	m_dataChristoffelXi.Deinitialize();
	m_dataElementArea.Deinitialize();
	m_dataElementAreaREdge.Deinitialize();
	m_dataTopography.Deinitialize();
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

	m_dataVorticity.Deinitialize();
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

	// Implement
	if (eDataType == DataType_Tracers) {
		_EXCEPTIONT("Not implemented.");
	}

	const GridData4D & dataNode  = m_datavecStateNode[iDataIndex];
	const GridData4D & dataREdge = m_datavecStateREdge[iDataIndex];

	// Variables on nodes
#pragma "Construct a Grid vector of vectors indicating all variables at each location"
	int nComponents = m_grid.GetModel().GetEquationSet().GetComponents();

	std::vector<int> nodevars;
	std::vector<int> redgevars;
	for (int c = 0; c < nComponents; c++) {
		if (m_grid.GetVarLocation(c) == DataLocation_Node) {
			nodevars.push_back(c);
		} else if (m_grid.GetVarLocation(c) == DataLocation_REdge) {
			redgevars.push_back(c);
		} else {
			_EXCEPTIONT("Not implemented.");
		}
	}

	// ChecksumType_Sum
	if (eChecksumType == ChecksumType_Sum) {
		for (c = 0; c < nodevars.size(); c++) {
		for (k = 0; k < m_grid.GetRElements(); k++) {
		for (i = m_box.GetAInteriorBegin(); i < m_box.GetAInteriorEnd(); i++) {
		for (j = m_box.GetBInteriorBegin(); j < m_box.GetBInteriorEnd(); j++) {
			dChecksums[nodevars[c]] +=
				  dataNode[nodevars[c]][k][i][j]
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
				  dataREdge[redgevars[c]][k][i][j]
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
			double dValue = fabs(dataNode[nodevars[c]][k][i][j]);
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
			double dValue = fabs(dataREdge[redgevars[c]][k][i][j]);
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
			double dValue = dataNode[nodevars[c]][k][i][j];
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
			double dValue = dataREdge[redgevars[c]][k][i][j];
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
			double dValue = dataNode[nodevars[c]][k][i][j];
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
			double dValue = dataREdge[redgevars[c]][k][i][j];
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

	// Vorticity data
	} else if (eDataType == DataType_Divergence) {
		m_connect.Pack(m_dataDivergence);
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
	const DataVector<double> & dCoeff,
	int ixDest,
	DataType eDataType
) {
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

