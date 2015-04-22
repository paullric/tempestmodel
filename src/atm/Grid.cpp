///////////////////////////////////////////////////////////////////////////////
///
///	\file    Grid.cpp
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

#include "Grid.h"
#include "Model.h"
#include "TestCase.h"
#include "ConsolidationStatus.h"
#include "VerticalStretch.h"

#include "Exception.h"

#include <cfloat>
#include <cmath>

#include <netcdfcpp.h>
#include "mpi.h"

///////////////////////////////////////////////////////////////////////////////

Grid::Grid(
	Model & model,
	int nABaseResolution,
	int nBBaseResolution,
	int nRefinementRatio,
	int nRElements,
	VerticalStaggering eVerticalStaggering
) :
	m_fInitialized(false),
	m_iGridStamp(0),
	m_model(model),
	m_fBlockParallelExchange(false),
	m_nABaseResolution(nABaseResolution),
	m_nBBaseResolution(nBBaseResolution),
	m_nRefinementRatio(nRefinementRatio),
	m_dReferenceLength(1.0),
	m_pVerticalStretchF(NULL),
	m_nRElements(nRElements),
	m_dZtop(1.0),
	m_eVerticalStaggering(eVerticalStaggering),
	m_nDegreesOfFreedomPerColumn(0),
	m_fHasReferenceState(false),
	m_fHasRayleighFriction(false)
{
	// Assign a default vertical stretching function
	m_pVerticalStretchF = new VerticalStretchUniform;
}

///////////////////////////////////////////////////////////////////////////////

Grid::~Grid() {
	if (m_pVerticalStretchF != NULL) {
		delete m_pVerticalStretchF;
	}

	for (int n = 0; n < m_vecGridPatches.size(); n++) {
		delete m_vecGridPatches[n];
	}
}

///////////////////////////////////////////////////////////////////////////////

void Grid::SetVerticalStretchFunction(
	VerticalStretchFunction * pVerticalStretchF
) {
	if (m_pVerticalStretchF != NULL) {
		delete m_pVerticalStretchF;
	}

	m_pVerticalStretchF = pVerticalStretchF;
}

///////////////////////////////////////////////////////////////////////////////

void Grid::EvaluateVerticalStretchF(
	double dREta,
	double & dREtaStretch,
	double & dDxREtaStretch
) {
	if (m_pVerticalStretchF == NULL) {
		_EXCEPTIONT("No VerticalStretchFunction defined in Grid");
	}

	(*m_pVerticalStretchF)(dREta, dREtaStretch, dDxREtaStretch);
}

///////////////////////////////////////////////////////////////////////////////

void Grid::InitializeVerticalCoordinate(
	const GridSpacing & gridspacing
) {
	// Number of components in the EquationSet
	int nComponents = m_model.GetEquationSet().GetComponents();

	// Initialize location and index for each variable
	bool fInitializeStaggering = false;
	if (m_vecVarLocation.IsInitialized()) {
		if (m_vecVarLocation.GetRows() != nComponents) {
			_EXCEPTIONT("Mismatch in staggered variable locations and "
				"number of free equations");
		}
	} else {
		m_vecVarLocation.Initialize(nComponents);
		fInitializeStaggering = true;
	}

	// If dimensionality is 2, then initialize a dummy REta
	if (m_model.GetEquationSet().GetDimensionality() == 2) {

		// Resize arrays of REta coordinates
		m_dREtaLevels.Initialize(1);
		m_dREtaInterfaces.Initialize(2);

		m_dREtaLevelsNormArea.Initialize(1);
		m_dREtaInterfacesNormArea.Initialize(2);

		// Uniform grid spacing
		m_dREtaLevels[0] = 0.5;
		m_dREtaInterfaces[0] = 0.0;
		m_dREtaInterfaces[1] = 1.0;

		m_dREtaLevelsNormArea[0] = 1.0;
		m_dREtaInterfacesNormArea[0] = 0.5;
		m_dREtaInterfacesNormArea[1] = 0.5;

		// Everything on model levels
		if (fInitializeStaggering) {
			m_vecVarLocation[0] = DataLocation_Node;
			m_vecVarLocation[1] = DataLocation_Node;
			m_vecVarLocation[2] = DataLocation_Node;
		}

		// Initialize number of degres of freedom per column
		m_nDegreesOfFreedomPerColumn = nComponents;

	// If dimensionality is 3 then initialize normally
	} else if (m_model.GetEquationSet().GetDimensionality() == 3) {

		// Check for agreement with grid spacing
		if (!gridspacing.DoesNodeCountAgree(m_nRElements)) {
			_EXCEPTIONT("Invalid node count for given vertical GridSpacing.");
		}

		// Zero point of GridSpacing object
		double dZeroCoord = gridspacing.GetZeroCoord();

		// Resize arrays of REta coordinates
		m_dREtaLevels.Initialize(m_nRElements);
		m_dREtaInterfaces.Initialize(m_nRElements+1);

		m_dREtaLevelsNormArea.Initialize(m_nRElements);
		m_dREtaInterfacesNormArea.Initialize(m_nRElements+1);

		// Get node/interface location from GridSpacing
		for (int k = 0; k < m_nRElements; k++) {
			m_dREtaLevels[k] = gridspacing.GetNode(k);
			m_dREtaLevelsNormArea[k] = gridspacing.GetNodeNormArea(k);
		}
		for (int k = 0; k <= m_nRElements; k++) {
			m_dREtaInterfaces[k] = gridspacing.GetEdge(k);
			m_dREtaInterfacesNormArea[k] = gridspacing.GetEdgeNormArea(k);
		}

		// Adjust normalized area on edges
		m_dREtaInterfacesNormArea[0] /= 2.0;
		m_dREtaInterfacesNormArea[m_nRElements] /= 2.0;

		if (m_eVerticalStaggering == VerticalStaggering_Interfaces) {
			m_dREtaLevelsNormArea[0] /= 2.0;
			m_dREtaLevelsNormArea[m_nRElements-1] /= 2.0;
		}

		// Location of variables
		if (fInitializeStaggering) {
			switch (m_eVerticalStaggering) {
				case VerticalStaggering_Levels:
					m_vecVarLocation[0] = DataLocation_Node;
					m_vecVarLocation[1] = DataLocation_Node;
					m_vecVarLocation[2] = DataLocation_Node;
					m_vecVarLocation[3] = DataLocation_Node;
					m_vecVarLocation[4] = DataLocation_Node;
					break;

				case VerticalStaggering_Interfaces:
					m_vecVarLocation[0] = DataLocation_Node;
					m_vecVarLocation[1] = DataLocation_Node;
					m_vecVarLocation[2] = DataLocation_Node;
					m_vecVarLocation[3] = DataLocation_Node;
					m_vecVarLocation[4] = DataLocation_Node;
					break;

				case VerticalStaggering_Lorenz:
					m_vecVarLocation[0] = DataLocation_Node;
					m_vecVarLocation[1] = DataLocation_Node;
					m_vecVarLocation[2] = DataLocation_Node;
					m_vecVarLocation[3] = DataLocation_REdge;
					m_vecVarLocation[4] = DataLocation_Node;
					break;

				case VerticalStaggering_CharneyPhillips:
					m_vecVarLocation[0] = DataLocation_Node;
					m_vecVarLocation[1] = DataLocation_Node;
					m_vecVarLocation[2] = DataLocation_REdge;
					m_vecVarLocation[3] = DataLocation_REdge;
					m_vecVarLocation[4] = DataLocation_Node;
					break;

				default:
					_EXCEPTIONT("Invalid VerticalStaggering specified");
			}
		}

		// Initialize number of degres of freedom per column
		m_nDegreesOfFreedomPerColumn = 0;
		for (int c = 0; c < nComponents; c++) {
			if (m_vecVarLocation[c] == DataLocation_Node) {
				m_nDegreesOfFreedomPerColumn += m_nRElements;
			} else if (m_vecVarLocation[c] == DataLocation_REdge) {
				m_nDegreesOfFreedomPerColumn += m_nRElements + 1;
			} else {
				_EXCEPTIONT("(UNIMPLEMENTED) Horizontal staggering");
			}
		}

	} else {
		_EXCEPTIONT("Invalid dimensionality");
	}

	// Calculate stretched values of REta
	if (m_pVerticalStretchF == NULL) {
		m_dREtaStretchLevels = m_dREtaLevels;
		m_dREtaStretchInterfaces = m_dREtaInterfaces;

	} else {
		double dDxREtaStretch;

		m_dREtaStretchLevels.Initialize(m_dREtaLevels.GetRows());
		m_dREtaStretchInterfaces.Initialize(m_dREtaInterfaces.GetRows());

		for (int k = 0; k < m_dREtaLevels.GetRows(); k++) {
			(*m_pVerticalStretchF)(
				m_dREtaLevels[k],
				m_dREtaStretchLevels[k],
				dDxREtaStretch);
		}

		for (int k = 0; k < m_dREtaInterfaces.GetRows(); k++) {
			(*m_pVerticalStretchF)(
				m_dREtaInterfaces[k],
				m_dREtaStretchInterfaces[k],
				dDxREtaStretch);
		}
	}

	// Convert node locations to indices in local arrays
	m_vecVarsAtLocation.resize((int)DataLocation_Count);
	for (int l = 0; l < (int)DataLocation_Count; l++) {
		m_vecVarsAtLocation[l] = 0;
	}

	m_vecVarIndex.resize(m_model.GetEquationSet().GetComponents());
	for (int c = 0; c < m_model.GetEquationSet().GetComponents(); c++) {
		if (m_vecVarLocation[c] == DataLocation_Node) {
			m_vecVarIndex[c] = m_vecVarsAtLocation[(int)DataLocation_Node]++;
		} else if (m_vecVarLocation[c] == DataLocation_AEdge) {
			m_vecVarIndex[c] = m_vecVarsAtLocation[(int)DataLocation_AEdge]++;
		} else if (m_vecVarLocation[c] == DataLocation_BEdge) {
			m_vecVarIndex[c] = m_vecVarsAtLocation[(int)DataLocation_BEdge]++;
		} else if (m_vecVarLocation[c] == DataLocation_REdge) {
			m_vecVarIndex[c] = m_vecVarsAtLocation[(int)DataLocation_REdge]++;
		} else {
			_EXCEPTIONT("Invalid variable location");
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void Grid::EvaluateTestCase(
	const TestCase & test,
	const Time & time,
	int iDataIndex
) {
	// Store the model cap
	m_dZtop = test.GetZtop();

	// Store the reference state flag
	m_fHasReferenceState = test.HasReferenceState();

	// Store the Rayleigh friction flag
	m_fHasRayleighFriction = test.HasRayleighFriction();

	// Evaluate the topography and state/tracer values of the test
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		m_vecActiveGridPatches[n]->
			EvaluateTestCase(test, time, iDataIndex);
	}
}

///////////////////////////////////////////////////////////////////////////////

void Grid::EvaluateTestCase_StateOnly(
	const TestCase & test,
	const Time & time,
	int iDataIndex
) {
	// Evaluate the pointwise values of the test
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		m_vecActiveGridPatches[n]->
			EvaluateTestCase_StateOnly(test, time, iDataIndex);
	}
}

///////////////////////////////////////////////////////////////////////////////

void Grid::ApplyBoundaryConditions(
	int iDataIndex,
	DataType eDataType
) {
	for (int n = 0; n < GetActivePatchCount(); n++) {
		m_vecActiveGridPatches[n]->
			ApplyBoundaryConditions(iDataIndex, eDataType);
	}
}

///////////////////////////////////////////////////////////////////////////////

void Grid::EvaluateGeometricTerms() {
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		m_vecActiveGridPatches[n]->EvaluateGeometricTerms();
	}
}

///////////////////////////////////////////////////////////////////////////////

void Grid::Checksum(
	DataType eDataType,
	DataVector<double> & dChecksums,
	int iDataIndex,
	ChecksumType eChecksumType
) const {

	// Identify root process
	int nRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

	// Initialize the local checksum array from DataType
	DataVector<double> dChecksumsLocal;
	if (eDataType == DataType_State) {
		dChecksumsLocal.Initialize(m_model.GetEquationSet().GetComponents());

	} else if (eDataType == DataType_Tracers) { 
		dChecksumsLocal.Initialize(m_model.GetEquationSet().GetTracers());
		if (m_model.GetEquationSet().GetTracers() == 0) {
			return;
		}

	} else {
		_EXCEPTIONT("Invalid DataType");
	}

	// Loop over all patches and calculate local checksums
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		m_vecActiveGridPatches[n]->Checksum(
			eDataType, dChecksumsLocal, iDataIndex, eChecksumType);
	}

	// Initialize global checksums array at root
	if (nRank == 0) {
		dChecksums.Initialize(dChecksumsLocal.GetRows());
	}

	// Compute sum over all processors and send to root node
	MPI_Op nMPIOperator;
	if (eChecksumType == ChecksumType_Linf) {
		nMPIOperator = MPI_MAX;
	} else {
		nMPIOperator = MPI_SUM;
	}

	MPI_Reduce(
		&(dChecksumsLocal[0]),
		&(dChecksums[0]),
		dChecksumsLocal.GetRows(),
		MPI_DOUBLE,
		nMPIOperator,
		0,
		MPI_COMM_WORLD);

	// Take the square root for the L2 norm sum
	if (nRank == 0) {
		if (eChecksumType == ChecksumType_L2) {
			for (int c = 0; c < dChecksums.GetRows(); c++) {
				dChecksums[c] = sqrt(dChecksums[c]);
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

double Grid::ComputeTotalEnergy(
	int iDataIndex
) const {
	// Compute local energy
	double dLocalEnergy = 0.0;
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		dLocalEnergy +=
			m_vecActiveGridPatches[n]->ComputeTotalEnergy(iDataIndex);
	}

	// Reduce to obtain global energy integral
	double dGlobalEnergy = 0.0;
	MPI_Reduce(
		&dLocalEnergy,
		&dGlobalEnergy,
		1,
		MPI_DOUBLE,
		MPI_SUM,
		0,
		MPI_COMM_WORLD);

	// Return global energy integral
	return dGlobalEnergy;
}

///////////////////////////////////////////////////////////////////////////////

double Grid::ComputeTotalPotentialEnstrophy(
	int iDataIndex
) {
	// Compute vorticity and divergence on the Grid
	ComputeVorticityDivergence(iDataIndex);

	// Compute local enstrophy
	double dLocalPotentialEnstrophy = 0.0;
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		dLocalPotentialEnstrophy +=
			m_vecActiveGridPatches[n]->
				ComputeTotalPotentialEnstrophy(iDataIndex);
	}

	// Reduce to obtain global energy integral
	double dGlobalPotentialEnstrophy = 0.0;
	MPI_Reduce(
		&dLocalPotentialEnstrophy,
		&dGlobalPotentialEnstrophy,
		1,
		MPI_DOUBLE,
		MPI_SUM,
		0,
		MPI_COMM_WORLD);

	// Return global energy integral
	return dGlobalPotentialEnstrophy;
}

///////////////////////////////////////////////////////////////////////////////

void Grid::Exchange(
	DataType eDataType,
	int iDataIndex
) {
	// Block parallel exchanges
	if (m_fBlockParallelExchange) {
		return;
	}

	// Verify all processors are prepared to exchange
	MPI_Barrier(MPI_COMM_WORLD);

	// Set up asynchronous recvs
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		m_vecActiveGridPatches[n]->PrepareExchange();
	}

	// Send data
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		m_vecActiveGridPatches[n]->Send(eDataType, iDataIndex);
	}

	// Receive data
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		m_vecActiveGridPatches[n]->Receive(eDataType, iDataIndex);
	}

	// Wait for send requests to complete
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		m_vecActiveGridPatches[n]->CompleteExchange();
	}

}

///////////////////////////////////////////////////////////////////////////////

void Grid::ExchangeBuffers() {

	// Block parallel exchanges
	if (m_fBlockParallelExchange) {
		return;
	}

	// Verify all processors are prepared to exchange
	MPI_Barrier(MPI_COMM_WORLD);

	// Set up asynchronous recvs
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		m_vecActiveGridPatches[n]->PrepareExchange();
	}

	// Send data
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		m_vecActiveGridPatches[n]->SendBuffers();
	}

	// Receive data
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		m_vecActiveGridPatches[n]->ReceiveBuffers();
	}

	// Wait for send requests to complete
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		m_vecActiveGridPatches[n]->CompleteExchange();
	}
}

///////////////////////////////////////////////////////////////////////////////

void Grid::ExchangeBuffersAndUnpack(
	DataType eDataType,
	int iDataIndex
) {

	// Block parallel exchanges
	if (m_fBlockParallelExchange) {
		return;
	}

	// Verify all processors are prepared to exchange
	MPI_Barrier(MPI_COMM_WORLD);

	// Set up asynchronous recvs
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		m_vecActiveGridPatches[n]->PrepareExchange();
	}

	// Send data
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		m_vecActiveGridPatches[n]->SendBuffers();
	}

	// Receive data
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		m_vecActiveGridPatches[n]->Receive(eDataType, iDataIndex);
	}

	// Wait for send requests to complete
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		m_vecActiveGridPatches[n]->CompleteExchange();
	}
}

///////////////////////////////////////////////////////////////////////////////

int Grid::GetLongestActivePatchPerimeter() const {

	// Longest perimeter
	int nLongestPerimeter = 0;

	// Loop over all active patches
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		int nPerimeter =
			m_vecActiveGridPatches[n]->GetPatchBox().GetInteriorPerimeter();

		if (nPerimeter > nLongestPerimeter) {
			nLongestPerimeter = nPerimeter;
		}
	}

	return nLongestPerimeter;
}

///////////////////////////////////////////////////////////////////////////////

int Grid::GetMaxNodeCount2D() const {

	// Total number of nodes over all patches of grid
	int nMaxNodes2D = 0;

	// Loop over all patches and obtain total node count
	for (int n = 0; n < m_vecGridPatches.size(); n++) {
		int nPatchNodes2D = m_vecGridPatches[n]->GetTotalNodeCount2D();

		if (nPatchNodes2D > nMaxNodes2D) {
			nMaxNodes2D = nPatchNodes2D;
		}
	}

	return nMaxNodes2D;
}

///////////////////////////////////////////////////////////////////////////////

int Grid::GetTotalNodeCount2D() const {

	// Total number of nodes over all patches of grid
	int nTotalNodes2D = 0;

	// Loop over all patches and obtain total node count
	for (int n = 0; n < m_vecGridPatches.size(); n++) {
		nTotalNodes2D += m_vecGridPatches[n]->GetTotalNodeCount2D();
	}

	return nTotalNodes2D;
}

///////////////////////////////////////////////////////////////////////////////

int Grid::GetMaxNodeCount(
	DataLocation loc
) const {

	// Total number of nodes over all patches of grid
	int nMaxNodes = 0;

	// Loop over all patches and obtain total node count
	for (int n = 0; n < m_vecGridPatches.size(); n++) {
		int nPatchNodes = m_vecGridPatches[n]->GetTotalNodeCount(loc);

		if (nPatchNodes > nMaxNodes) {
			nMaxNodes = nPatchNodes;
		}
	}

	return nMaxNodes;
}

///////////////////////////////////////////////////////////////////////////////

int Grid::GetTotalNodeCount(
	DataLocation loc
) const {

	// Total number of nodes over all patches of grid
	int nTotalNodes = 0;

	// Loop over all patches and obtain total node count
	for (int n = 0; n < m_vecGridPatches.size(); n++) {
		nTotalNodes += m_vecGridPatches[n]->GetTotalNodeCount(loc);
	}

	return nTotalNodes;
}

///////////////////////////////////////////////////////////////////////////////

int Grid::GetMaxDegreesOfFreedom() const {

	// Most nodes per patch
	int nMaxDOFs = 0;

	// Loop over all patches and obtain max DOFs from state and tracer data
	for (int n = 0; n < m_vecGridPatches.size(); n++) {
		int nStateDOFs =
			m_vecGridPatches[n]->GetTotalDegreesOfFreedom(
				DataType_State, DataLocation_REdge);

		int nTracersDOFs =
			m_vecGridPatches[n]->GetTotalDegreesOfFreedom(
				DataType_Tracers, DataLocation_REdge);

		if (nTracersDOFs > nMaxDOFs) {
			nMaxDOFs = nTracersDOFs;
		}
		if (nStateDOFs > nMaxDOFs) {
			nMaxDOFs = nStateDOFs;
		}
	}

	return nMaxDOFs;
}

///////////////////////////////////////////////////////////////////////////////

void Grid::ConsolidateDataAtRoot(
	ConsolidationStatus & status,
	DataVector<double> & dataRecvBuffer,
	int & nRecvCount,
	int & ixRecvPatch,
	DataType & eRecvDataType,
	DataLocation & eRecvDataLocation
) const {

	// Get process id
	int nRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

	// Non-root processes should not call this function
	if (nRank != 0) {
		_EXCEPTIONT("Non-root process calling ConsolidateDataAtRoot");
	}

	// Check if done
	if (status.Done()) {
		_EXCEPTIONT("Attempting to consolidate data after completion");
	}

	// Receive a consolidation message
	MPI_Status mpistatus;

	MPI_Recv(
		&(dataRecvBuffer[0]),
		dataRecvBuffer.GetRows(),
		MPI_DOUBLE,
		MPI_ANY_SOURCE,
		MPI_ANY_TAG,
		MPI_COMM_WORLD,
		&mpistatus);

	// Process tag for DataType and global patch index
	ConsolidationStatus::ParseTag(
		mpistatus.MPI_TAG,
		ixRecvPatch,
		eRecvDataType,
		eRecvDataLocation);

	if ((ixRecvPatch < 0) || (ixRecvPatch >= m_vecGridPatches.size())) {
		_EXCEPTIONT("Panel tag index out of range");
	}

	status.SetReceiveStatus(ixRecvPatch, eRecvDataType, eRecvDataLocation);

	// Verify consistency of patch information
	GridPatch * pPatch = m_vecGridPatches[ixRecvPatch];

	MPI_Get_count(&mpistatus, MPI_DOUBLE, &nRecvCount);

	int nExpectedRecvCount =
		pPatch->GetTotalDegreesOfFreedom(eRecvDataType, eRecvDataLocation);

	if (nExpectedRecvCount != nRecvCount) {
		_EXCEPTION2("State dimension mismatch (%i %i)",
			nExpectedRecvCount, nRecvCount);
	}
}

///////////////////////////////////////////////////////////////////////////////

void Grid::ConsolidateDataToRoot(
	ConsolidationStatus & status
) const {
	
	// Get process id
	int nRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

	// If no tracers, Tracer data should not be consolidated
	if ((status.Contains(DataType_Tracers)) &&
		(m_model.GetEquationSet().GetTracers() == 0)
	) {
		_EXCEPTIONT("Attempting to consolidate empty tracer data");
	}

	// Loop over all patches
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		const GridPatch * pPatch = m_vecActiveGridPatches[n];

		// Data
		const GridData4D & dataStateNode =
			pPatch->GetDataState(0, DataLocation_Node);
		const GridData4D & dataStateREdge =
			pPatch->GetDataState(0, DataLocation_REdge);

		const GridData4D & dataRefStateNode =
			pPatch->GetReferenceState(DataLocation_Node);
		const GridData4D & dataRefStateREdge =
			pPatch->GetReferenceState(DataLocation_REdge);

		const GridData4D & dataTracers =
			pPatch->GetDataTracers(0);

		const DataMatrix3D<double> & dataJacobian   = pPatch->GetJacobian();
		const DataMatrix<double>   & dataTopography = pPatch->GetTopography();
		const DataMatrix<double>   & dataLongitude  = pPatch->GetLongitude();
		const DataMatrix<double>   & dataLatitude   = pPatch->GetLatitude();
		const DataMatrix3D<double> & dataZLevels    = pPatch->GetZLevels();

		const GridData3D & dataRayleighStrengthNode =
			pPatch->GetRayleighStrength(DataLocation_Node);
		const GridData3D & dataRayleighStrengthREdge =
			pPatch->GetRayleighStrength(DataLocation_REdge);

		// Send state data on nodes to root process
		if (status.Contains(DataType_State, DataLocation_Node)) {
			MPI_Isend(
				(void*)(dataStateNode[0][0][0]),
				dataStateNode.GetTotalElements(),
				MPI_DOUBLE,
				0,
				ConsolidationStatus::GenerateTag(
					pPatch->GetPatchIndex(),
					DataType_State,
					DataLocation_Node),
				MPI_COMM_WORLD,
				status.GetNextSendRequest());
		}

		// Send state data on radial edges to root process
		if (status.Contains(DataType_State, DataLocation_REdge)) {
			MPI_Isend(
				(void*)(dataStateREdge[0][0][0]),
				dataStateREdge.GetTotalElements(),
				MPI_DOUBLE,
				0,
				ConsolidationStatus::GenerateTag(
					pPatch->GetPatchIndex(),
					DataType_State,
					DataLocation_REdge),
				MPI_COMM_WORLD,
				status.GetNextSendRequest());
		}

		// Send reference state data on nodes to root process
		if (status.Contains(DataType_RefState, DataLocation_Node)) {
			MPI_Isend(
				(void*)(dataRefStateNode[0][0][0]),
				dataRefStateNode.GetTotalElements(),
				MPI_DOUBLE,
				0,
				ConsolidationStatus::GenerateTag(
					pPatch->GetPatchIndex(),
					DataType_RefState,
					DataLocation_Node),
				MPI_COMM_WORLD,
				status.GetNextSendRequest());
		}

		// Send reference state data on radial edges to root process
		if (status.Contains(DataType_RefState, DataLocation_REdge)) {
			MPI_Isend(
				(void*)(dataRefStateREdge[0][0][0]),
				dataRefStateREdge.GetTotalElements(),
				MPI_DOUBLE,
				0,
				ConsolidationStatus::GenerateTag(
					pPatch->GetPatchIndex(),
					DataType_RefState,
					DataLocation_REdge),
				MPI_COMM_WORLD,
				status.GetNextSendRequest());
		}

		// Send tracer data to root process
		if (status.Contains(DataType_Tracers)) {
			MPI_Isend(
				(void*)(dataTracers[0][0][0]),
				dataTracers.GetTotalElements(),
				MPI_DOUBLE,
				0,
				ConsolidationStatus::GenerateTag(
					pPatch->GetPatchIndex(), DataType_Tracers),
				MPI_COMM_WORLD,
				status.GetNextSendRequest());
		}

		// Send Jacobian data to root process
		if (status.Contains(DataType_Jacobian)) {
			MPI_Isend(
				(void*)(dataJacobian[0][0]),
				dataJacobian.GetTotalElements(),
				MPI_DOUBLE,
				0,
				ConsolidationStatus::GenerateTag(
					pPatch->GetPatchIndex(), DataType_Jacobian),
				MPI_COMM_WORLD,
				status.GetNextSendRequest());
		}

		// Send topography data to root process
		if (status.Contains(DataType_Topography)) {
			MPI_Isend(
				(void*)(dataTopography[0]),
				dataTopography.GetTotalElements(),
				MPI_DOUBLE,
				0,
				ConsolidationStatus::GenerateTag(
					pPatch->GetPatchIndex(), DataType_Topography),
				MPI_COMM_WORLD,
				status.GetNextSendRequest());
		}

		// Send longitude data to root process
		if (status.Contains(DataType_Longitude)) {
			MPI_Isend(
				(void*)(dataLongitude[0]),
				dataLongitude.GetTotalElements(),
				MPI_DOUBLE,
				0,
				ConsolidationStatus::GenerateTag(
					pPatch->GetPatchIndex(), DataType_Longitude),
				MPI_COMM_WORLD,
				status.GetNextSendRequest());
		}

		// Send latitude data to root process
		if (status.Contains(DataType_Latitude)) {
			MPI_Isend(
				(void*)(dataLatitude[0]),
				dataLatitude.GetTotalElements(),
				MPI_DOUBLE,
				0,
				ConsolidationStatus::GenerateTag(
					pPatch->GetPatchIndex(), DataType_Latitude),
				MPI_COMM_WORLD,
				status.GetNextSendRequest());
		}

		// Send z-coordinate data to root process
		if (status.Contains(DataType_Z)) {
			MPI_Isend(
				(void*)(dataZLevels[0][0]),
				dataZLevels.GetTotalElements(),
				MPI_DOUBLE,
				0,
				ConsolidationStatus::GenerateTag(
					pPatch->GetPatchIndex(), DataType_Z),
				MPI_COMM_WORLD,
				status.GetNextSendRequest());
		}

		// Send Rayleigh strength data (on nodes) to root process
		if (status.Contains(DataType_RayleighStrength, DataLocation_Node)) {
			MPI_Isend(
				(void*)(dataRayleighStrengthNode[0][0]),
				dataRayleighStrengthNode.GetTotalElements(),
				MPI_DOUBLE,
				0,
				ConsolidationStatus::GenerateTag(
					pPatch->GetPatchIndex(),
					DataType_RayleighStrength,
					DataLocation_Node),
				MPI_COMM_WORLD,
				status.GetNextSendRequest());
		}

		// Send Rayleigh strength data (on interfaces) to root process
		if (status.Contains(DataType_RayleighStrength, DataLocation_REdge)) {
			MPI_Isend(
				(void*)(dataRayleighStrengthREdge[0][0]),
				dataRayleighStrengthREdge.GetTotalElements(),
				MPI_DOUBLE,
				0,
				ConsolidationStatus::GenerateTag(
					pPatch->GetPatchIndex(),
					DataType_RayleighStrength,
					DataLocation_REdge),
				MPI_COMM_WORLD,
				status.GetNextSendRequest());
		}

	}
}

///////////////////////////////////////////////////////////////////////////////

void Grid::ComputeVorticityDivergence(
	int iDataIndex
) {
	// Loop over all grid patches
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		m_vecActiveGridPatches[n]->ComputeVorticityDivergence(iDataIndex);
	}
}

///////////////////////////////////////////////////////////////////////////////

void Grid::ComputeTemperature(
	int iDataIndex
) {
	// Loop over all grid patches
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		m_vecActiveGridPatches[n]->ComputeTemperature(iDataIndex);
	}
}

///////////////////////////////////////////////////////////////////////////////

void Grid::InterpolateNodeToREdge(
	int iVar,
	int iDataIndex
) {
	// Loop over all grid patches
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		m_vecActiveGridPatches[n]->InterpolateNodeToREdge(iVar, iDataIndex);
	}
}

///////////////////////////////////////////////////////////////////////////////

void Grid::InterpolateREdgeToNode(
	int iVar,
	int iDataIndex
) {
	// Loop over all grid patches
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		m_vecActiveGridPatches[n]->InterpolateREdgeToNode(iVar, iDataIndex);
	}
}

///////////////////////////////////////////////////////////////////////////////

void Grid::ReduceInterpolate(
	const DataVector<double> & dAlpha,
	const DataVector<double> & dBeta,
	const DataVector<int> & iPatch,
	DataType eDataType,
	DataLocation eDataLocation,
	bool fInterpAllVariables,
	DataMatrix3D<double> & dInterpData,
	bool fIncludeReferenceState,
	bool fConvertToPrimitive
) const {
	// Check interpolation data array size
	if ((dAlpha.GetRows() != dBeta.GetRows()) ||
		(dAlpha.GetRows() != iPatch.GetRows())
	) {
		_EXCEPTIONT("Inconsistency in vector lengths.");
	}

	if ((eDataType == DataType_Tracers) &&
		(m_model.GetEquationSet().GetTracers() == 0)
	) {
		_EXCEPTIONT("Unable to Interpolate with no tracers.");
	}
	
	// Check interpolation data array size
	if ((eDataType == DataType_State) &&
		(dInterpData.GetRows() != m_model.GetEquationSet().GetComponents())
	) {
		_EXCEPTIONT("InterpData dimension mismatch (0)");
	}

	if ((eDataType == DataType_Tracers) &&
		(dInterpData.GetRows() != m_model.GetEquationSet().GetTracers())
	) {
		_EXCEPTIONT("InterpData dimension mismatch (0)");
	}

	if ((eDataType == DataType_Topography) &&
		(dInterpData.GetRows() != 1)
	) {
		_EXCEPTIONT("InterpData dimension mismatch (0)");
	}

	if ((eDataType == DataType_Vorticity) &&
		(dInterpData.GetRows() != 1)
	) {
		_EXCEPTIONT("InterpData dimension mismatch (0)");
	}

	if ((eDataType == DataType_Divergence) &&
		(dInterpData.GetRows() != 1)
	) {
		_EXCEPTIONT("InterpData dimension mismatch (0)");
	}

	if ((eDataType == DataType_Temperature) &&
		(dInterpData.GetRows() != 1)
	) {
		_EXCEPTIONT("InterpData dimension mismatch (0)");
	}

	if ((eDataLocation == DataLocation_None) &&
		(dInterpData.GetColumns() != 1)
	) {
		_EXCEPTIONT("InterpData dimension mismatch (1)");
	}

	if ((eDataLocation == DataLocation_Node) &&
		(dInterpData.GetColumns() != GetRElements())
	) {
		_EXCEPTIONT("InterpData dimension mismatch (1)");
	}

	if ((eDataLocation == DataLocation_REdge) &&
		(dInterpData.GetColumns() != GetRElements() + 1)
	) {
		_EXCEPTIONT("InterpData dimension mismatch (1)");
	}

	if (dInterpData.GetSubColumns() != dAlpha.GetRows()) {
		_EXCEPTIONT("InterpData dimension mismatch (2)");
	}

	// Zero the interpolated data
	dInterpData.Zero();

	// Interpolate state data
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		m_vecActiveGridPatches[n]->InterpolateData(
			dAlpha, dBeta, iPatch,
			eDataType,
			eDataLocation,
			fInterpAllVariables,
			dInterpData,
			fIncludeReferenceState,
			fConvertToPrimitive);
	}

	// Perform an Reduce operation to combine all data
	int nRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

	if (nRank == 0) {
		MPI_Reduce(
			MPI_IN_PLACE,
			&(dInterpData[0][0][0]),
			dInterpData.GetRows()
				* dInterpData.GetColumns()
				* dInterpData.GetSubColumns(),
			MPI_DOUBLE,
			MPI_SUM,
			0,
			MPI_COMM_WORLD);

	} else {
		MPI_Reduce(
			&(dInterpData[0][0][0]),
			NULL,
			dInterpData.GetRows()
				* dInterpData.GetColumns()
				* dInterpData.GetSubColumns(),
			MPI_DOUBLE,
			MPI_SUM,
			0,
			MPI_COMM_WORLD);
	}
}

///////////////////////////////////////////////////////////////////////////////

void Grid::ConvertReferenceToPatchCoord(
	const DataVector<double> & dXReference,
	const DataVector<double> & dYReference,
	DataVector<double> & dAlpha,
	DataVector<double> & dBeta,
	DataVector<int> & iPatch
) const {
	_EXCEPTIONT("Unimplemented.");
}

///////////////////////////////////////////////////////////////////////////////

GridPatch * Grid::AddPatch(
	GridPatch * pPatch
) {
	int ixNextPatch = m_vecGridPatches.size();

	// Add the patch to the vector of patches
	m_vecGridPatches.push_back(pPatch);

	// Set the patch index
	pPatch->m_ixPatch = ixNextPatch;

	// Update the cumulative 2D index
	if (ixNextPatch == 0) {
		m_vecCumulativePatch2DNodeIndex.push_back(0);
	}

	m_vecCumulativePatch2DNodeIndex.push_back(
		m_vecCumulativePatch2DNodeIndex[ixNextPatch]
		+ pPatch->GetTotalNodeCount2D());

	return pPatch;
}

///////////////////////////////////////////////////////////////////////////////

void Grid::ToFile(
	NcFile & ncfile
) {

	// Output physical constants
	const PhysicalConstants & phys = m_model.GetPhysicalConstants();

	ncfile.add_att("earth_radius", phys.GetEarthRadius());
	ncfile.add_att("g", phys.GetG());
	ncfile.add_att("omega", phys.GetOmega());
	ncfile.add_att("alpha", phys.GetAlpha());
	ncfile.add_att("Rd", phys.GetR());
	ncfile.add_att("Cp", phys.GetCp());
	ncfile.add_att("T0", phys.GetT0());
	ncfile.add_att("P0", phys.GetP0());
	ncfile.add_att("rho_water", phys.GetRhoWater());
	ncfile.add_att("Rvap", phys.GetRvap());
	ncfile.add_att("Mvap", phys.GetMvap());
	ncfile.add_att("Lvap", phys.GetLvap());

	// Create patch index dimension
	NcDim * dimPatchIndex =
		ncfile.add_dim("patch_index", GetPatchCount());

	// Get the length of each dimension array
	int iANodeIndex = 0;
	int iBNodeIndex = 0;
	int iAEdgeIndex = 0;
	int iBEdgeIndex = 0;

	int nANodeCount = 0;
	int nBNodeCount = 0;
	int nAEdgeCount = 0;
	int nBEdgeCount = 0;

	for (int n = 0; n < GetPatchCount(); n++) {
		const PatchBox & box = GetPatch(n)->GetPatchBox();

		nANodeCount += box.GetANodes().GetRows();
		nBNodeCount += box.GetBNodes().GetRows();
		nAEdgeCount += box.GetAEdges().GetRows();
		nBEdgeCount += box.GetBEdges().GetRows();
	}

	NcDim * dimANodeCount =
		ncfile.add_dim("alpha_node_index", nANodeCount);
	NcDim * dimBNodeCount =
		ncfile.add_dim("beta_node_index", nBNodeCount);
	NcDim * dimAEdgeCount =
		ncfile.add_dim("alpha_edge_index", nAEdgeCount);
	NcDim * dimBEdgeCount =
		ncfile.add_dim("beta_edge_index", nBEdgeCount);

	NcVar * varANodeCoord =
		ncfile.add_var(
			"alpha_node_coord", ncDouble, dimANodeCount);
	NcVar * varBNodeCoord =
		ncfile.add_var(
			"beta_node_coord", ncDouble, dimBNodeCount);
	NcVar * varAEdgeCoord =
		ncfile.add_var(
			"alpha_edge_coord", ncDouble, dimAEdgeCount);
	NcVar * varBEdgeCoord =
		ncfile.add_var(
			"beta_edge_coord", ncDouble, dimBEdgeCount);

	// Output integer grid information
	int iGridInfo[6];
	iGridInfo[0] = m_iGridStamp;
	iGridInfo[1] = m_nABaseResolution;
	iGridInfo[2] = m_nBBaseResolution;
	iGridInfo[3] = m_nRefinementRatio;
	iGridInfo[4] = (int)(m_fHasReferenceState);
	iGridInfo[5] = (int)(m_fHasRayleighFriction);

	NcDim * dimGridInfoCount =
		ncfile.add_dim("grid_info_count", 6);

	NcVar * varGridInfo =
		ncfile.add_var("grid_info", ncInt, dimGridInfoCount);

	varGridInfo->put(iGridInfo, 6);

	// Output floating point grid information
	double dGridInfo[2];
	dGridInfo[0] = m_dReferenceLength;
	dGridInfo[1] = m_dZtop;

	NcDim * dimGridInfoFloat =
		ncfile.add_dim("grid_info_float_count", 2);

	NcVar * varGridInfoFloat =
		ncfile.add_var("grid_info_float", ncDouble, dimGridInfoFloat);

	varGridInfoFloat->put(dGridInfo, 2);

	// Output REta coordinate of levels
	NcDim * dimREtaLevels =
		ncfile.add_dim("reta_levels", m_dREtaLevels.GetRows());

	NcVar * varREtaLevels =
		ncfile.add_var("reta_levels", ncDouble, dimREtaLevels);

	varREtaLevels->put(
		&(m_dREtaStretchLevels[0]), m_dREtaStretchLevels.GetRows());

	// Output REta coordinate of interfaces
	NcDim * dimREtaInterfaces =
		ncfile.add_dim("reta_interfaces", m_dREtaInterfaces.GetRows());

	NcVar * varREtaInterfaces =
		ncfile.add_var("reta_interfaces", ncDouble, dimREtaInterfaces);

	varREtaInterfaces->put(
		&(m_dREtaStretchInterfaces[0]), m_dREtaStretchInterfaces.GetRows());
/*
	// Output location of each variable
	NcDim * dimVarLocation =
		ncfile.add_dim("component", m_vecVarLocation.size());

	NcVar * varVarLocation =
		ncfile.add_var("var_location", ncInt, dimVarLocation);

	varVarLocation->put(
		(int*)(&(m_vecVarLocation[0])), m_vecVarLocation.size());
*/
	// Output PatchBox for each patch
	NcDim * dimPatchInfoCount =
		ncfile.add_dim("patch_info_count", 7);

	NcVar * varPatchInfo =
		ncfile.add_var("patch_info", ncInt, dimPatchIndex, dimPatchInfoCount);

	for (int n = 0; n < GetPatchCount(); n++) {
		int iPatchInfo[7];

		const PatchBox & box = GetPatch(n)->GetPatchBox();

		iPatchInfo[0] = box.GetPanel();
		iPatchInfo[1] = box.GetRefinementLevel();
		iPatchInfo[2] = box.GetHaloElements();
		iPatchInfo[3] = box.GetAGlobalInteriorBegin();
		iPatchInfo[4] = box.GetAGlobalInteriorEnd();
		iPatchInfo[5] = box.GetBGlobalInteriorBegin();
		iPatchInfo[6] = box.GetBGlobalInteriorEnd();

		varPatchInfo->set_cur(n, 0);
		varPatchInfo->put(iPatchInfo, 1, 7);

		varANodeCoord->set_cur(iANodeIndex);
		varBNodeCoord->set_cur(iBNodeIndex);
		varAEdgeCoord->set_cur(iAEdgeIndex);
		varBEdgeCoord->set_cur(iBEdgeIndex);

		varANodeCoord->put(box.GetANodes(), box.GetANodes().GetRows());
		varAEdgeCoord->put(box.GetAEdges(), box.GetAEdges().GetRows());
		varBNodeCoord->put(box.GetBNodes(), box.GetBNodes().GetRows());
		varBEdgeCoord->put(box.GetBEdges(), box.GetBEdges().GetRows());

		iANodeIndex += box.GetANodes().GetRows();
		iAEdgeIndex += box.GetAEdges().GetRows();
		iBNodeIndex += box.GetBNodes().GetRows();
		iBEdgeIndex += box.GetBEdges().GetRows();
	}
}

///////////////////////////////////////////////////////////////////////////////

void Grid::FromFile(
	const std::string & strGridFile
) {
	// Check for existing patches
	if (GetPatchCount() != 0) {
		_EXCEPTIONT("Trying to load over non-empty grid");
	}

	// Open the NetCDF file
	NcFile ncfile(strGridFile.c_str(), NcFile::ReadOnly);
	if (!ncfile.is_valid()) {
		_EXCEPTION1("Invalid restart file: %s", strGridFile.c_str());
	}

	// Input physical constants
	PhysicalConstants & phys = m_model.GetPhysicalConstants();

	NcAtt * attPhys;

	attPhys = ncfile.get_att("earth_radius");
	phys.SetEarthRadius(attPhys->as_double(0));

	attPhys = ncfile.get_att("g");
	phys.SetG(attPhys->as_double(0));

	attPhys = ncfile.get_att("omega");
	phys.SetOmega(attPhys->as_double(0));

	attPhys = ncfile.get_att("alpha");
	phys.SetAlpha(attPhys->as_double(0));

	attPhys = ncfile.get_att("Rd");
	phys.SetR(attPhys->as_double(0));

	attPhys = ncfile.get_att("Cp");
	phys.SetCp(attPhys->as_double(0));

	attPhys = ncfile.get_att("T0");
	phys.SetT0(attPhys->as_double(0));

	attPhys = ncfile.get_att("P0");
	phys.SetP0(attPhys->as_double(0));

	attPhys = ncfile.get_att("rho_water");
	phys.SetRhoWater(attPhys->as_double(0));

	attPhys = ncfile.get_att("Rvap");
	phys.SetRvap(attPhys->as_double(0));

	attPhys = ncfile.get_att("Mvap");
	phys.SetMvap(attPhys->as_double(0));

	attPhys = ncfile.get_att("Lvap");
	phys.SetLvap(attPhys->as_double(0));

	// Load in grid info
	int iGridInfo[6];

	NcVar * varGridInfo = ncfile.get_var("grid_info");

	varGridInfo->get(iGridInfo, 6);

	m_iGridStamp = iGridInfo[0];
	m_nABaseResolution = iGridInfo[1];
	m_nBBaseResolution = iGridInfo[2];
	m_nRefinementRatio = iGridInfo[3];
	m_fHasReferenceState = (bool)(iGridInfo[4]);
	m_fHasRayleighFriction = (bool)(iGridInfo[5]);

	// Load in floating point grid info
	double dGridInfo[2];

	NcVar * varGridInfoFloat = ncfile.get_var("grid_info_float");

	varGridInfoFloat->get(dGridInfo, 2);

	m_dReferenceLength = dGridInfo[0];
	m_dZtop = dGridInfo[1];

	// Load in location of REta levels
	NcVar * varREtaLevels = ncfile.get_var("reta_levels");

	NcDim * dimREtaLevels = varREtaLevels->get_dim(0);
	int nREtaLevels = static_cast<int>(dimREtaLevels->size());

	m_dREtaStretchLevels.Initialize(nREtaLevels);
	varREtaLevels->get(m_dREtaStretchLevels, nREtaLevels);

	// Load in location of REta interfaces
	NcVar * varREtaInterfaces = ncfile.get_var("reta_interfaces");

	NcDim * dimREtaInterfaces = varREtaInterfaces->get_dim(0);
	int nREtaInterfaces = static_cast<int>(dimREtaInterfaces->size());

	m_dREtaStretchInterfaces.Initialize(nREtaInterfaces);
	varREtaInterfaces->get(m_dREtaStretchInterfaces, nREtaInterfaces);
/*
	// Load in location of variables
	NcVar * varVarLocation = ncfile.get_var("var_location");

	NcDim * dimVarLocation = varVarLocation->get_dim(0);
	int nComponents = static_cast<int>(dimVarLocation->size());

	m_vecVarLocation.resize(nComponents);
	varVarLocation->get((int*)(&(m_vecVarLocation[0])), nComponents);
*/
	// Coordinate arrays
	int iANodeIndex = 0;
	int iBNodeIndex = 0;
	int iAEdgeIndex = 0;
	int iBEdgeIndex = 0;

	DataVector<double> dANodes;
	DataVector<double> dBNodes;
	DataVector<double> dAEdges;
	DataVector<double> dBEdges;

	// Load in all PatchBoxes
	NcVar * varPatchInfo = ncfile.get_var("patch_info");
	if (varPatchInfo == NULL) {
		_EXCEPTIONT("Invalid GridFile; variable patch_info required");
	}

	NcVar * varANodeCoord = ncfile.get_var("alpha_node_coord");
	if (varANodeCoord == NULL) {
		_EXCEPTIONT("Invalid GridFile; variable alpha_node_coord required");
	}

	NcVar * varBNodeCoord = ncfile.get_var("beta_node_coord");
	if (varBNodeCoord == NULL) {
		_EXCEPTIONT("Invalid GridFile; variable beta_node_coord required");
	}

	NcVar * varAEdgeCoord = ncfile.get_var("alpha_edge_coord");
	if (varAEdgeCoord == NULL) {
		_EXCEPTIONT("Invalid GridFile; variable alpha_edge_coord required");
	}

	NcVar * varBEdgeCoord = ncfile.get_var("beta_edge_coord");
	if (varBEdgeCoord == NULL) {
		_EXCEPTIONT("Invalid GridFile; variable beta_edge_coord required");
	}

	NcDim * dimPatchInfoCount = varPatchInfo->get_dim(0);
	int nPatches = static_cast<int>(dimPatchInfoCount->size());

	for (int ix = 0; ix < nPatches; ix++) {
		int iPatchInfo[7];
		varPatchInfo->set_cur(ix, 0);
		varPatchInfo->get(iPatchInfo, 1, 7);

		int nANodes = iPatchInfo[4] - iPatchInfo[3] + 2 * iPatchInfo[2];
		int nBNodes = iPatchInfo[6] - iPatchInfo[5] + 2 * iPatchInfo[2];

		dANodes.Initialize(nANodes);
		dBNodes.Initialize(nBNodes);
		dAEdges.Initialize(nANodes+1);
		dBEdges.Initialize(nBNodes+1);

		varANodeCoord->set_cur(iANodeIndex);
		varBNodeCoord->set_cur(iBNodeIndex);
		varAEdgeCoord->set_cur(iAEdgeIndex);
		varBEdgeCoord->set_cur(iBEdgeIndex);

		varANodeCoord->get(dANodes, nANodes);
		varBNodeCoord->get(dBNodes, nBNodes);
		varAEdgeCoord->get(dAEdges, nANodes+1);
		varBEdgeCoord->get(dBEdges, nBNodes+1);

		PatchBox box(
			iPatchInfo[0],
			iPatchInfo[1],
			iPatchInfo[2],
			iPatchInfo[3],
			iPatchInfo[4],
			iPatchInfo[5],
			iPatchInfo[6],
			dANodes,
			dBNodes,
			dAEdges,
			dBEdges);

		AddPatch(ix, box);

		iANodeIndex += box.GetANodes().GetRows();
		iBNodeIndex += box.GetBNodes().GetRows();
		iAEdgeIndex += box.GetAEdges().GetRows();
		iBEdgeIndex += box.GetBEdges().GetRows();
	}
}

///////////////////////////////////////////////////////////////////////////////

void Grid::DistributePatches() {

	// Number of processors
	int nSize;
	MPI_Comm_size(MPI_COMM_WORLD, &nSize);

	// Current processor
	int nRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

	// Loop over all patches and initialize data
	for (int n = 0; n < m_vecGridPatches.size(); n++) {
		int nPatchProcessor = n % nSize;
		if (nPatchProcessor == nRank) {
			m_vecGridPatches[n]->InitializeDataLocal();
			m_vecActiveGridPatches.push_back(m_vecGridPatches[n]);
		} else {
			m_vecGridPatches[n]->InitializeDataRemote(nPatchProcessor);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void Grid::InitializeConnectivity() {

	// Vector of nodal points around element
	DataVector<int> vecIxA;
	DataVector<int> vecIxB;
	DataVector<int> vecPanel;
	DataVector<int> vecPatch;

	// Determine longest perimeter
	int nLongestActivePerimeter = GetLongestActivePatchPerimeter() + 4;
	vecIxA.Initialize(nLongestActivePerimeter);
	vecIxB.Initialize(nLongestActivePerimeter);
	vecPanel.Initialize(nLongestActivePerimeter);
	vecPatch.Initialize(nLongestActivePerimeter);

	// Loop over all active patches
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {

		GridPatch * pPatch = m_vecActiveGridPatches[n];

		const PatchBox & box = pPatch->GetPatchBox();

		int ix = 0;

		// Bottom-left corner
		vecIxA[ix] = box.GetAGlobalInteriorBegin()-1;
		vecIxB[ix] = box.GetBGlobalInteriorBegin()-1;
		vecPanel[ix] = box.GetPanel();
		ix++;

		// Bottom edge
		for (int i = box.GetAGlobalInteriorBegin();
		         i < box.GetAGlobalInteriorEnd(); i++
		) {
			vecIxA[ix] = i;
			vecIxB[ix] = box.GetBGlobalInteriorBegin()-1;
			vecPanel[ix] = box.GetPanel();
			ix++;
		}

		// Bottom-right corner
		vecIxA[ix] = box.GetAGlobalInteriorEnd();
		vecIxB[ix] = box.GetBGlobalInteriorBegin()-1;
		vecPanel[ix] = box.GetPanel();
		ix++;

		// Right edge
		for (int j = box.GetBGlobalInteriorBegin();
		         j < box.GetBGlobalInteriorEnd(); j++
		) {
			vecIxA[ix] = box.GetAGlobalInteriorEnd();
			vecIxB[ix] = j;
			vecPanel[ix] = box.GetPanel();
			ix++;
		}

		// Top-right corner
		vecIxA[ix] = box.GetAGlobalInteriorEnd();
		vecIxB[ix] = box.GetBGlobalInteriorEnd();
		vecPanel[ix] = box.GetPanel();
		ix++;

		// Top edge
		for (int i = box.GetAGlobalInteriorEnd()-1;
		         i >= box.GetAGlobalInteriorBegin(); i--
		) {
			vecIxA[ix] = i;
			vecIxB[ix] = box.GetBGlobalInteriorEnd();
			vecPanel[ix] = box.GetPanel();
			ix++;
		}

		// Top-left corner
		vecIxA[ix] = box.GetAGlobalInteriorBegin()-1;
		vecIxB[ix] = box.GetBGlobalInteriorEnd();
		vecPanel[ix] = box.GetPanel();
		ix++;

		// Left edge
		for (int j = box.GetBGlobalInteriorEnd()-1;
		         j >= box.GetBGlobalInteriorBegin(); j--
		) {
			vecIxA[ix] = box.GetAGlobalInteriorBegin()-1;
			vecIxB[ix] = j;
			vecPanel[ix] = box.GetPanel();
			ix++;
		}

		// Get neighboring patches at each halo node
		GetPatchFromCoordinateIndex(
			box.GetRefinementLevel(),
			vecIxA,
			vecIxB,
			vecPanel,
			vecPatch,
			ix);

		// Verify index length
		if (ix != box.GetInteriorPerimeter() + 4) {
			_EXCEPTIONT("Index mismatch");
		}

		// Reset index
		ix = 0;

		// Add connectivity to bottom-left corner
		if (vecPatch[ix] != GridPatch::InvalidIndex) {
			pPatch->ExteriorConnect(
				Direction_BottomLeft,
				m_vecGridPatches[vecPatch[ix]]);
		}

		ix++;

		// Add connectivity to bottom edge: Look for segments along each
		// edge that connect to distinct patches and construct corresponding
		// ExteriorNeighbors.
		{
			int ixFirstBegin = box.GetAInteriorBegin();
			int iCurrentPatch = vecPatch[ix];

			for (int i = ixFirstBegin; i <= box.GetAInteriorEnd(); i++) {
				if ((i == box.GetAInteriorEnd()) ||
					(vecPatch[ix] != iCurrentPatch)
				) {
					const GridPatch * pPatchBottom
						= m_vecGridPatches[iCurrentPatch];

					pPatch->ExteriorConnect(
						Direction_Bottom,
						pPatchBottom,
						ixFirstBegin,
						i);

					if (i != box.GetAInteriorEnd()) {
						pPatch->ExteriorConnect(
							Direction_BottomLeft,
							pPatchBottom,
							i,
							box.GetBInteriorBegin());

						ixFirstBegin = i;
						iCurrentPatch = vecPatch[ix];
						pPatchBottom = m_vecGridPatches[iCurrentPatch];

						pPatch->ExteriorConnect(
							Direction_BottomRight,
							pPatchBottom,
							i - 1,
							box.GetBInteriorBegin());
					}
				}
				if (i != box.GetAInteriorEnd()) {
					ix++;
				}
			}
		}

		// Add connectivity to bottom-right corner
		if (vecPatch[ix] != GridPatch::InvalidIndex) {
			pPatch->ExteriorConnect(
				Direction_BottomRight,
				m_vecGridPatches[vecPatch[ix]]);
		}

		ix++;

		// Add connectivity to right edge
		{
			int ixFirstBegin = box.GetBInteriorBegin();
			int iCurrentPatch = vecPatch[ix];

			for (int j = ixFirstBegin; j <= box.GetBInteriorEnd(); j++) {
				if ((j == box.GetBInteriorEnd()) ||
					(vecPatch[ix] != iCurrentPatch)
				) {
					const GridPatch * pPatchRight
						= m_vecGridPatches[iCurrentPatch];

					pPatch->ExteriorConnect(
						Direction_Right,
						pPatchRight,
						ixFirstBegin,
						j);

					if (j != box.GetBInteriorEnd()) {
						pPatch->ExteriorConnect(
							Direction_BottomRight,
							pPatchRight,
							box.GetAInteriorEnd()-1,
							j);

						ixFirstBegin = j;
						iCurrentPatch = vecPatch[ix];
						pPatchRight = m_vecGridPatches[iCurrentPatch];

						pPatch->ExteriorConnect(
							Direction_TopRight,
							pPatchRight,
							box.GetAInteriorEnd()-1,
							j - 1);
					}
				}
				if (j != box.GetBInteriorEnd()) {
					ix++;
				}
			}
		}

		// Add connectivity to top-right corner
		if (vecPatch[ix] != GridPatch::InvalidIndex) {
			pPatch->ExteriorConnect(
				Direction_TopRight,
				m_vecGridPatches[vecPatch[ix]]);
		}

		ix++;

		// Add connectivity to top edge
		{
			int ixFirstEnd = box.GetAInteriorEnd();
			int iCurrentPatch = vecPatch[ix];

			for (int i = ixFirstEnd-1; i >= box.GetAInteriorBegin()-1; i--) {
				if ((i == box.GetAInteriorBegin()-1) ||
					(vecPatch[ix] != iCurrentPatch)
				) {
					const GridPatch * pPatchTop =
						m_vecGridPatches[iCurrentPatch];

					pPatch->ExteriorConnect(
						Direction_Top,
						pPatchTop,
						i + 1,
						ixFirstEnd);

					if (i != box.GetAInteriorBegin()-1) {
						pPatch->ExteriorConnect(
							Direction_TopRight,
							pPatchTop,
							i,
							box.GetBInteriorEnd()-1);

						ixFirstEnd = i + 1;
						iCurrentPatch = vecPatch[ix];
						pPatchTop = m_vecGridPatches[iCurrentPatch];

						pPatch->ExteriorConnect(
							Direction_TopLeft,
							pPatchTop,
							i + 1,
							box.GetBInteriorEnd()-1);
					}
				}
				if (i != box.GetAInteriorBegin()-1) {
					ix++;
				}
			}
		}

		// Add connectivity to top-left corner
		if (vecPatch[ix] != GridPatch::InvalidIndex) {
			pPatch->ExteriorConnect(
				Direction_TopLeft,
				m_vecGridPatches[vecPatch[ix]]);
		}

		ix++;

		// Add connectivity to top edge
		{
			int ixFirstEnd = box.GetBInteriorEnd();
			int iCurrentPatch = vecPatch[ix];

			for (int j = ixFirstEnd-1; j >= box.GetBInteriorBegin()-1; j--) {
				if ((j == box.GetBInteriorBegin()-1) ||
					(vecPatch[ix] != iCurrentPatch)
				) {
					const GridPatch * pPatchLeft =
						m_vecGridPatches[iCurrentPatch];

					pPatch->ExteriorConnect(
						Direction_Left,
						pPatchLeft,
						j + 1,
						ixFirstEnd);

					if (j != box.GetBInteriorBegin()-1) {
						pPatch->ExteriorConnect(
							Direction_TopLeft,
							pPatchLeft,
							box.GetAInteriorBegin(),
							j);

						ixFirstEnd = j + 1;
						iCurrentPatch = vecPatch[ix];
						pPatchLeft = m_vecGridPatches[iCurrentPatch];

						pPatch->ExteriorConnect(
							Direction_BottomLeft,
							pPatchLeft,
							box.GetAInteriorBegin(),
							j + 1);
					}
				}
				if (j != box.GetBInteriorBegin()-1) {
					ix++;
				}
			}
		}

		// Verify index length
		if (ix != box.GetInteriorPerimeter() + 4) {
			_EXCEPTIONT("Index mismatch");
		}

		// Initialize flux connectivity for patch
		pPatch->GetConnectivity().BuildFluxConnectivity();
	}
}

///////////////////////////////////////////////////////////////////////////////

void Grid::CopyData(
	int ixSource,
	int ixDest,
	DataType eDataType
) {
	// Loop over all grid patches
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		m_vecActiveGridPatches[n]->
			CopyData(ixSource, ixDest, eDataType);
	}
}

///////////////////////////////////////////////////////////////////////////////

void Grid::LinearCombineData(
	const DataVector<double> & dCoeff,
	int ixDest,
	DataType eDataType
) {
	// Loop over all grid patches
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		m_vecActiveGridPatches[n]->
			LinearCombineData(dCoeff, ixDest, eDataType);
	}
}

///////////////////////////////////////////////////////////////////////////////

void Grid::ZeroData(
	int ixData,
	DataType eDataType
) {
	// Loop over all grid patches
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		m_vecActiveGridPatches[n]->ZeroData(ixData, eDataType);
	}
}

///////////////////////////////////////////////////////////////////////////////

void Grid::AddReferenceState(
	int ix
) {
	// Loop over all grid patches
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		m_vecActiveGridPatches[n]->AddReferenceState(ix);
	}
}

///////////////////////////////////////////////////////////////////////////////

