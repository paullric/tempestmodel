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
#include "GridSpacing.h"
#include "VerticalStretch.h"
#include "ConsolidationStatus.h"

#include "Exception.h"

#include <cfloat>
#include <cmath>

#ifndef NO_NETCDF
#include <netcdfcpp.h>
#endif

#ifdef USE_MPI
#include "mpi.h"
#endif

///////////////////////////////////////////////////////////////////////////////

Grid::Grid(
	Model & model
) :
	m_fInitialized(false),
	m_model(model),
	m_fBlockParallelExchange(false),
	m_pVerticalStretchF(NULL)
{ }

///////////////////////////////////////////////////////////////////////////////

void Grid::DefineParameters() {
	if (m_dcGridParameters.IsAttached()) {
		_EXCEPTIONT("Attempting to recall DefineParameters");
	}
#pragma message "May have to modify DataContainer to incorporate padding for DataChunks"

	// Four lateral boundaries (rectangular mesh)
	m_eBoundaryCondition.SetSize(4);

	// Initialize the GridParameters DataContainer
	m_dcGridParameters.PushDataChunk(&m_nMaxPatchCount);
	m_dcGridParameters.PushDataChunk(&m_nRElements);
	m_dcGridParameters.PushDataChunk(&m_iGridStamp);
	m_dcGridParameters.PushDataChunk(&m_eVerticalStaggering);
	m_dcGridParameters.PushDataChunk(&m_nABaseResolution);
	m_dcGridParameters.PushDataChunk(&m_nBBaseResolution);
	m_dcGridParameters.PushDataChunk(&m_nRefinementRatio);
	m_dcGridParameters.PushDataChunk(&m_eBoundaryCondition);
	m_dcGridParameters.PushDataChunk(&m_dReferenceLength);
	m_dcGridParameters.PushDataChunk(&m_dZtop);
	m_dcGridParameters.PushDataChunk(&m_fHasReferenceState);
	m_dcGridParameters.PushDataChunk(&m_fHasRayleighFriction);

	m_dcGridParameters.Allocate();
}

///////////////////////////////////////////////////////////////////////////////

void Grid::SetParameters(
	int nRElements,
	int nMaxPatchCount,
	int nABaseResolution,
	int nBBaseResolution,
	int nRefinementRatio,
	VerticalStaggering eVerticalStaggering
) {
	if (!m_dcGridParameters.IsAttached()) {
		_EXCEPTIONT("DefineParameters() must be called before SetParameters()");
	}

	// Default GridParameters values
	m_nMaxPatchCount = nMaxPatchCount;
	m_nRElements = nRElements;
	m_iGridStamp = 0;
	m_eVerticalStaggering = eVerticalStaggering;
	m_nABaseResolution = nABaseResolution;
	m_nBBaseResolution = nBBaseResolution;
	m_nRefinementRatio = nRefinementRatio;
	m_dReferenceLength = 1.0;
	m_dZtop = 1.0;
	m_fHasReferenceState = false;
	m_fHasRayleighFriction = false;

	// Set boundary conditions
	for (int i = 0; i < 4; i++) {
		m_eBoundaryCondition[i] = BoundaryCondition_Default;
	}
}

///////////////////////////////////////////////////////////////////////////////

void Grid::InitializeDataLocal() {

	if (m_dcGridPatchData.IsAttached()) {
		_EXCEPTIONT("Attempting to call InitializeDataLocal on"
			" initialized Grid");
	}

	// Number of components
	int nComponents = m_model.GetEquationSet().GetComponents();

	// Set the size of all DataArray objects
	m_aPatchBoxes.SetSize(m_nMaxPatchCount);
	m_dREtaLevels.SetSize(m_nRElements);
	m_dREtaInterfaces.SetSize(m_nRElements+1);
	m_dREtaLevelsNormArea.SetSize(m_nRElements);
	m_dREtaInterfacesNormArea.SetSize(m_nRElements+1);
	m_dREtaStretchLevels.SetSize(m_nRElements);
	m_dREtaStretchInterfaces.SetSize(m_nRElements+1);
	m_vecVarLocation.SetSize(nComponents);
	m_vecVarIndex.SetSize(nComponents);
	m_vecVarsAtLocation.SetSize((size_t)DataLocation_Count);

	// Initialize the GridData DataContainer
	m_dcGridPatchData.PushDataChunk(&m_nInitializedPatchBoxes);
	m_dcGridPatchData.PushDataChunk(&m_aPatchBoxes);
	m_dcGridPatchData.PushDataChunk(&m_dREtaLevels);
	m_dcGridPatchData.PushDataChunk(&m_dREtaInterfaces);
	m_dcGridPatchData.PushDataChunk(&m_dREtaLevelsNormArea);
	m_dcGridPatchData.PushDataChunk(&m_dREtaInterfacesNormArea);
	m_dcGridPatchData.PushDataChunk(&m_dREtaStretchLevels);
	m_dcGridPatchData.PushDataChunk(&m_dREtaStretchInterfaces);
	m_dcGridPatchData.PushDataChunk(&m_vecVarLocation);
	m_dcGridPatchData.PushDataChunk(&m_vecVarIndex);
	m_dcGridPatchData.PushDataChunk(&m_vecVarsAtLocation);

	m_dcGridPatchData.Allocate();

	// Default GridData values
	m_nInitializedPatchBoxes = 0;
}

///////////////////////////////////////////////////////////////////////////////

Grid::~Grid() {
	if (m_pVerticalStretchF != NULL) {
		delete m_pVerticalStretchF;
	}

	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		delete m_vecActiveGridPatches[n];
	}
}

///////////////////////////////////////////////////////////////////////////////

void Grid::SetBoundaryCondition(
	Direction eDir,
	BoundaryCondition eBoundaryCondition
) {
	int iDir = static_cast<int>(eDir);

	if ((iDir < 0) || (iDir > 3)) {
		_EXCEPTIONT("Invalid Direction specified");
	}

	m_eBoundaryCondition[iDir] = eBoundaryCondition;
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

void Grid::InitializeVerticalCoordinate() {

	// Number of components in the EquationSet
	int nComponents = m_model.GetEquationSet().GetComponents();

	// Check array allocation
	if (!m_vecVarLocation.IsAttached()) {
		_EXCEPTIONT("vecVarLocation not initializead");
	}
	if (!m_vecVarIndex.IsAttached()) {
		_EXCEPTIONT("vecVarIndex not initialized");
	}
	if (!m_dREtaLevels.IsAttached()) {
		_EXCEPTIONT("dREtaLevels not initialized");
	}
	if (!m_dREtaInterfaces.IsAttached()) {
		_EXCEPTIONT("dREtaInterfaces not initialized");
	}
	if (!m_dREtaStretchLevels.IsAttached()) {
		_EXCEPTIONT("dREtaStretchLevels not initialized");
	}
	if (!m_dREtaStretchInterfaces.IsAttached()) {
		_EXCEPTIONT("dREtaStretchInterfaces not initialized");
	}

	// If dimensionality is 2, then initialize a dummy REta
	if (m_model.GetEquationSet().GetDimensionality() == 2) {

		// Uniform grid spacing
		m_dREtaLevels[0] = 0.5;
		m_dREtaInterfaces[0] = 0.0;
		m_dREtaInterfaces[1] = 1.0;

		m_dREtaLevelsNormArea[0] = 1.0;
		m_dREtaInterfacesNormArea[0] = 0.5;
		m_dREtaInterfacesNormArea[1] = 0.5;

		// Everything on model levels
		m_vecVarLocation[0] = DataLocation_Node;
		m_vecVarLocation[1] = DataLocation_Node;
		m_vecVarLocation[2] = DataLocation_Node;

	// If dimensionality is 3 then initialize normally
	} else if (m_model.GetEquationSet().GetDimensionality() == 3) {

		// Location of variables
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

		// Default to uniform GridSpacing
		double dAvgDeltaElement = 1.0 / static_cast<double>(m_nRElements);
		double dZeroCoord = 0.0;

		GridSpacingUniform gridspacing(dAvgDeltaElement, dZeroCoord);

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

	} else {
		_EXCEPTIONT("Invalid dimensionality");
	}
/*
	// Calculate stretched values of REta
	if (m_pVerticalStretchF == NULL) {
		m_dREtaStretchLevels = m_dREtaLevels;
		m_dREtaStretchInterfaces = m_dREtaInterfaces;

	} else {
		double dDxREtaStretch;

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
*/
	// Convert node locations to indices in local arrays
	for (size_t l = 0; l < (size_t)DataLocation_Count; l++) {
		m_vecVarsAtLocation[l] = 0;
	}

	//m_vecVarIndex.Allocate(m_model.GetEquationSet().GetComponents());
	for (size_t c = 0; c < m_model.GetEquationSet().GetComponents(); c++) {
		if (m_vecVarLocation[c] == DataLocation_Node) {
			m_vecVarIndex[c] =
				m_vecVarsAtLocation[(size_t)DataLocation_Node]++;

		} else if (m_vecVarLocation[c] == DataLocation_AEdge) {
			m_vecVarIndex[c] =
				m_vecVarsAtLocation[(size_t)DataLocation_AEdge]++;

		} else if (m_vecVarLocation[c] == DataLocation_BEdge) {
			m_vecVarIndex[c] =
				m_vecVarsAtLocation[(size_t)DataLocation_BEdge]++;

		} else if (m_vecVarLocation[c] == DataLocation_REdge) {
			m_vecVarIndex[c] =
				m_vecVarsAtLocation[(size_t)DataLocation_REdge]++;

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

void Grid::EvaluateGeometricTerms() {
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		m_vecActiveGridPatches[n]->EvaluateGeometricTerms();
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

void Grid::Checksum(
	DataType eDataType,
	DataArray1D<double> & dChecksums,
	int iDataIndex,
	ChecksumType eChecksumType
) const {
#ifdef USE_MPI
	// Identify root process
	int nRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

	// Initialize the local checksum array from DataType
	DataArray1D<double> dChecksumsLocal;
	if (eDataType == DataType_State) {
		dChecksumsLocal.Allocate(m_model.GetEquationSet().GetComponents());

	} else if (eDataType == DataType_Tracers) { 
		int nTracers = m_model.GetEquationSet().GetTracers();
		if (nTracers == 0) {
			return;
		}

		dChecksumsLocal.Allocate(nTracers);

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
		dChecksums.Allocate(dChecksumsLocal.GetRows());
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
#endif
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

	// Global energy
	double dGlobalEnergy = 0.0;

#ifdef USE_MPI
	// Reduce to obtain global energy integral
	MPI_Reduce(
		&dLocalEnergy,
		&dGlobalEnergy,
		1,
		MPI_DOUBLE,
		MPI_SUM,
		0,
		MPI_COMM_WORLD);
#endif

	// Return global energy integral
	return dGlobalEnergy;
}

///////////////////////////////////////////////////////////////////////////////

double Grid::ComputeTotalPotentialEnstrophy(
	int iDataIndex
) {
	// Compute vorticity and divergence on the Grid
	ComputeVorticityDivergence(iDataIndex);

	// Compute local potential enstrophy
	double dLocalPotentialEnstrophy = 0.0;
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		dLocalPotentialEnstrophy +=
			m_vecActiveGridPatches[n]->
				ComputeTotalPotentialEnstrophy(iDataIndex);
	}

	// Global potential enstrophy
	double dGlobalPotentialEnstrophy = 0.0;

#ifdef USE_MPI
	// Reduce to obtain global energy integral
	MPI_Reduce(
		&dLocalPotentialEnstrophy,
		&dGlobalPotentialEnstrophy,
		1,
		MPI_DOUBLE,
		MPI_SUM,
		0,
		MPI_COMM_WORLD);
#endif

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

#ifdef USE_MPI
	// Verify all processors are prepared to exchange
	MPI_Barrier(MPI_COMM_WORLD);
#endif

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

#ifdef USE_MPI
	// Verify all processors are prepared to exchange
	MPI_Barrier(MPI_COMM_WORLD);
#endif

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

#ifdef USE_MPI
	// Verify all processors are prepared to exchange
	MPI_Barrier(MPI_COMM_WORLD);
#endif

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

int Grid::GetMaxNodeCount2D() const {

	// Total number of nodes over all patches of grid
	int nMaxNodes2D = 0;

	// Loop over all patches and obtain total node count
	for (int n = 0; n < m_aPatchBoxes.GetRows(); n++) {
		int nPatchNodes2D = m_aPatchBoxes[n].GetTotalNodeCount2D();

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
	for (int n = 0; n < m_aPatchBoxes.GetRows(); n++) {
		nTotalNodes2D += m_aPatchBoxes[n].GetTotalNodeCount2D();
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
	for (int n = 0; n < m_aPatchBoxes.GetRows(); n++) {
		const PatchBox & box = GetPatchBox(n);

		int nPatchNodes;
		if (loc == DataLocation_Node) {
			nPatchNodes = box.GetTotalNodes() * GetRElements();
		} else if (loc == DataLocation_REdge) {
			nPatchNodes = box.GetTotalNodes() * (GetRElements() + 1);
		} else {
			_EXCEPTIONT("Invalid location");
		}

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
	for (int n = 0; n < m_aPatchBoxes.GetRows(); n++) {
		const PatchBox & box = GetPatchBox(n);

		int nPatchNodes;
		if (loc == DataLocation_Node) {
			nPatchNodes = box.GetTotalNodes() * GetRElements();
		} else if (loc == DataLocation_REdge) {
			nPatchNodes = box.GetTotalNodes() * (GetRElements() + 1);
		} else {
			_EXCEPTIONT("Invalid location");
		}

		nTotalNodes += nPatchNodes;
	}

	return nTotalNodes;
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
	const DataArray1D<double> & dAlpha,
	const DataArray1D<double> & dBeta,
	const DataArray1D<int> & iPatch,
	DataType eDataType,
	DataLocation eDataLocation,
	bool fInterpAllVariables,
	DataArray3D<double> & dInterpData,
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

#ifdef USE_MPI
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
#endif
}

///////////////////////////////////////////////////////////////////////////////

void Grid::ConvertReferenceToPatchCoord(
	const DataArray1D<double> & dXReference,
	const DataArray1D<double> & dYReference,
	DataArray1D<double> & dAlpha,
	DataArray1D<double> & dBeta,
	DataArray1D<int> & iPatch
) const {
	_EXCEPTIONT("Unimplemented.");
}

///////////////////////////////////////////////////////////////////////////////

GridPatch * Grid::ActivateEmptyPatch(
	int ixPatch
) {
	GridPatch * pPatch = NewPatch(ixPatch);

	m_vecActiveGridPatches.push_back(pPatch);
	m_vecActiveGridPatchIndices.push_back(ixPatch);

	return pPatch;
}

///////////////////////////////////////////////////////////////////////////////

void Grid::DeactivatePatch(
	int ixPatch
) {
	for (int i = 0; i < m_vecActiveGridPatchIndices.size(); i++) {
		if (m_vecActiveGridPatchIndices[i] == ixPatch) {
			delete m_vecActiveGridPatches[i];

			m_vecActiveGridPatchIndices.erase(
				m_vecActiveGridPatchIndices.begin() + i);
			m_vecActiveGridPatches.erase(
				m_vecActiveGridPatches.begin() + i);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void Grid::DistributePatches() {
#ifdef USE_MPI
	// Number of processors
	int nSize;
	MPI_Comm_size(MPI_COMM_WORLD, &nSize);

	// Current processor
	int nRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

	// Loop over all patches and initialize data
	m_vecPatchProcessor.resize(m_aPatchBoxes.GetRows());
	for (int n = 0; n < m_aPatchBoxes.GetRows(); n++) {
		int iPatchProcessor = n % nSize;
		m_vecPatchProcessor[n] = iPatchProcessor;

		if (iPatchProcessor == nRank) {
			GridPatch * pPatch = NewPatch(n);
			pPatch->InitializeDataLocal();
			m_vecActiveGridPatches.push_back(pPatch);
			m_vecActiveGridPatchIndices.push_back(n);
		}
	}
#endif
}

///////////////////////////////////////////////////////////////////////////////

void Grid::RegisterExchangeBuffer(
	int ixSourcePatch,
	int ixTargetPatch,
	Direction dir
) {
	const PatchBox & box = GetPatchBox(ixSourcePatch);

	// First patch coordinate index
	int ixFirst;
	int ixSecond;

	if ((dir == Direction_Right) ||
		(dir == Direction_Left)
	) {
		ixFirst  = box.GetBInteriorBegin();
		ixSecond = box.GetBInteriorEnd();

	} else if (
		(dir == Direction_Top) ||
		(dir == Direction_Bottom)
	) {
		ixFirst  = box.GetAInteriorBegin();
		ixSecond = box.GetAInteriorEnd();

	} else if (dir == Direction_TopRight) {
		ixFirst  = box.GetAInteriorEnd()-1;
		ixSecond = box.GetBInteriorEnd()-1;

	} else if (dir == Direction_TopLeft) {
		ixFirst  = box.GetAInteriorBegin();
		ixSecond = box.GetBInteriorEnd()-1;

	} else if (dir == Direction_BottomLeft) {
		ixFirst  = box.GetAInteriorBegin();
		ixSecond = box.GetBInteriorBegin();

	} else if (dir == Direction_BottomRight) {
		ixFirst  = box.GetAInteriorEnd()-1;
		ixSecond = box.GetBInteriorBegin();

	} else {
		_EXCEPTIONT("Invalid direction");
	}

	// Exterior connect
	RegisterExchangeBuffer(
		ixSourcePatch,
		ixTargetPatch,
		dir,
		ixFirst,
		ixSecond);
}

///////////////////////////////////////////////////////////////////////////////

void Grid::RegisterExchangeBuffer(
	int ixSourcePatch,
	int ixTargetPatch,
	Direction dir,
	int ixFirst,
	int ixSecond
) {
	if (ixSourcePatch == GridPatch::InvalidIndex) {
		_EXCEPTIONT("Invalid ixSourcePatch");
	}

	// Check for NULL patches (do nothing)
	if (ixTargetPatch == GridPatch::InvalidIndex) {
		return;
	}

	// Get the GridPatch's PatchBox
	const PatchBox & boxSource = GetPatchBox(ixSourcePatch);
	const PatchBox & boxTarget = GetPatchBox(ixTargetPatch);

	// Build key
	ExchangeBufferInfo info;
	info.ixSourcePatch = ixSourcePatch;
	info.ixTargetPatch = ixTargetPatch;
	info.dir = dir;

	// Build exhange buffer metadata
	info.ixFirst = ixFirst;
	info.ixSecond = ixSecond;

	// Get number of components
	const Model & model = GetModel();

	const EquationSet & eqn = model.GetEquationSet();

	size_t sStateTracerMaxVariables;
	if (eqn.GetComponents() > eqn.GetTracers()) {
		sStateTracerMaxVariables = eqn.GetComponents();
	} else {
		sStateTracerMaxVariables = eqn.GetTracers();
	}

	info.sHaloElements = model.GetHaloElements();
	info.sComponents = sStateTracerMaxVariables;
	info.sMaxRElements = GetRElements() + 1;

	// Get the opposing direction
	GetOpposingDirection(
		boxSource.GetPanel(),
		boxTarget.GetPanel(),
		dir,
		info.dirOpposing,
		info.fReverseDirection,
		info.fFlippedCoordinate);

	// Determine the size of the boundary (number of elements along exterior
	// edge).  Used in computing the size of the send/recv buffers.
	if ((dir == Direction_Right) ||
		(dir == Direction_Top) ||
		(dir == Direction_Left) ||
		(dir == Direction_Bottom)
	) {
		info.sBoundarySize = ixSecond - ixFirst;
	} else {
		info.sBoundarySize = info.sHaloElements;
	}

	info.CalculateByteSize();

	if ((dir == Direction_TopRight) && (
		(ixFirst < boxSource.GetAInteriorBegin() + info.sBoundarySize - 1) ||
		(ixSecond < boxSource.GetBInteriorBegin() + info.sBoundarySize - 1)
	)) {
		_EXCEPTIONT("Insufficient interior elements to build "
			"diagonal connection.");
	}

	if ((dir == Direction_TopLeft) && (
		(ixFirst > boxSource.GetAInteriorEnd() - info.sBoundarySize) ||
		(ixSecond < boxSource.GetBInteriorBegin() + info.sBoundarySize - 1)
	)) {
		_EXCEPTIONT("Insufficient interior elements to build "
			"diagonal connection.");
	}

	if ((dir == Direction_BottomLeft) && (
		(ixFirst > boxSource.GetAInteriorEnd() - info.sBoundarySize) ||
		(ixSecond > boxSource.GetBInteriorEnd() - info.sBoundarySize)
	)) {
		_EXCEPTIONT("Insufficient interior elements to build "
			"diagonal connection.");
	}

	if ((dir == Direction_BottomRight) && (
		(ixFirst < boxSource.GetAInteriorBegin() + info.sBoundarySize - 1) ||
		(ixSecond > boxSource.GetBInteriorEnd() - info.sBoundarySize)
	)) {
		_EXCEPTIONT("Insufficient interior elements to build "
			"diagonal connection.");
	}

	// Add the exchange buffer information to the registry
	m_aExchangeBufferRegistry.Register(info);
}

///////////////////////////////////////////////////////////////////////////////

void Grid::InitializeExchangeBuffersFromPatch(
	int ixSourcePatch
) {
	const PatchBox & box = GetPatchBox(ixSourcePatch);

	// Vector of nodal points around element
	int nPerimeter = box.GetInteriorPerimeter() + 4;

	DataArray1D<int> vecIxA(nPerimeter);
	DataArray1D<int> vecIxB(nPerimeter);
	DataArray1D<int> vecPanel(nPerimeter);
	DataArray1D<int> vecPatch(nPerimeter);

	// Perimeter node index
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
		RegisterExchangeBuffer(
			ixSourcePatch,
			vecPatch[ix],
			Direction_BottomLeft);
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
				RegisterExchangeBuffer(
					ixSourcePatch,
					iCurrentPatch,
					Direction_Bottom,
					ixFirstBegin,
					i);

				if (i != box.GetAInteriorEnd()) {
					RegisterExchangeBuffer(
						ixSourcePatch,
						iCurrentPatch,
						Direction_BottomLeft,
						i,
						box.GetBInteriorBegin());

					ixFirstBegin = i;
					iCurrentPatch = vecPatch[ix];

					RegisterExchangeBuffer(
						ixSourcePatch,
						iCurrentPatch,
						Direction_BottomRight,
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
		RegisterExchangeBuffer(
			ixSourcePatch,
			vecPatch[ix],
			Direction_BottomRight);
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
				RegisterExchangeBuffer(
					ixSourcePatch,
					iCurrentPatch,
					Direction_Right,
					ixFirstBegin,
					j);

				if (j != box.GetBInteriorEnd()) {
					RegisterExchangeBuffer(
						ixSourcePatch,
						iCurrentPatch,
						Direction_BottomRight,
						box.GetAInteriorEnd()-1,
						j);

					ixFirstBegin = j;
					iCurrentPatch = vecPatch[ix];

					RegisterExchangeBuffer(
						ixSourcePatch,
						iCurrentPatch,
						Direction_TopRight,
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
		RegisterExchangeBuffer(
			ixSourcePatch,
			vecPatch[ix],
			Direction_TopRight);
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
				RegisterExchangeBuffer(
					ixSourcePatch,
					iCurrentPatch,
					Direction_Top,
					i + 1,
					ixFirstEnd);

				if (i != box.GetAInteriorBegin()-1) {
					RegisterExchangeBuffer(
						ixSourcePatch,
						iCurrentPatch,
						Direction_TopRight,
						i,
						box.GetBInteriorEnd()-1);

					ixFirstEnd = i + 1;
					iCurrentPatch = vecPatch[ix];

					RegisterExchangeBuffer(
						ixSourcePatch,
						iCurrentPatch,
						Direction_TopLeft,
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
		RegisterExchangeBuffer(
			ixSourcePatch,
			vecPatch[ix],
			Direction_TopLeft);
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
				RegisterExchangeBuffer(
					ixSourcePatch,
					iCurrentPatch,
					Direction_Left,
					j + 1,
					ixFirstEnd);

				if (j != box.GetBInteriorBegin()-1) {
					RegisterExchangeBuffer(
						ixSourcePatch,
						iCurrentPatch,
						Direction_TopLeft,
						box.GetAInteriorBegin(),
						j);

					ixFirstEnd = j + 1;
					iCurrentPatch = vecPatch[ix];

					RegisterExchangeBuffer(
						ixSourcePatch,
						iCurrentPatch,
						Direction_BottomLeft,
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
}

///////////////////////////////////////////////////////////////////////////////

void Grid::InitializeExchangeBuffersFromActivePatches() {
	for (int n = 0; n < m_vecActiveGridPatchIndices.size(); n++) {
		InitializeExchangeBuffersFromPatch(m_vecActiveGridPatchIndices[n]);
	}
}

///////////////////////////////////////////////////////////////////////////////

void Grid::InitializeAllExchangeBuffers() {
	for (int n = 0; n < m_nInitializedPatchBoxes; n++) {
		InitializeExchangeBuffersFromPatch(n);
	}
}

///////////////////////////////////////////////////////////////////////////////

void Grid::InitializeConnectivity(
	bool fAllocate
) {

	// Clear all existing neighbors fro all active patches
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		m_vecActiveGridPatches[n]->m_connect.ClearNeighbors();
	}

	// Loop through all exchange buffers
	const ExchangeBufferRegistry::ExchangeBufferInfoVector & vecRegistry =
		m_aExchangeBufferRegistry.GetRegistry();

	for (int m = 0; m < vecRegistry.size(); m++) {

		const ExchangeBufferInfo & info = vecRegistry[m];

		// If the GridPatch associated with this exchange buffer is active,
		// add neighbors
		GridPatch * pPatch = NULL;
		for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
			if (m_vecActiveGridPatchIndices[n] == info.ixSourcePatch) {
				pPatch = m_vecActiveGridPatches[n];
				break;
			}
		}

		if (pPatch != NULL) {

			const PatchBox & boxSource = GetPatchBox(info.ixSourcePatch);
			const PatchBox & boxTarget = GetPatchBox(info.ixTargetPatch);

			// Build the new ExteriorNeighbor
			ExteriorNeighbor * pNeighbor =
				new ExteriorNeighbor(info);

			if (fAllocate) {
				m_aExchangeBufferRegistry.Allocate(m);
			}

			pNeighbor->AttachBuffers(
				m_aExchangeBufferRegistry.GetRecvBufferPtr(m),
				m_aExchangeBufferRegistry.GetSendBufferPtr(m));

			pPatch->m_connect.AddExteriorNeighbor(pNeighbor);
		}
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
	const DataArray1D<double> & dCoeff,
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

