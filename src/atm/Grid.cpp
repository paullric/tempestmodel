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

#include "Exception.h"

#include <cmath>

#include "mpi.h"

///////////////////////////////////////////////////////////////////////////////

Grid::Grid(
	const Model & model,
	int nABaseResolution,
	int nBBaseResolution,
	int nRefinementRatio,
	int nRElements
) :
	m_fInitialized(false),
	m_model(model),
	m_nABaseResolution(nABaseResolution),
	m_nBBaseResolution(nBBaseResolution),
	m_nRefinementRatio(nRefinementRatio),
	m_nRElements(nRElements),
	m_dZtop(1.0)
{
}

///////////////////////////////////////////////////////////////////////////////

Grid::~Grid() {
	for (int n = 0; n < m_vecGridPatches.size(); n++) {
		delete m_vecGridPatches[n];
	}
}

///////////////////////////////////////////////////////////////////////////////

void Grid::InitializeVerticalCoordinate(
	const GridSpacing & gridspacing
) {
	// Initialize location and index for each variable
	m_vecVarLocation.resize(m_model.GetEquationSet().GetComponents());

	// If dimensionality is 2, then initialize a dummy REta
	if (m_model.GetEquationSet().GetDimensionality() == 2) {

		// Resize arrays of REta coordinates
		m_dREtaLevels.Initialize(1);
		m_dREtaInterfaces.Initialize(2);

		// Uniform grid spacing
		m_dREtaInterfaces[0] = 0.0;
		m_dREtaInterfaces[1] = 1.0;
		m_dREtaLevels[0] = 0.5;

#pragma message "Do not hardcode this information"
		// Everything on model levels
		m_vecVarLocation[0] = DataLocation_Node;
		m_vecVarLocation[1] = DataLocation_Node;
		m_vecVarLocation[2] = DataLocation_Node;

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

		// Uniform grid spacing
		for (int k = 0; k <= m_nRElements; k++) {
			m_dREtaInterfaces[k] = gridspacing.GetEdge(k);
		}
		for (int k = 0; k < m_nRElements; k++) {
			m_dREtaLevels[k] = gridspacing.GetNode(k);
		}

		// Location of variables
		m_vecVarLocation[0] = DataLocation_Node;
		m_vecVarLocation[1] = DataLocation_Node;
		m_vecVarLocation[2] = DataLocation_Node;
		m_vecVarLocation[3] = DataLocation_Node;
		m_vecVarLocation[4] = DataLocation_Node;

	} else {
		_EXCEPTIONT("Invalid dimensionality");
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
	double dTime,
	int iDataIndex
) {
	// Store the model cap
	m_dZtop = test.GetZtop();

	// Store the reference state flag
	m_fHasReferenceState = test.HasReferenceState();

	// Evaluate the pointwise values of the test
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		m_vecActiveGridPatches[n]->
			EvaluateTestCase(test, dTime, iDataIndex);
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

void Grid::Exchange(
	DataType eDataType,
	int iDataIndex
) {
	// Set up asynchronous recvs
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		m_vecActiveGridPatches[n]->PrepareExchange();
	}

	// Send data
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		m_vecActiveGridPatches[n]->Send(eDataType, iDataIndex);
	}
/*
	// Postprocess data
	int nExpectedMessageCount = 0;

	for (int n = 0; n < m_vecGridPatches.size(); n++) {
		const GridPatch * pPatch = m_vecGridPatches[n];
		if (pPatch->ContainsData()) {
			nExpectedMessageCount +=
				pPatch->GetConnectivity().GetExpectedMessageCount();
		}
	}

	// Receive data from exterior neighbors
	int nRecvMessageCount = 0;
	for (; nRecvMessageCount != nExpectedMessageCount; nRecvMessageCount++) {
	}
*/
	// Receive data
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		m_vecActiveGridPatches[n]->Receive(eDataType, iDataIndex);
	}
}

///////////////////////////////////////////////////////////////////////////////

int Grid::GetLargestGridPatchNodes() const {

	// Most nodes per patch
	int nMaxNodes = 0;

	// Loop over all patches
	for (int n = 0; n < m_vecGridPatches.size(); n++) {
		int nPatchNodes =
			m_vecGridPatches[n]->GetPatchBox().GetTotalNodes();

		if (nPatchNodes > nMaxNodes) {
			nMaxNodes = nPatchNodes;
		}
	}

	return nMaxNodes;
}

///////////////////////////////////////////////////////////////////////////////

int Grid::GetTotalNodeCount() const {

	// Total number of nodes over all patches of grid
	int nTotalNodes = 0;

	// Loop over all patches and obtain total node count
	for (int n = 0; n < m_vecGridPatches.size(); n++) {
		nTotalNodes += m_vecGridPatches[n]->GetTotalNodeCount();
	}

	return nTotalNodes;
}

///////////////////////////////////////////////////////////////////////////////

int Grid::GetMaximumDegreesOfFreedom() const {

	// Most nodes per patch
	int nMaxDOFs = 0;

	// Loop over all patches and obtain max DOFs from state and tracer data
	for (int n = 0; n < m_vecGridPatches.size(); n++) {
		int nStateDOFs =
			m_model.GetEquationSet().GetComponents()
			* m_nRElements
			* m_vecGridPatches[n]->GetPatchBox().GetTotalNodes();

		int nTracersDOFs =
			m_model.GetEquationSet().GetTracers()
			* m_nRElements
			* m_vecGridPatches[n]->GetPatchBox().GetTotalNodes();

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
	DataType & eRecvDataType
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
		mpistatus.MPI_TAG, ixRecvPatch, eRecvDataType);

	if ((ixRecvPatch < 0) || (ixRecvPatch >= m_vecGridPatches.size())) {
		_EXCEPTIONT("Panel tag index out of range");
	}

	status.SetReceiveStatus(ixRecvPatch, eRecvDataType);

	// Verify consistency of patch information
	GridPatch * pPatch = m_vecGridPatches[ixRecvPatch];

	MPI_Get_count(&mpistatus, MPI_DOUBLE, &nRecvCount);

	if (eRecvDataType == DataType_State) {
		int nExpectedRecvCount =
			m_model.GetEquationSet().GetComponents()
			* m_nRElements
			* pPatch->GetPatchBox().GetTotalNodes();

		if (nExpectedRecvCount != nRecvCount) {
			_EXCEPTIONT("State dimension mismatch");
		}

	} else if (eRecvDataType == DataType_Tracers) {
		int nExpectedRecvCount =
			m_model.GetEquationSet().GetTracers()
			* m_nRElements
			* pPatch->GetPatchBox().GetTotalNodes();

		if (nExpectedRecvCount != nRecvCount) {
			_EXCEPTIONT("Tracers dimension mismatch");
		}

	} else if (eRecvDataType == DataType_Jacobian) {
		int nExpectedRecvCount =
			m_nRElements
			* pPatch->GetPatchBox().GetTotalNodes();

		int nDiff =
			GetCumulativePatch3DNodeIndex(ixRecvPatch+1)
			- GetCumulativePatch3DNodeIndex(ixRecvPatch);

		if (nExpectedRecvCount != nRecvCount) {
			_EXCEPTIONT("Jacobian dimension mismatch");
		}
	}
#pragma message "Perform check for other data types"
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
		const GridData4D & dataState    = pPatch->GetDataState(0);
		const GridData4D & dataTracers  = pPatch->GetDataTracers(0);

		const DataMatrix3D<double> & dataJacobian   = pPatch->GetJacobian();
		const DataMatrix<double>   & dataTopography = pPatch->GetTopography();
		const DataMatrix<double>   & dataLongitude  = pPatch->GetLongitude();
		const DataMatrix<double>   & dataLatitude   = pPatch->GetLatitude();
		const DataMatrix3D<double> & dataZLevels    = pPatch->GetZLevels();

		// Send state data to root process
		if (status.Contains(DataType_State)) {
			MPI_Isend(
				(void*)(dataState[0][0][0]),
				dataState.GetTotalElements(),
				MPI_DOUBLE,
				0,
				ConsolidationStatus::GenerateTag(
					pPatch->GetPatchIndex(), DataType_State),
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
	const DataVector<int> & iPanel,
	DataType eDataType,
	DataMatrix3D<double> & dInterpData,
	bool fIncludeReferenceState,
	bool fConvertToPrimitive
) const {
	if ((dAlpha.GetRows() != dBeta.GetRows()) ||
		(dAlpha.GetRows() != iPanel.GetRows())
	) {
		_EXCEPTIONT("Inconsistency in vector lengths.");
	}
	if ((eDataType == DataType_Tracers) &&
		(m_model.GetEquationSet().GetTracers() == 0)
	) {
		_EXCEPTIONT("Unable to Interpolate with no tracers.");
	}

	// Initialize interpolated data
	if (eDataType == DataType_State) {
		dInterpData.Initialize(
			m_model.GetEquationSet().GetComponents(),
			GetRElements(),
			dAlpha.GetRows());

	} else if (eDataType == DataType_Tracers) {
		dInterpData.Initialize(
			m_model.GetEquationSet().GetTracers(),
			GetRElements(),
			dAlpha.GetRows());

	} else if (eDataType == DataType_Vorticity) {
		dInterpData.Initialize(
			1,
			GetRElements(),
			dAlpha.GetRows());

	} else if (eDataType == DataType_Divergence) {
		dInterpData.Initialize(
			1,
			GetRElements(),
			dAlpha.GetRows());

	} else {
		_EXCEPTIONT("Invalid DataType / Not implemented.");
	}

	// Interpolate state data
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		m_vecActiveGridPatches[n]->InterpolateData(
			dAlpha, dBeta, iPanel,
			eDataType,
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
		+ pPatch->GetTotalNodeCount());

	return pPatch;
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

void Grid::AddReferenceState(
	int ix
) {
	// Loop over all grid patches
	for (int n = 0; n < m_vecActiveGridPatches.size(); n++) {
		m_vecActiveGridPatches[n]->AddReferenceState(ix);
	}
}

///////////////////////////////////////////////////////////////////////////////

void Grid::ConvertReferenceToABP(
	const DataVector<double> & dXReference,
	const DataVector<double> & dYReference,
	DataVector<double> & dAlpha,
	DataVector<double> & dBeta,
	DataVector<int> & iPanel
) const {
	_EXCEPTIONT("Unimplemented.");
}

///////////////////////////////////////////////////////////////////////////////

