///////////////////////////////////////////////////////////////////////////////
///
///	\file    OutputManagerComposite.cpp
///	\author  Paul Ullrich
///	\version August 11, 2010
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

#include "OutputManagerComposite.h"

#include "Model.h"
#include "Grid.h"
#include "ConsolidationStatus.h"

#include "TimeObj.h"
#include "Announce.h"

#ifdef USE_MPI
#include "mpi.h"
#endif

#include <iostream>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <sys/stat.h>

///////////////////////////////////////////////////////////////////////////////

OutputManagerComposite::OutputManagerComposite(
	Grid & grid,
	const Time & timeOutputFrequency,
	std::string strOutputDir,
	std::string strOutputFormat,
	std::string strRestartFile
) :
	OutputManager(
		grid,
		timeOutputFrequency,
		strOutputDir,
		strOutputFormat,
		1),
	m_strRestartFile(strRestartFile)
{
}

///////////////////////////////////////////////////////////////////////////////

OutputManagerComposite::~OutputManagerComposite() {
	CloseFile();
}

///////////////////////////////////////////////////////////////////////////////

bool OutputManagerComposite::OpenFile(
	const std::string & strFileName
) {
#ifdef USE_MPI
	// Determine processor rank
	int nRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

	// Determine space allocation for each GridPatch
	m_vecGridPatchByteSize.Allocate(m_grid.GetPatchCount(), 2);

	// Determine maximum recv buffer size
	for (int i = 0; i < m_grid.GetActivePatchCount(); i++) {
		const GridPatch * pPatch = m_grid.GetActivePatch(i);

		int iPatchIx = pPatch->GetPatchIndex();
		if (iPatchIx > m_grid.GetPatchCount()) {
			_EXCEPTION2("PatchIndex (%i) out of range [0,%i)",
				iPatchIx, m_grid.GetPatchCount());
		}

		const DataContainer & dcGeometric =
			pPatch->GetDataContainerGeometric();
		m_vecGridPatchByteSize[iPatchIx][0] =
			dcGeometric.GetTotalByteSize();

		const DataContainer & dcActiveState =
			pPatch->GetDataContainerActiveState();
		m_vecGridPatchByteSize[iPatchIx][1] =
			dcActiveState.GetTotalByteSize();
	}

	// Reduce GridPatch byte size
	if (nRank != 0) {
		MPI_Reduce(
			&(m_vecGridPatchByteSize[0][0]),
			&(m_vecGridPatchByteSize[0][0]),
			m_vecGridPatchByteSize.GetTotalSize(),
			MPI_INT,
			MPI_MAX,
			0,
			MPI_COMM_WORLD);

	// Open file at root
	} else if (nRank == 0) {

		// Reduce
		MPI_Reduce(
			MPI_IN_PLACE,
			&(m_vecGridPatchByteSize[0][0]),
			m_vecGridPatchByteSize.GetTotalSize(),
			MPI_INT,
			MPI_MAX,
			0,
			MPI_COMM_WORLD);

		// The active Model
		const Model & model = m_grid.GetModel();

		// Get maximum byte size for exchange
		int sMaxRecvBufferByteSize = 0;
		for (int i = 0; i < m_vecGridPatchByteSize.GetRows(); i++) {
			if (m_vecGridPatchByteSize[i][0] > sMaxRecvBufferByteSize) {
				sMaxRecvBufferByteSize = m_vecGridPatchByteSize[i][0];
			}
			if (m_vecGridPatchByteSize[i][1] > sMaxRecvBufferByteSize) {
				sMaxRecvBufferByteSize = m_vecGridPatchByteSize[i][1];
			}
		}

		// Allocate Recv buffer at root
		if (m_vecRecvBuffer.GetRows() < sMaxRecvBufferByteSize) {
			m_vecRecvBuffer.Allocate(sMaxRecvBufferByteSize);
		}

		// Initialize byte location for each GridPatch
		m_vecGridPatchByteLoc.Allocate(m_grid.GetPatchCount(), 2);
		m_vecGridPatchByteLoc[0][0] = 0;
		m_vecGridPatchByteLoc[0][1] = m_vecGridPatchByteSize[0][0];
		for (int i = 1; i < m_vecGridPatchByteSize.GetRows(); i++) { 
			m_vecGridPatchByteLoc[i][0] =
				m_vecGridPatchByteLoc[i-1][1]
				+ m_vecGridPatchByteSize[i-1][1];

			m_vecGridPatchByteLoc[i][1] =
				m_vecGridPatchByteLoc[i][0]
				+ m_vecGridPatchByteSize[i][0];
		}

		// Allocate array of received messages
		m_iReceivedMessages.Allocate(m_grid.GetPatchCount(), 2);

		// Check for existing file
		if (m_ofsActiveOutput.is_open()) {
			_EXCEPTIONT("Restart file already open");
		}

		// Open new binary output stream
		std::string strRestartFileName = strFileName + ".restart.dat";
		m_ofsActiveOutput.open(
			strRestartFileName.c_str(), std::ios::binary | std::ios::out);

		if (!m_ofsActiveOutput) {
			_EXCEPTION1("Error opening output file \"%s\"",
				strRestartFileName.c_str());
		}
/*
		for (int i = 0; i < m_vecGridPatchByteSize.GetRows(); i++) {
			printf("%i %i\n",
				m_vecGridPatchByteLoc[i][0],
				m_vecGridPatchByteLoc[i][1]);
		}
*/
	}

	// Wait for all processes to complete
	MPI_Barrier(MPI_COMM_WORLD);

	return true;
#else
	_EXCEPTIONT("Not implemented without USE_MPI");
#endif

}	

///////////////////////////////////////////////////////////////////////////////

void OutputManagerComposite::CloseFile() {
	m_ofsActiveOutput.close();
}

///////////////////////////////////////////////////////////////////////////////

void OutputManagerComposite::Output(
	const Time & time
) {
#ifdef USE_MPI
	// Check for open file
	if (!IsFileOpen()) {
		_EXCEPTIONT("No file available for output");
	}

	// Verify that only one output has been performed
	if (m_ixOutputTime != 0) {
		_EXCEPTIONT("Only one Composite output allowed per file");
	}

	// Determine processor rank
	int nRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

	// Data types
	const int ExchangeDataType_Geometric = 0;
	const int ExchangeDataType_ActiveState = 1;
	const int ExchangeDataType_Count = 2;

	// Send data from GridPatches to root
	if (nRank != 0) {

		int nActivePatches = m_grid.GetActivePatchCount();

		DataArray1D<MPI_Request> vecSendReqGeo(nActivePatches);
		DataArray1D<MPI_Request> vecSendReqAcS(nActivePatches);

		for (int i = 0; i < nActivePatches; i++) {
			const GridPatch * pPatch = m_grid.GetActivePatch(i);

			const DataContainer & dcGeometric =
				pPatch->GetDataContainerGeometric();
			int nGeometricDataByteSize =
				dcGeometric.GetTotalByteSize();
			const unsigned char * pGeometricData =
				dcGeometric.GetPointer();

			MPI_Isend(
				pGeometricData,
				nGeometricDataByteSize,
				MPI_BYTE,
				0,
				ExchangeDataType_Geometric,
				MPI_COMM_WORLD,
				&(vecSendReqGeo[i]));

			const DataContainer & dcActiveState =
				pPatch->GetDataContainerActiveState();
			int nActiveStateDataByteSize =
				dcActiveState.GetTotalByteSize();
			const unsigned char * pActiveStateData =
				dcActiveState.GetPointer();

			MPI_Isend(
				pActiveStateData,
				nActiveStateDataByteSize,
				MPI_BYTE,
				0,
				ExchangeDataType_ActiveState,
				MPI_COMM_WORLD,
				&(vecSendReqAcS[i]));
/*
			printf("Send Geo %i %i\n", pPatch->GetPatchIndex(), nGeometricDataByteSize);
			printf("Send AcS %i %i\n", pPatch->GetPatchIndex(), nActiveStateDataByteSize);
*/
		}

		int iWaitAllMsgGeo =
			MPI_Waitall(
				nActivePatches,
				&(vecSendReqGeo[0]),
				MPI_STATUSES_IGNORE);

		if (iWaitAllMsgGeo == MPI_ERR_IN_STATUS) {
			_EXCEPTIONT("MPI_Waitall returned MPI_ERR_IN_STATUS");
		}

		int iWaitAllMsgAcS =
			MPI_Waitall(
				nActivePatches,
				&(vecSendReqAcS[0]),
				MPI_STATUSES_IGNORE);

		if (iWaitAllMsgAcS == MPI_ERR_IN_STATUS) {
			_EXCEPTIONT("MPI_Waitall returned MPI_ERR_IN_STATUS");
		}

		//std::cout << "WAIT DONE" << std::endl;

	// Receive data and output to file 
	} else if (nRank == 0) {

		// Active Model
		const Model & model = m_grid.GetModel();

		// Output start time as string
		std::string strStartTime = model.GetStartTime().ToLongString();
		int nStartTimeLen = strStartTime.length();
		m_ofsActiveOutput.write(
			(const char *)(&nStartTimeLen), sizeof(int));
		m_ofsActiveOutput.write(
			strStartTime.c_str(), sizeof(char)*nStartTimeLen);

		// Output equation set name
		std::string strEquationSetName = model.GetEquationSet().GetName();
		int nEquationSetNameLen = strEquationSetName.length();
		m_ofsActiveOutput.write(
			(const char *)(&nEquationSetNameLen), sizeof(int));
		m_ofsActiveOutput.write(
			strEquationSetName.c_str(), sizeof(char)*nEquationSetNameLen);

		// Output the grid
		//m_grid.ToFile(*m_pActiveNcOutput);

		// Reference position
		std::streampos posRefFile = m_ofsActiveOutput.tellp();
		if (posRefFile == (-1)) {
			_EXCEPTIONT("ActiveOutput::tellp() fail");
		}

		// Write my GridPatch data to file
		for (int i = 0; i < m_grid.GetActivePatchCount(); i++) {
			const GridPatch * pPatch = m_grid.GetActivePatch(i);

			int iPatchIx = pPatch->GetPatchIndex();
			if (iPatchIx > m_grid.GetPatchCount()) {
				_EXCEPTION2("PatchIndex (%i) out of range [0,%i)",
					iPatchIx, m_grid.GetPatchCount());
			}
/*
			printf("Write %i %i %i %i %i\n",
				iPatchIx,
				m_vecGridPatchByteLoc[iPatchIx][0],
				m_vecGridPatchByteSize[iPatchIx][0],
				m_vecGridPatchByteLoc[iPatchIx][1],
				m_vecGridPatchByteSize[iPatchIx][1]);
*/
			const DataContainer & dcGeometric =
				pPatch->GetDataContainerGeometric();
			int nGeometricDataByteSize =
				dcGeometric.GetTotalByteSize();
			const char * pGeometricData =
				(const char *)(dcGeometric.GetPointer());

			m_ofsActiveOutput.seekp(
				posRefFile + m_vecGridPatchByteLoc[iPatchIx][0]);
			m_ofsActiveOutput.write(
				pGeometricData, m_vecGridPatchByteSize[iPatchIx][0]);

			const DataContainer & dcActiveState =
				pPatch->GetDataContainerActiveState();
			int nActiveStateDataByteSize =
				dcActiveState.GetTotalByteSize();
			const char * pActiveStateData =
				(const char *)(dcActiveState.GetPointer());

			m_ofsActiveOutput.seekp(
				posRefFile + m_vecGridPatchByteLoc[iPatchIx][1]);
			m_ofsActiveOutput.write(
				pActiveStateData, m_vecGridPatchByteSize[iPatchIx][1]);
		}

		// Recieve all data objects from neighbors
		int nRemainingMessages =
			2 * (m_grid.GetPatchCount() - m_grid.GetActivePatchCount());

		for (; nRemainingMessages > 0; nRemainingMessages--) {

			// Receive a consolidation message
			MPI_Status status;

			MPI_Recv(
				&(m_vecRecvBuffer[0]),
				m_vecRecvBuffer.GetRows(),
				MPI_CHAR,
				MPI_ANY_SOURCE,
				MPI_ANY_TAG,
				MPI_COMM_WORLD,
				&status);

			// Data type from TAG
			int iDataType = status.MPI_TAG;
			if ((iDataType < 0) || (iDataType >= ExchangeDataType_Count)) {
				_EXCEPTION1("MPI_TAG (%i) out of range", iDataType);
			}

			// Check patch index
			DataArray1D<int> pPatchIx;
			pPatchIx.AttachToData(&(m_vecRecvBuffer[0]));

			if (pPatchIx[0] > m_vecGridPatchByteLoc.GetRows()) {
				_EXCEPTION2("PatchIndex (%i) out of range [0,%i)",
					pPatchIx[0], m_vecGridPatchByteLoc.GetRows());
			}
/*
			if (iDataType == 0) {
				printf("Recv Geo %i %i\n",
					pPatchIx[0],
					m_vecGridPatchByteSize[pPatchIx[0]][iDataType]);
			} else {
				printf("Recv AcS %i %i\n",
					pPatchIx[0],
					m_vecGridPatchByteSize[pPatchIx[0]][iDataType]);
			}
*/
			m_ofsActiveOutput.seekp(
				posRefFile + m_vecGridPatchByteLoc[pPatchIx[0]][iDataType]);
			m_ofsActiveOutput.write(
				(const char *)(&(m_vecRecvBuffer[0])),
				m_vecGridPatchByteSize[pPatchIx[0]][iDataType]);
		}
	}

	// Barrier
	MPI_Barrier(MPI_COMM_WORLD);
#else
	_EXCEPTIONT("Not implemented without USE_MPI");
#endif
}

///////////////////////////////////////////////////////////////////////////////

Time OutputManagerComposite::Input(
	const std::string & strFileName
) {

	// Set the flag indicating that output came from a restart file
	m_fFromRestartFile = true;

	// Determine processor rank
	int nRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

	// The active model
	const Model & model = m_grid.GetModel();
/*
	// Open NetCDF file
	pNcFile = new NcFile(strFileName.c_str(), NcFile::ReadOnly);
	if (pNcFile == NULL) {
		_EXCEPTIONT("Error opening NetCDF file");
	}
	if (! pNcFile->is_valid()) {
		_EXCEPTIONT("Error opening NetCDF file");
	}

	// Get the time
	NcAtt * attCurrentTime = pNcFile->get_att("current_time");
	if (attCurrentTime == NULL) {
		_EXCEPTIONT("Attribute \'current_time\' not found in restart file");
	}
	std::string strCurrentTime = attCurrentTime->as_string(0);

	Time timeCurrent(strCurrentTime);

	// Equation set
	const EquationSet & eqn = m_grid.GetModel().GetEquationSet();

	for (int n = 0; n < m_grid.GetActivePatchCount(); n++) {
		GridPatch * pPatch = m_grid.GetActivePatch(n);

		int nCumulative2DNodeIx =
			m_grid.GetCumulativePatch2DNodeIndex(pPatch->GetPatchIndex());

		int nPatchNodeCount = pPatch->GetTotalNodeCount2D();

		// Input topography here
		DataArray2D<double> & dataTopography =
			pPatch->GetTopography();

		NcVar * varTopography = pNcFile->get_var("ZS");
		if (varTopography == NULL) {
			_EXCEPTIONT("Cannot find variable \'ZS\' in file");
		}
		varTopography->set_cur(nCumulative2DNodeIx);

		varTopography->get(dataTopography[0], nPatchNodeCount);

		// Input Rayleigh strength here
		DataArray3D<double> & dataRayleighStrengthNode =
			pPatch->GetRayleighStrength(DataLocation_Node);
		DataArray3D<double> & dataRayleighStrengthREdge =
			pPatch->GetRayleighStrength(DataLocation_REdge);

		NcVar * varRayleighNode = pNcFile->get_var("Rayleigh_Node");
		NcVar * varRayleighREdge = pNcFile->get_var("Rayleigh_REdge");

		varRayleighNode->set_cur(
			nCumulative2DNodeIx * m_grid.GetRElements());
		varRayleighREdge->set_cur(
			nCumulative2DNodeIx * (m_grid.GetRElements()+1));

		varRayleighNode->get(
			dataRayleighStrengthNode[0][0],
			pPatch->GetTotalNodeCount(DataLocation_Node));

		varRayleighREdge->get(
			dataRayleighStrengthREdge[0][0],
			pPatch->GetTotalNodeCount(DataLocation_REdge));

		// Input state
		DataArray4D<double> & dataStateNode =
			pPatch->GetDataState(0, DataLocation_Node);
		DataArray4D<double> & dataStateREdge =
			pPatch->GetDataState(0, DataLocation_REdge);

		DataArray4D<double> & dataRefStateNode =
			pPatch->GetReferenceState(DataLocation_Node);
		DataArray4D<double> & dataRefStateREdge =
			pPatch->GetReferenceState(DataLocation_REdge);

		for (int c = 0; c < eqn.GetComponents(); c++) {

			DataLocation loc = m_grid.GetVarLocation(c);

			// Number of radial elements for this component
			int nRadialDegreesOfFreedom;
			if (loc == DataLocation_Node) {
				nRadialDegreesOfFreedom = m_grid.GetRElements();
			} else if (loc == DataLocation_REdge) {
				nRadialDegreesOfFreedom = m_grid.GetRElements()+1;
			} else {
				_EXCEPTION();
			}

			// Number of nodes on this patch for this component
			int nComponentSize = pPatch->GetTotalNodeCount(loc);

			// Load state variable
			std::string strComponentName = eqn.GetComponentShortName(c);

			NcVar * var = pNcFile->get_var(strComponentName.c_str());
			if (var == NULL) {
				_EXCEPTION1("Cannot find variable \'%s\' in file",
					strComponentName.c_str());
			}
			var->set_cur(nCumulative2DNodeIx * nRadialDegreesOfFreedom);

			if (loc == DataLocation_Node) {
				var->get(dataStateNode[c][0][0], nComponentSize);
			} else if (loc == DataLocation_REdge) {
				var->get(dataStateREdge[c][0][0], nComponentSize);
			} else {
				_EXCEPTION();
			}

			// Load reference state
			if (m_grid.HasReferenceState()) {

				// Load reference state on nodes
				std::string strRefComponentNameNode =
					strComponentName + "_RefNode";

				NcVar * varRefNode =
					pNcFile->get_var(strRefComponentNameNode.c_str());
				if (varRefNode == NULL) {
					_EXCEPTION1("Cannot find variable \'%s\' in file",
						strRefComponentNameNode.c_str());
				}
				varRefNode->set_cur(
					nCumulative2DNodeIx * m_grid.GetRElements());

				varRefNode->get(
					dataRefStateNode[c][0][0],
					pPatch->GetTotalNodeCount(DataLocation_Node));

				// Load reference state on radial edges
				std::string strRefComponentNameREdge =
					strComponentName + "_RefREdge";

				NcVar * varRefREdge =
					pNcFile->get_var(strRefComponentNameREdge.c_str());
				if (varRefREdge == NULL) {
					_EXCEPTION1("Cannot find variable \'%s\' in file",
						strRefComponentNameREdge.c_str());
				}
				varRefREdge->set_cur(
					nCumulative2DNodeIx * (m_grid.GetRElements()+1));

				varRefREdge->get(
					dataRefStateREdge[c][0][0],
					pPatch->GetTotalNodeCount(DataLocation_REdge));
			}
		}

		// Input tracers

	}

	// Close the file
	if (pNcFile != NULL) {
		delete pNcFile;
	}

	return timeCurrent;
*/
	_EXCEPTION();
}

///////////////////////////////////////////////////////////////////////////////

