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

#ifdef TEMPEST_MPIOMP
#include <mpi.h>
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
		1)
{
	m_iCheck = 171456;
}

///////////////////////////////////////////////////////////////////////////////

OutputManagerComposite::~OutputManagerComposite() {
	CloseFile();
}

///////////////////////////////////////////////////////////////////////////////

bool OutputManagerComposite::OpenFile(
	const std::string & strFileName
) {
#ifdef TEMPEST_MPIOMP
	// Determine processor rank
	int nRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

	// Open file at root
	if (nRank == 0) {

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
	}

	// Wait for all processes to complete
	MPI_Barrier(MPI_COMM_WORLD);

	return true;
#else
	_EXCEPTIONT("Not implemented without TEMPEST_MPIOMP");
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
#ifdef TEMPEST_MPIOMP
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

	// Determine space allocation for each GridPatch
	m_vecGridPatchByteSize.Allocate(m_grid.GetPatchCount(), 2);

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

	// Reduce GridPatch byte size and write Grid information
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

		// Write check bits
		m_ofsActiveOutput.write((const char *)(&m_iCheck), sizeof(int));

		// Write current time
		const Time & timeCurrent = model.GetCurrentTime();
		m_ofsActiveOutput.write((const char *)(&(timeCurrent)), sizeof(Time));

		// Write Grid information to file
		const DataContainer & dcGridParameters =
			m_grid.GetDataContainerParameters();
		int nGridParametersByteSize =
			dcGridParameters.GetTotalByteSize();
		const char * pGridParameters =
			(const char *)(dcGridParameters.GetPointer());

		m_ofsActiveOutput.write(
			pGridParameters, nGridParametersByteSize);

		const DataContainer & dcGridPatchData =
			m_grid.GetDataContainerPatchData();
		int nGridPatchDataByteSize =
			dcGridPatchData.GetTotalByteSize();
		const char * pGridPatchData =
			(const char *)(dcGridPatchData.GetPointer());

		m_ofsActiveOutput.write(
			pGridPatchData, nGridPatchDataByteSize);
	}

	// Send data from GridPatches to root
	if (nRank != 0) {

		int nActivePatches = m_grid.GetActivePatchCount();
		if (nActivePatches == 0) {
			_EXCEPTIONT("No GridPatches on processor");
		}

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
				const_cast<unsigned char *>(pGeometricData),
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
				const_cast<unsigned char *>(pActiveStateData),
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
	} else {

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
			int iPatchIx = *((int*)(&(m_vecRecvBuffer[0])));

			if (iPatchIx > m_vecGridPatchByteLoc.GetRows()) {
				_EXCEPTION2("PatchIndex (%i) out of range [0,%i)",
					iPatchIx, m_vecGridPatchByteLoc.GetRows());
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
				posRefFile + m_vecGridPatchByteLoc[iPatchIx][iDataType]);
			m_ofsActiveOutput.write(
				(const char *)(&(m_vecRecvBuffer[0])),
				m_vecGridPatchByteSize[iPatchIx][iDataType]);
		}
	}

	// Barrier
	MPI_Barrier(MPI_COMM_WORLD);

#else
	_EXCEPTIONT("Not implemented without TEMPEST_MPIOMP");
#endif
}

///////////////////////////////////////////////////////////////////////////////

Time OutputManagerComposite::Input(
	const std::string & strFileName
) {
#ifdef TEMPEST_MPIOMP
	// Set the flag indicating that output came from a restart file
	m_fFromRestartFile = true;

	// Determine processor rank
	int nRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

	// The active model
	const Model & model = m_grid.GetModel();

	// Open binary input stream
	std::ifstream ifsActiveInput;
	ifsActiveInput.open(
		strFileName.c_str(), std::ios::binary | std::ios::in);

	if (!ifsActiveInput) {
		_EXCEPTION1("Unable to open input file \"%s\"",
			strFileName.c_str());
	}

	// Read check bits
	int iCheckInput;
	ifsActiveInput.read((char *)(&iCheckInput), sizeof(int));
	if (iCheckInput != m_iCheck) {
		_EXCEPTION1("Invalid or incompatible input file \"%s\"",
			strFileName.c_str());
	}

	// Read current time
	Time timeCurrent;
	ifsActiveInput.read((char *)(&(timeCurrent)), sizeof(Time));

	// Read Grid parameters from file
	DataContainer & dcGridParameters = m_grid.GetDataContainerParameters();
	int nGridParametersByteSize =
		dcGridParameters.GetTotalByteSize();
	char * pGridParameters =
		(char *)(dcGridParameters.GetPointer());

	ifsActiveInput.read(pGridParameters, nGridParametersByteSize);

	// Initialize the Grid from specified parameters
	m_grid.InitializeDataLocal();

	// Load Grid data from file
	DataContainer & dcGridPatchData = m_grid.GetDataContainerPatchData();
	int nGridPatchDataByteSize =
		dcGridPatchData.GetTotalByteSize();
	char * pGridPatchData =
		(char *)(dcGridPatchData.GetPointer());

	ifsActiveInput.read(pGridPatchData, nGridPatchDataByteSize);

	// Distribute GridPatches to processors
	m_grid.DistributePatches();

	// Determine space allocation for each GridPatch
	m_vecGridPatchByteSize.Allocate(m_grid.GetPatchCount(), 2);

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

	MPI_Allreduce(
		MPI_IN_PLACE,
		&(m_vecGridPatchByteSize[0][0]),
		m_vecGridPatchByteSize.GetTotalSize(),
		MPI_INT,
		MPI_MAX,
		MPI_COMM_WORLD);

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

	// Reference position
	std::streampos posRefFile = ifsActiveInput.tellg();
	if (posRefFile == (-1)) {
		_EXCEPTIONT("ActiveInput::tellg() fail");
	}

	// Load in GridPatch data from file
	for (int i = 0; i < m_grid.GetActivePatchCount(); i++) {
		GridPatch * pPatch = m_grid.GetActivePatch(i);

		int iPatchIx = pPatch->GetPatchIndex();
		if (iPatchIx > m_grid.GetPatchCount()) {
			_EXCEPTION2("PatchIndex (%i) out of range [0,%i)",
				iPatchIx, m_grid.GetPatchCount());
		}

		DataContainer & dcGeometric =
			pPatch->GetDataContainerGeometric();
		int nGeometricDataByteSize =
			dcGeometric.GetTotalByteSize();
		char * pGeometricData =
			(char *)(dcGeometric.GetPointer());

		ifsActiveInput.seekg(
			posRefFile + m_vecGridPatchByteLoc[iPatchIx][0]);
		ifsActiveInput.read(
			pGeometricData, m_vecGridPatchByteSize[iPatchIx][0]);

		DataContainer & dcActiveState =
			pPatch->GetDataContainerActiveState();
		int nActiveStateDataByteSize =
			dcActiveState.GetTotalByteSize();
		char * pActiveStateData =
			(char *)(dcActiveState.GetPointer());

		ifsActiveInput.seekg(
			posRefFile + m_vecGridPatchByteLoc[iPatchIx][1]);
		ifsActiveInput.read(
			pActiveStateData, m_vecGridPatchByteSize[iPatchIx][1]);
	}

	// Close the file
	ifsActiveInput.close();

	// Barrier
	MPI_Barrier(MPI_COMM_WORLD);

#else
	_EXCEPTIONT("Not implemented without TEMPEST_MPIOMP");
#endif

	return timeCurrent;
}

///////////////////////////////////////////////////////////////////////////////

