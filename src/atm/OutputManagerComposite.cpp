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

#include "Announce.h"

#include "mpi.h"

#include <iostream>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <sys/stat.h>

///////////////////////////////////////////////////////////////////////////////

OutputManagerComposite::OutputManagerComposite(
	Grid & grid,
	double dOutputDeltaT,
	std::string strOutputDir,
	std::string strOutputFormat,
	std::string strRestartFile
) :
	OutputManager(
		grid,
		dOutputDeltaT,
		strOutputDir,
		strOutputFormat,
		1),
	m_pActiveNcOutput(NULL),
	m_strRestartFile(strRestartFile)
{
}

///////////////////////////////////////////////////////////////////////////////

OutputManagerComposite::~OutputManagerComposite() {
	CloseFile();
}

///////////////////////////////////////////////////////////////////////////////

void OutputManagerComposite::ResizeDataBuffer() {
	m_vecLocalData.Initialize(
		m_grid.GetMaximumDegreesOfFreedom());
}

///////////////////////////////////////////////////////////////////////////////

bool OutputManagerComposite::OpenFile(
	const std::string & strFileName
) {
	// Determine processor rank; only proceed if root node
	int nRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

	// The active model
	const Model & model = m_grid.GetModel();

	// Open new NetCDF file
	NcVar * varZs;

	if (nRank == 0) {

		// Check for existing NetCDF file
		if (m_pActiveNcOutput != NULL) {
			_EXCEPTIONT("NetCDF file already open");
		}

		// Open new NetCDF file
		std::string strNcFileName = strFileName + ".restart.nc";
		m_pActiveNcOutput = new NcFile(strNcFileName.c_str(), NcFile::Replace);
		if (m_pActiveNcOutput == NULL) {
			_EXCEPTIONT("Error opening NetCDF file");
		}

		// Create nodal time dimension
		NcDim * dimTime =
			m_pActiveNcOutput->add_dim("time");

		m_varTime = m_pActiveNcOutput->add_var("time", ncDouble, dimTime);

		// Create patch index dimension
		NcDim * dimPatchIndex =
			m_pActiveNcOutput->add_dim(
				"patch_index", m_grid.GetPatchCount());

		// Create nodal index dimension
		NcDim * dimIndex =
			m_pActiveNcOutput->add_dim(
				"node_index", m_grid.GetTotalNodeCount());

		// Output physical constants
		const PhysicalConstants & phys = model.GetPhysicalConstants();

		m_pActiveNcOutput->add_att("earth_radius", phys.GetEarthRadius());
		m_pActiveNcOutput->add_att("g", phys.GetG());
		m_pActiveNcOutput->add_att("omega", phys.GetOmega());
		m_pActiveNcOutput->add_att("alpha", phys.GetAlpha());
		m_pActiveNcOutput->add_att("Rd", phys.GetR());
		m_pActiveNcOutput->add_att("Cp", phys.GetCp());
		m_pActiveNcOutput->add_att("T0", phys.GetT0());
		m_pActiveNcOutput->add_att("P0", phys.GetP0());
		m_pActiveNcOutput->add_att("rho_water", phys.GetRhoWater());
		m_pActiveNcOutput->add_att("Rvap", phys.GetRvap());
		m_pActiveNcOutput->add_att("Mvap", phys.GetMvap());
		m_pActiveNcOutput->add_att("Lvap", phys.GetLvap());

		// Output equation set
		const EquationSet & eqn = model.GetEquationSet();

		m_pActiveNcOutput->add_att("equation_set", eqn.GetName().c_str());

		// Create variables
		for (int c = 0; c < eqn.GetComponents(); c++) {
			m_vecComponentVar.push_back(
				m_pActiveNcOutput->add_var(
					eqn.GetComponentShortName(c).c_str(),
					ncDouble, dimTime, dimIndex));
		}

		for (int c = 0; c < eqn.GetTracers(); c++) {
			m_vecTracersVar.push_back(
				m_pActiveNcOutput->add_var(
					eqn.GetTracerShortName(c).c_str(),
					ncDouble, dimTime, dimIndex));
		}

		// Output PatchBox for each patch
		NcDim * dimPatchBox = m_pActiveNcOutput->add_dim("patchbox_info", 7);

		NcVar * pPanel =
			m_pActiveNcOutput->add_var(
				"PatchBox", ncInt, dimPatchIndex, dimPatchBox);

		for (int n = 0; n < m_grid.GetPatchCount(); n++) {
			int nBox[7];

			const PatchBox & box = m_grid.GetPatch(n)->GetPatchBox();

			nBox[0] = box.GetPanel();
			nBox[1] = box.GetRefinementLevel();
			nBox[2] = box.GetHaloElements();
			nBox[3] = box.GetAGlobalInteriorBegin();
			nBox[4] = box.GetAGlobalInteriorEnd();
			nBox[5] = box.GetBGlobalInteriorBegin();
			nBox[6] = box.GetBGlobalInteriorEnd();

			pPanel->set_cur(n, 0);
			pPanel->put(nBox, 1, 7);
		}

		// Output topography for each patch
		varZs = m_pActiveNcOutput->add_var("ZS", ncDouble, dimIndex);
	}

    // Begin data consolidation
    std::vector<DataType> vecDataTypes;
    vecDataTypes.push_back(DataType_Topography);

    ConsolidationStatus status(m_grid, vecDataTypes);

    m_grid.ConsolidateDataToRoot(status);

    DataVector<double> dataRecvBuffer;
    dataRecvBuffer.Initialize(
        m_grid.GetLargestGridPatchNodes()
        * m_grid.GetRElements());

    while ((nRank == 0) && (!status.Done())) {

        int nRecvCount;
        int ixRecvPatch;
        DataType eRecvDataType;

        m_grid.ConsolidateDataAtRoot(
            status,
            dataRecvBuffer,
            nRecvCount,
            ixRecvPatch,
            eRecvDataType);

        // Store topography data
        if (eRecvDataType == DataType_Topography) {
            int nPatchIx = m_grid.GetCumulativePatch2DNodeIndex(ixRecvPatch);
            varZs->set_cur(nPatchIx);
            varZs->put(&(dataRecvBuffer[0]), nRecvCount);

		} else {
			_EXCEPTIONT("Invalid DataType");
		}
	}

	// Wait for all processes to complete
	MPI_Barrier(MPI_COMM_WORLD);

	return true;
}	

///////////////////////////////////////////////////////////////////////////////

void OutputManagerComposite::CloseFile() {
	if (m_pActiveNcOutput != NULL) {
		delete(m_pActiveNcOutput);
		m_pActiveNcOutput = NULL;

		m_vecComponentVar.clear();
		m_vecTracersVar.clear();
	}
}

///////////////////////////////////////////////////////////////////////////////

void OutputManagerComposite::Output(
	double dTime
) {
	// Check for open file
	if (!IsFileOpen()) {
		_EXCEPTIONT("No file available for output");
	}

	// Determine processor rank
	int nRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

	// Equation set
	const EquationSet & eqn = m_grid.GetModel().GetEquationSet();

	// Add new time
	if (nRank == 0) {
		m_varTime->set_cur(m_ixOutputTime);
		m_varTime->put(&dTime, 1);
	}

	// Begin data consolidation
	std::vector<DataType> vecDataTypes;
	vecDataTypes.push_back(DataType_State);
	if (eqn.GetTracers() != 0) {
		vecDataTypes.push_back(DataType_Tracers);
	}

	ConsolidationStatus status(m_grid, vecDataTypes);

	m_grid.ConsolidateDataToRoot(status);

	// Data buffer
	DataVector<double> dataRecvBuffer;
	if (nRank == 0) {
		dataRecvBuffer.Initialize(m_grid.GetMaximumDegreesOfFreedom());
	}

	// Receive all data objects from neighbors
	while ((nRank == 0) && (!status.Done())) {
		int nRecvCount;
		int ixRecvPatch;
		DataType eRecvDataType;

		m_grid.ConsolidateDataAtRoot(
			status,
			dataRecvBuffer,
			nRecvCount,
			ixRecvPatch,
			eRecvDataType);

		// Store state variable data
		if (eRecvDataType == DataType_State) {
			int nComponentSize =
				nRecvCount / eqn.GetComponents();

			if (nRecvCount % eqn.GetComponents() != 0) {
				_EXCEPTIONT("Invalid message length");
			}

			int nPatchIx = m_grid.GetCumulativePatch3DNodeIndex(ixRecvPatch);
			for (int c = 0; c < eqn.GetComponents(); c++) {
				m_vecComponentVar[c]->set_cur(m_ixOutputTime, nPatchIx);
				m_vecComponentVar[c]->put(
					&(dataRecvBuffer[nComponentSize * c]), 1, nComponentSize);
			}

		// Store tracer variable data
		} else if (eRecvDataType == DataType_Tracers) {
			int nComponentSize =
				nRecvCount / eqn.GetTracers();

			if (nRecvCount % eqn.GetTracers() != 0) {
				_EXCEPTIONT("Invalid message length");
			}

			int nPatchIx = m_grid.GetCumulativePatch3DNodeIndex(ixRecvPatch);
			for (int c = 0; c < eqn.GetTracers(); c++) {
				m_vecTracersVar[c]->set_cur(m_ixOutputTime, nPatchIx);
				m_vecTracersVar[c]->put(
					&(dataRecvBuffer[nComponentSize * c]), 1, nComponentSize);
			}
		}
	}

	// Barrier
	MPI_Barrier(MPI_COMM_WORLD);
}

///////////////////////////////////////////////////////////////////////////////

void OutputManagerComposite::Input(
	Grid & grid
) const {
	// Determine processor rank; only proceed if root node
	int nRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

	// The active model
	const Model & model = m_grid.GetModel();

	// Open new NetCDF file
	NcFile * pNcFile = NULL;
	NcVar * varZs;

	if (nRank == 0) {

		// Open NetCDF file
		pNcFile = new NcFile(m_strRestartFile.c_str(), NcFile::ReadOnly);
		if (m_pActiveNcOutput == NULL) {
			_EXCEPTIONT("Error opening NetCDF file");
		}

	}

	// Close the file
	if (pNcFile != NULL) {
		delete pNcFile;
	}
}

///////////////////////////////////////////////////////////////////////////////

