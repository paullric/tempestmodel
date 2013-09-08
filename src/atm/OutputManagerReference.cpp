///////////////////////////////////////////////////////////////////////////////
///
///	\file    OutputManagerReference.cpp
///	\author  Paul Ullrich
///	\version May 7, 2013
///
///	<remarks>
///		Copyright 2000-2010 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks

#include "OutputManagerReference.h"

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

OutputManagerReference::OutputManagerReference(
	Grid & grid,
	double dOutputDeltaT,
	std::string strOutputDir,
	std::string strOutputFormat,
	int nXReference,
	int nYReference,
	bool fXNodal,
	bool fYNodal,
	bool fOutputVorticity,
	bool fOutputDivergence
) :
	OutputManager(
		grid,
		dOutputDeltaT,
		strOutputDir,
		strOutputFormat),
	m_pActiveNcOutput(NULL),
	m_fOutputVorticity(fOutputVorticity),
	m_fOutputDivergence(fOutputDivergence)
{
	// Create the coordinate arrays
	m_dXCoord.Initialize(nXReference);
	m_dYCoord.Initialize(nYReference);

	double dDeltaX = 360.0 / static_cast<double>(nXReference);

	// Initialize the coordinate arrays
	if (fXNodal) {
		for (int i = 0; i < nXReference; i++) {
			m_dXCoord[i] =
				dDeltaX * static_cast<double>(i) - 180.0;
		}
	} else {
		for (int i = 0; i < nXReference; i++) {
			m_dXCoord[i] =
				dDeltaX * (static_cast<double>(i) + 0.5) - 180.0;
		}
	}

	if (fYNodal) {
		double dDeltaY = 180.0 / static_cast<double>(nYReference-1);
		for (int i = 0; i < nYReference; i++) {
			m_dYCoord[i] =
				dDeltaY * static_cast<double>(i) - 90.0;
		}

	} else {
		double dDeltaY = 180.0 / static_cast<double>(nYReference);
		for (int i = 0; i < nYReference; i++) {
			m_dYCoord[i] =
				dDeltaY * (static_cast<double>(i) + 0.5) - 90.0;
		}
	}

	// Construct array of reference coordinates
	DataVector<double> dXReference;
	dXReference.Initialize(nXReference * nYReference);

	DataVector<double> dYReference;
	dYReference.Initialize(nXReference * nYReference);

	int ix = 0;
	for (int j = 0; j < nYReference; j++) {
	for (int i = 0; i < nXReference; i++) {
		dXReference[ix] = m_dXCoord[i] * M_PI / 180.0;
		dYReference[ix] = m_dYCoord[j] * M_PI / 180.0;
		ix++;
	}
	}

	grid.ConvertReferenceToABP(
		dXReference,
		dYReference,
		m_dAlpha,
		m_dBeta,
		m_iPanel);

	// Allocate data arrays
	m_dataRefState.Initialize(
		grid.GetModel().GetEquationSet().GetComponents(),
		grid.GetRElements(),
		nXReference * nYReference);

	m_dataState.Initialize(
		grid.GetModel().GetEquationSet().GetComponents(),
		grid.GetRElements(),
		nXReference * nYReference);

	if (grid.GetModel().GetEquationSet().GetTracers() != 0) {
		m_dataTracers.Initialize(
			grid.GetModel().GetEquationSet().GetTracers(),
			grid.GetRElements(),
			nXReference * nYReference);
	}
}

///////////////////////////////////////////////////////////////////////////////

OutputManagerReference::~OutputManagerReference() {
	DeinitializeNcOutput();
}

///////////////////////////////////////////////////////////////////////////////

void OutputManagerReference::InitializeNcOutput(
	const std::string & strFileName
) {
	// Determine processor rank; only proceed if root node
	int nRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

	// The active model
	const Model & model = m_grid.GetModel();

	// Open NetCDF file on root process
	if (nRank == 0) {

		// Check for existing NetCDF file
		if (m_pActiveNcOutput != NULL) {
			_EXCEPTIONT("NetCDF file already open");
		}

		// Open new NetCDF file
		m_pActiveNcOutput = new NcFile(strFileName.c_str(), NcFile::Replace);
		if (m_pActiveNcOutput == NULL) {
			_EXCEPTIONT("Error opening NetCDF file");
		}

		// Create nodal time dimension
		NcDim * dimTime =
			m_pActiveNcOutput->add_dim("time");

		m_varTime = m_pActiveNcOutput->add_var("time", ncDouble, dimTime);

		// Create levels dimension
		NcDim * dimLev =
			m_pActiveNcOutput->add_dim("lev", m_grid.GetRElements());

		// Create interfaces dimension
		NcDim * dimILev =
			m_pActiveNcOutput->add_dim("ilev", m_grid.GetRElements()+1);

		// Create latitude dimension
		NcDim * dimLat =
			m_pActiveNcOutput->add_dim("lat", m_dYCoord.GetRows());

		// Create longitude dimension
		NcDim * dimLon =
			m_pActiveNcOutput->add_dim("lon", m_dXCoord.GetRows());

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
					ncDouble, dimTime, dimLev, dimLat, dimLon));
		}

		for (int c = 0; c < eqn.GetTracers(); c++) {
			m_vecTracersVar.push_back(
				m_pActiveNcOutput->add_var(
					eqn.GetTracerShortName(c).c_str(),
					ncDouble, dimTime, dimLev, dimLat, dimLon));
		}

		// Divergence variable
		if (m_fOutputVorticity) {
			m_varVorticity =
				m_pActiveNcOutput->add_var(
					"ZETA", ncDouble, dimTime, dimLev, dimLat, dimLon);
		}

		// Divergence variable
		if (m_fOutputDivergence) {
			m_varDivergence =
				m_pActiveNcOutput->add_var(
					"DELTA", ncDouble, dimTime, dimLev, dimLat, dimLon);
		}

		// Output longitudes and latitudes
		NcVar * varLon = m_pActiveNcOutput->add_var("lon", ncDouble, dimLon);
		NcVar * varLat = m_pActiveNcOutput->add_var("lat", ncDouble, dimLat);

		varLon->put(m_dXCoord, m_dXCoord.GetRows());
		varLat->put(m_dYCoord, m_dYCoord.GetRows());

		// Output levels
		NcVar * varLev =
			m_pActiveNcOutput->add_var("lev", ncDouble, dimLev);

		varLev->put(
			m_grid.GetREtaLevels(),
			m_grid.GetREtaLevels().GetRows());

		// Output levels
		NcVar * varILev =
			m_pActiveNcOutput->add_var("ilev", ncDouble, dimILev);

		varILev->put(
			m_grid.GetREtaInterfaces(),
			m_grid.GetREtaInterfaces().GetRows());
	}

	// Wait for all processes to complete
	MPI_Barrier(MPI_COMM_WORLD);
}	

///////////////////////////////////////////////////////////////////////////////

void OutputManagerReference::DeinitializeNcOutput() {
	if (m_pActiveNcOutput != NULL) {
		delete(m_pActiveNcOutput);
		m_pActiveNcOutput = NULL;

		m_vecComponentVar.clear();
		m_vecTracersVar.clear();
	}
}

///////////////////////////////////////////////////////////////////////////////

void OutputManagerReference::Output(
	double dTime
) {
	int nRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

	// Call up the stack
	OutputManager::Output(dTime);

	// Equation set
	const EquationSet & eqn = m_grid.GetModel().GetEquationSet();

	// Add new time
	if (nRank == 0) {
		m_varTime->set_cur(m_nOutputFileIx);
		m_varTime->put(&dTime, 1);
	}

	// Vertically interpolate data to model levels
	for (int c = 0; c < eqn.GetComponents(); c++) {
		if (m_grid.GetVarLocation(c) == DataLocation_REdge) {
			m_grid.InterpolateREdgeToNode(c, 0);
		}
	}

	// Perform Interpolate / Reduction on state data
	m_dataState.Zero();

	m_grid.ReduceInterpolate(
		m_dAlpha, m_dBeta, m_iPanel,
		DataType_State, m_dataState, true);

	// Perform Interpolate / Reduction on tracers data
	if (m_grid.GetModel().GetEquationSet().GetTracers() != 0) {
		m_dataTracers.Zero();

		m_grid.ReduceInterpolate(
			m_dAlpha, m_dBeta, m_iPanel,
			DataType_Tracers, m_dataTracers, true);
	}

	// Perform Interpolate / Reduction on computed vorticity
	if (m_fOutputVorticity || m_fOutputDivergence) {
		if (m_dataState.GetRows() < 2) {
			_EXCEPTIONT("Insufficient components");
		}

		m_grid.ComputeVorticityDivergence(0);

		if (m_fOutputVorticity) {
			m_grid.ReduceInterpolate(
				m_dAlpha, m_dBeta, m_iPanel,
				DataType_Vorticity, m_dataVorticity);
		}
		if (m_fOutputDivergence) {
			m_grid.ReduceInterpolate(
				m_dAlpha, m_dBeta, m_iPanel,
				DataType_Divergence, m_dataDivergence);
		}
	}

	// Store state variable data
	if (nRank == 0) {
		for (int c = 0; c < eqn.GetComponents(); c++) {
			m_vecComponentVar[c]->set_cur(m_nOutputFileIx, 0, 0, 0);
			m_vecComponentVar[c]->put(
				&(m_dataState[c][0][0]),
				1,
				m_dataState.GetColumns(),
				m_dYCoord.GetRows(),
				m_dXCoord.GetRows());
		}

		// Store tracer variable data
		if (m_grid.GetModel().GetEquationSet().GetTracers() != 0) {
			for (int c = 0; c < eqn.GetTracers(); c++) {
				m_vecTracersVar[c]->set_cur(m_nOutputFileIx, 0, 0, 0);
				m_vecTracersVar[c]->put(
					&(m_dataTracers[c][0][0]),
					1,
					m_dataTracers.GetColumns(),
					m_dYCoord.GetRows(),
					m_dXCoord.GetRows());
			}
		}

		// Store vorticity data
		if (m_fOutputVorticity) {
			m_varVorticity->set_cur(m_nOutputFileIx, 0, 0, 0);
			m_varVorticity->put(
				&(m_dataVorticity[0][0][0]),
				1,
				m_dataVorticity.GetColumns(),
				m_dYCoord.GetRows(),
				m_dXCoord.GetRows());
		}

		// Store divergence data
		if (m_fOutputDivergence) {
			m_varDivergence->set_cur(m_nOutputFileIx, 0, 0, 0);
			m_varDivergence->put(
				&(m_dataDivergence[0][0][0]),
				1,
				m_dataDivergence.GetColumns(),
				m_dYCoord.GetRows(),
				m_dXCoord.GetRows());
		}
	}

	// Barrier
	MPI_Barrier(MPI_COMM_WORLD);
}

///////////////////////////////////////////////////////////////////////////////

