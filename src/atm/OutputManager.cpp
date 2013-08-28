///////////////////////////////////////////////////////////////////////////////
///
///	\file    OutputManager.cpp
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

#include "OutputManager.h"

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

OutputManager::OutputManager(
	Grid & grid,
	double dOutputDeltaT,
	std::string strOutputDir,
	std::string strOutputFormat
) :
	m_grid(grid),
	m_fFromRecoveryFile(false),
	m_nOutputFileIx(-1),
	m_dLastOutputTime(0.0),
	m_dNextOutputTime(0.0),
	m_dOutputDeltaT(dOutputDeltaT),
	m_strOutputDir(strOutputDir),
	m_strOutputFormat(strOutputFormat)
{
	// Create the output directory
	int nRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

	if (nRank == 0) {
		mkdir(m_strOutputDir.c_str(), 0777);
	}
}

///////////////////////////////////////////////////////////////////////////////

bool OutputManager::IsOutputNeeded(
	double dTime
) {
	// Only output initial and final results if no output time is set
	if (m_dOutputDeltaT == 0.0) {
		return false;
	}

	// Check current time against next output time
	if (dTime < (1.0 - 1.0e-12) * m_dNextOutputTime) {
		return false;
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////

void OutputManager::ManageOutput(
	double dTime
) {
	// Notification
	Announce("Output (%i): %f", m_nOutputFileIx+1, dTime);

	// Output
	Output(dTime);

	// Update last output time
	m_dLastOutputTime = dTime;

	// Update next output time
	m_dNextOutputTime += m_dOutputDeltaT;
}

///////////////////////////////////////////////////////////////////////////////

void OutputManager::InitialOutput(
	double dTime
) {
	// Check if we were initialized from a recovery file
	if (m_fFromRecoveryFile) {
		Announce("Output (%i): %f (Initial; Output Suppressed)",
			m_nOutputFileIx+1, dTime);

		m_dNextOutputTime = dTime + m_dOutputDeltaT;

	// Output initial conditions
	} else {
		Announce("Output (%i): %f (Initial)", m_nOutputFileIx+1, dTime);

		Output(dTime);

		m_dNextOutputTime = m_dOutputDeltaT;
	}
}

///////////////////////////////////////////////////////////////////////////////

void OutputManager::FinalOutput(
	double dTime
) {
	// Verify that the time has changed sufficiently since last output
	if (fabs(dTime - m_dLastOutputTime) < dTime * DBL_EPSILON) {
		return;
	}

	// Notification
	Announce("Output (%i): %f (Final)", m_nOutputFileIx+1, dTime);

	// Output a new file if necessary
	Output(dTime);
}

///////////////////////////////////////////////////////////////////////////////

void OutputManager::Output(
	double dTime
) {
	// Increment the current file index
	m_nOutputFileIx++;
}

///////////////////////////////////////////////////////////////////////////////

