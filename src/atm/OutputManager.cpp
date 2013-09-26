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
	std::string strOutputPrefix,
	int nOutputsPerFile
) :
	m_grid(grid),
	m_fFromRecoveryFile(false),
	m_fIsFileOpen(false),
	m_ixOutputTime(0),
	m_ixOutputFile(0),
	m_dLastOutputTime(0.0),
	m_dNextOutputTime(0.0),
	m_dOutputDeltaT(dOutputDeltaT),
	m_strOutputDir(strOutputDir),
	m_strOutputPrefix(strOutputPrefix),
	m_nOutputsPerFile(nOutputsPerFile)
{
	// Create the output directory
	int nRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

	if ((nRank == 0) && (m_strOutputDir != "")) {
		mkdir(m_strOutputDir.c_str(), 0777);
	}
}

///////////////////////////////////////////////////////////////////////////////

void OutputManager::GetFileName(std::string & strFileName) const {
	char szIndex[128];
	sprintf(szIndex, "%06i", m_ixOutputFile);
	
	strFileName = m_strOutputDir;
	strFileName += "/";
	strFileName += m_strOutputPrefix;
	strFileName += szIndex;
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

void OutputManager::PerformOutput(
	double dTime
) {
	// Open the file
	if (!m_fIsFileOpen) {
		std::string strActiveFileName;
		GetFileName(strActiveFileName);
		m_fIsFileOpen = OpenFile(strActiveFileName);

		if (!m_fIsFileOpen) {
			_EXCEPTION1("Unable to open file \"%s\"",
				strActiveFileName.c_str());
		}
	}

	// Output
	Output(dTime);
	m_ixOutputTime ++;

	// Check if time limit is reached
	if (m_ixOutputTime == m_nOutputsPerFile) {
		CloseFile();
		m_fIsFileOpen = false;
		m_ixOutputFile ++;
		m_ixOutputTime = 0;
	}
}

///////////////////////////////////////////////////////////////////////////////

void OutputManager::ManageOutput(
	double dTime
) {
	// Notification
	Announce("Output (%i): %f", m_ixOutputFile+1, dTime);

	// Perform the output
	PerformOutput(dTime);

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
			m_ixOutputFile+1, dTime);

		m_dNextOutputTime = dTime + m_dOutputDeltaT;

	// Output initial conditions
	} else {
		Announce("Output (%i/%i): %f (Initial)",
			m_ixOutputFile+1, m_ixOutputTime+1, dTime);

		PerformOutput(dTime);

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
	Announce("Output (%i/%i): %f (Final)",
		m_ixOutputFile+1, m_ixOutputTime+1, dTime);

	// Output
	PerformOutput(dTime);
}

///////////////////////////////////////////////////////////////////////////////

