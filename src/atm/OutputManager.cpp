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

#include "TimeObj.h"
#include "Model.h"
#include "Grid.h"
#include "ConsolidationStatus.h"

#include "Announce.h"

#include <mpi.h>

#include <iostream>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <sys/stat.h>

///////////////////////////////////////////////////////////////////////////////

OutputManager::OutputManager(
	Grid & grid,
	const Time & timeOutputFrequency,
	std::string strOutputDir,
	std::string strOutputPrefix,
	int nOutputsPerFile
) :
	m_grid(grid),
	m_fFromRestartFile(false),
	m_fIsFileOpen(false),
	m_ixOutputTime(0),
	m_ixOutputFile(0),
	m_timeLastOutput(),
	m_timeNextOutput(),
	m_timeOutputFrequency(timeOutputFrequency),
	m_strOutputDir(strOutputDir),
	m_strOutputPrefix(strOutputPrefix),
	m_nOutputsPerFile(nOutputsPerFile)
{
#ifdef TEMPEST_MPIOMP
	// Create the output directory
	int nRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

	if ((nRank == 0) && (m_strOutputDir != "")) {
		mkdir(m_strOutputDir.c_str(), 0777);
	}
#endif
}

///////////////////////////////////////////////////////////////////////////////

void OutputManager::GetFileName(
	const Time & time,
	std::string & strFileName
) const {
	strFileName = m_strOutputDir;
	strFileName += "/";
	strFileName += m_strOutputPrefix;
	strFileName += ".";
	strFileName += time.ToShortString();
}

///////////////////////////////////////////////////////////////////////////////

bool OutputManager::IsOutputNeeded(
	const Time & time
) {
	// Only output initial and final results if no output time is set
	if (m_timeOutputFrequency.IsZero()) {
		return false;
	}

	// Check current time against next output time
	if (time < m_timeNextOutput) {
		return false;
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////

void OutputManager::PerformOutput(
	const Time & time
) {
	// Open the file
	if (!m_fIsFileOpen) {
		std::string strActiveFileName;
		GetFileName(time, strActiveFileName);
		m_fIsFileOpen = OpenFile(strActiveFileName);

		if (!m_fIsFileOpen) {
			_EXCEPTION1("Unable to open file \"%s\"",
				strActiveFileName.c_str());
		}
	}

	// Output
	Output(time);
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
	const Time & time
) {
	// Notification
	Announce("%s (%i/%i): %s",
		GetName(),
		m_ixOutputFile+1,
		m_ixOutputTime+1,
		time.ToString().c_str());

	// Perform the output
	PerformOutput(time);

	// Update last output time
	m_timeLastOutput = time;

	// Update next output time
	m_timeNextOutput += m_timeOutputFrequency;
}

///////////////////////////////////////////////////////////////////////////////

void OutputManager::InitialOutput(
	const Time & time
) {
	// Check if we were initialized from a recovery file
	if (m_fFromRestartFile) {
		Announce("%s (%i): %s (Initial; Output Suppressed)",
			GetName(),
			m_ixOutputFile+1,
			time.ToString().c_str());

	// Output initial conditions
	} else {
		Announce("%s (%i/%i): %s (Initial)",
			GetName(),
			m_ixOutputFile+1,
			m_ixOutputTime+1,
			time.ToString().c_str());

		PerformOutput(time);
	}

	m_timeNextOutput = time;
	m_timeNextOutput += m_timeOutputFrequency;
}

///////////////////////////////////////////////////////////////////////////////

void OutputManager::FinalOutput(
	const Time & time
) {
	// Verify that the time has changed sufficiently since last output
	if (time == m_timeLastOutput) {
		return;
	}

	// Notification
	Announce("%s (%i/%i): %s (Final)",
		GetName(),
		m_ixOutputFile+1,
		m_ixOutputTime+1,
		time.ToString().c_str());

	// Output
	PerformOutput(time);
}

///////////////////////////////////////////////////////////////////////////////

