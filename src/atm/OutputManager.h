///////////////////////////////////////////////////////////////////////////////
///
///	\file    OutputManager.h
///	\author  Paul Ullrich
///	\version March 8, 2013
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

#ifndef _OUTPUTMANAGER_H_
#define _OUTPUTMANAGER_H_

#include "TimeObj.h"
#include "DataArray1D.h"
#include "DataArray4D.h"

#include "netcdfcpp.h"

#include <vector>

class Grid;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Interface class for handling output during the model simulation to
///		the user, either via standard output or to disk.
///	</summary>
class OutputManager {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	OutputManager(
		Grid & grid,
		const Time & timeOutputFrequency,
		std::string strOutputDir,
		std::string strOutputPrefix,
		int nOutputsPerFile
	);

	///	<summary>
	///		Destructor.
	///	</summary>
	virtual ~OutputManager()
	{ }

public:
	///	<summary>
	///		Get the name of this OutputManager.
	///	</summary>
	virtual const char * GetName() const {
		return "Output";
	}

protected:
	///	<summary>
	///		Get the active file name.
	///	</summary>
	void GetFileName(
		const Time & time,
		std::string & strFileName
	) const;

	///	<summary>
	///		Perform the output.
	///	</summary>
	void PerformOutput(const Time & time);

public:
	///	<summary>
	///		Determine if an output is needed.
	///	</summary>
	bool IsOutputNeeded(const Time & time);

	///	<summary>
	///		Perform an output.
	///	</summary>
	void ManageOutput(const Time & time);

	///	<summary>
	///		Write the initial system state to a file.
	///	</summary>
	void InitialOutput(const Time & time);

	///	<summary>
	///		Write the final system state to a file.
	///	</summary>
	void FinalOutput(const Time & time);

public:
	///	<summary>
	///		Returns true if this OutputManager supports the input operation.
	///	</summary>
	virtual bool SupportsInput() const {
		return false;
	}

	///	<summary>
	///		Load in the state from a file.
	///	</summary>
	virtual Time Input(
		const std::string & strFileName
	) {
		_EXCEPTIONT("OutputManager does not support Input operation");
	}

protected:
	///	<summary>
	///		Check if a file is currently open.
	///	</summary>
	inline bool IsFileOpen() const {
		return m_fIsFileOpen;
	}

	///	<summary>
	///		Open a new file.
	///	</summary>
	virtual bool OpenFile(
		const std::string & strFileName
	) {
		return true;
	}

	///	<summary>
	///		Close an existing file.
	///	</summary>
	virtual void CloseFile()
	{ }

	///	<summary>
	///		Handle manager-specific file output.
	///	</summary>
	virtual void Output(
		const Time & time
	) = 0;

protected:
	///	<summary>
	///		Grid associated with this OutputManager.
	///	</summary>
	Grid & m_grid;

	///	<summary>
	///		Flag indicating that the initial conditions came from a
	///		recovery file and that output should be supressed.
	///	</summary>
	bool m_fFromRestartFile;

	///	<summary>
	///		Flag indicating that a file is currently open.
	///	</summary>
	bool m_fIsFileOpen;

	///	<summary>
	///		Current output time index.
	///	</summary>
	int m_ixOutputTime;

	///	<summary>
	///		Current output file index.
	///	</summary>
	int m_ixOutputFile;

	///	<summary>
	///		Last output time.
	///	</summary>
	Time m_timeLastOutput;

	///	<summary>
	///		Next output time.
	///	</summary>
	Time m_timeNextOutput;

	///	<summary>
	///		Time between successive outputs (in seconds).
	///	</summary>
	Time m_timeOutputFrequency;

	///	<summary>
	///		Output directory.
	///	</summary>
	std::string m_strOutputDir;

	///	<summary>
	///		Format of the output filename.
	///	</summary>
	std::string m_strOutputPrefix;

	///	<summary>
	///		Number of time steps output per file.
	///	</summary>
	int m_nOutputsPerFile;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Vector of pointers to OutputManagers.
///	</summary>
typedef std::vector<OutputManager *> OutputManagerVector;

///////////////////////////////////////////////////////////////////////////////

#endif

