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

#include "DataVector.h"
#include "DataMatrix4D.h"

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
	///		Types of data that can be handled by ConsolidateData.
	///	</summary>
	enum OutputDataType {
		OutputDataType_Default = 0,
		OutputDataType_Data = OutputDataType_Default,
		OutputDataType_Tracers = 1,
		OutputDataType_CentroidZ = 2
	};

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	OutputManager(
		Grid & grid,
		double dOutputDeltaT,
		std::string strOutputDir,
		std::string strOutputFormat
	);

	///	<summary>
	///		Destructor.
	///	</summary>
	virtual ~OutputManager()
	{ }

public:
	///	<summary>
	///		Determine if an output is needed.
	///	</summary>
	bool IsOutputNeeded(double dTime);

	///	<summary>
	///		Perform an output.
	///	</summary>
	void ManageOutput(double dTime);

	///	<summary>
	///		Write the initial system state to a file.
	///	</summary>
	void InitialOutput(double dTime);

	///	<summary>
	///		Write the final system state to a file.
	///	</summary>
	void FinalOutput(double dTime);

protected:
	///	<summary>
	///		Handle manager-specific file output.
	///	</summary>
	virtual void Output(
		double dTime
	);

protected:
	///	<summary>
	///		Grid associated with this OutputManager.
	///	</summary>
	Grid & m_grid;

	///	<summary>
	///		Flag indicating if we were initialized from a recovery file.
	///	</summary>
	bool m_fFromRecoveryFile;

	///	<summary>
	///		Current output file index.
	///	</summary>
	int m_nOutputFileIx;

	///	<summary>
	///		Last output time.
	///	</summary>
	double m_dLastOutputTime;

	///	<summary>
	///		Next output time.
	///	</summary>
	double m_dNextOutputTime;

	///	<summary>
	///		Time between successive outputs.
	///	</summary>
	double m_dOutputDeltaT;

	///	<summary>
	///		Output directory.
	///	</summary>
	std::string m_strOutputDir;

	///	<summary>
	///		Format of the output filename.
	///	</summary>
	std::string m_strOutputFormat;
};

///////////////////////////////////////////////////////////////////////////////

#endif

