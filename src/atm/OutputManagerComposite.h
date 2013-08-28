///////////////////////////////////////////////////////////////////////////////
///
///	\file    OutputManagerComposite.h
///	\author  Paul Ullrich
///	\version May 6, 2013
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

#ifndef _OUTPUTMANAGERCOMPOSITE_H_
#define _OUTPUTMANAGERCOMPOSITE_H_

#include "OutputManager.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		An OutputManager which implements direct dumping of data structures
///		to a NetCDF file for later recovery.
///	</summary>
class OutputManagerComposite : public OutputManager {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	OutputManagerComposite(
		Grid & grid,
		double dOutputDeltaT,
		std::string strOutputDir,
		std::string strOutputFormat
	);

	///	<summary>
	///		Destructor.
	///	</summary>
	virtual ~OutputManagerComposite();

public:
	///	<summary>
	///		Resize the data buffer to store maximum patch data.
	///	</summary>
	void ResizeDataBuffer();

public:
	///	<summary>
	///		Load in data from a recovery file.
	///	</summary>
	bool InitialInput(
		const std::string & strFilename,
		double & dTime
	);

public:
	///	<summary>
	///		Open a new NetCDF file.
	///	</summary>
	void InitializeNcOutput(
		const std::string & strFileName
	);

	///	<summary>
	///		Close an existing NetCDF file.
	///	</summary>
	void DeinitializeNcOutput();

protected:
	///	<summary>
	///		Write output to a file.
	///	</summary>
	void Output(
		double dTime
	);

protected:
	///	<summary>
	///		Active output file.
	///	</summary>
	NcFile * m_pActiveNcOutput;

	///	<summary>
	///		Time variable.
	///	</summary>
	NcVar * m_varTime;

	///	<summary>
	///		Vector of component variables.
	///	</summary>
	std::vector<NcVar *> m_vecComponentVar;

	///	<summary>
	///		Vector of tracer variables.
	///	</summary>
	std::vector<NcVar *> m_vecTracersVar;

private:
	///	<summary>
	///		Vectors used for temporary storage of data.
	///	</summary>
	DataVector<double> m_vecLocalData;
};

///////////////////////////////////////////////////////////////////////////////

#endif

