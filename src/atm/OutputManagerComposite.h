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
#include "DataArray1D.h"
#include "DataArray2D.h"
#include <fstream>

class Time;

///////////////////////////////////////////////////////////////////////////////

#ifndef TEMPEST_NETCDF
typedef int NcFile;
typedef int NcVar;
#endif

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		An OutputManager which implements direct dumping of data structures
///		to a NetCDF file for later recovery.
///	</summary>
class OutputManagerComposite :
	public OutputManager
{

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	OutputManagerComposite(
		Grid & grid,
		const Time & timeOutputFrequency,
		std::string strOutputDir,
		std::string strOutputPrefix,
		std::string strRestartFile = ""
	);

	///	<summary>
	///		Destructor.
	///	</summary>
	virtual ~OutputManagerComposite();

	///	<summary>
	///		Get the name of the OutputManager.
	///	</summary>
	virtual const char * GetName() const {
		return "Composite";
	}

public:
	///	<summary>
	///		Open a new NetCDF file.
	///	</summary>
	virtual bool OpenFile(
		const std::string & strFileName
	);

	///	<summary>
	///		Close an existing NetCDF file.
	///	</summary>
	virtual void CloseFile();

protected:
	///	<summary>
	///		Write output to a file.
	///	</summary>
	void Output(
		const Time & time
	);

protected:
	///	<summary>
	///		Returns true if this OutputManager supports the input operation.
	///	</summary>
	virtual bool SupportsInput() const {
		return true;
	}

public:
	///	<summary>
	///		Initialize the grid data from a file.
	///	</summary>
	virtual Time Input(
		const std::string & strFileName
	);

protected:
	///	<summary>
	///		Check bits.
	///	</summary>
	int m_iCheck;

protected:
	///	<summary>
	///		Active output file.
	///	</summary>
	std::ofstream m_ofsActiveOutput;

private:
	///	<summary>
	///		Byte size for each GridPatch.
	///	</summary>
	DataArray2D<int> m_vecGridPatchByteSize;

	///	<summary>
	///		Byte location for each GridPatch.
	///	</summary>
	DataArray2D<std::streamoff> m_vecGridPatchByteLoc;

	///	<summary>
	///		Vector used at root node for storage of received data.
	///	</summary>
	DataArray1D<double> m_vecRecvBuffer;
};

///////////////////////////////////////////////////////////////////////////////

#endif

