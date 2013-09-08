///////////////////////////////////////////////////////////////////////////////
///
///	\file    OutputManagerReference.h
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

#ifndef _OUTPUTMANAGERREFERENCE_H_
#define _OUTPUTMANAGERREFERENCE_H_

#include "OutputManager.h"

#include "DataMatrix3D.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		OutputManager that handles interpolated output to the reference grid
///		(RLL grid on the sphere, Cartesian grid on the plane).
///	</summary>
class OutputManagerReference : public OutputManager {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	OutputManagerReference(
		Grid & grid,
		double dOutputDeltaT,
		std::string strOutputDir,
		std::string strOutputFormat,
		int nXReference,
		int nYReference,
		bool fXNodal,
		bool fYNodal,
		bool fOutputVorticity = false,
		bool fOutputDivergence = false
	);

	///	<summary>
	///		Modify the flag which indicates whether vorticity should be
	///		computed and output.
	///	</summary>
	void OutputVorticity(
		bool fOutputVorticity = true
	) {
		m_fOutputVorticity = fOutputVorticity;
	}

	///	<summary>
	///		Modify the flag which indicates whether divergence should be
	///		computed and output.
	///	</summary>
	void OutputDivergence(
		bool fOutputDivergence = true
	) {
		m_fOutputDivergence = fOutputDivergence;
	}

	///	<summary>
	///		Destructor.
	///	</summary>
	virtual ~OutputManagerReference();

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
	///		Vector of longitudes.
	///	</summary>
	DataVector<double> m_dXCoord;

	///	<summary>
	///		Vector of latitudes.
	///	</summary>
	DataVector<double> m_dYCoord;

	///	<summary>
	///		Vector of alpha reference points.
	///	</summary>
	DataVector<double> m_dAlpha;

	///	<summary>
	///		Vector of reference points.
	///	</summary>
	DataVector<double> m_dBeta;

	///	<summary>
	///		Vector of panel indices.
	///	</summary>
	DataVector<int> m_iPanel;

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
	///		Reference state data.
	///	</summary>
	DataMatrix3D<double> m_dataRefState;

	///	<summary>
	///		Interpolated state data on the reference grid.
	///	</summary>
	DataMatrix3D<double> m_dataState;

	///	<summary>
	///		Interpolated tracers data on the reference grid.
	///	</summary>
	DataMatrix3D<double> m_dataTracers;

private:
	///	<summary>
	///		Flag indicating whether vorticity should be computed and output.
	///	</summary>
	bool m_fOutputVorticity;

	///	<summary>
	///		Vorticity output variable.
	///	</summary>
	NcVar * m_varVorticity;

	///	<summary>
	///		Computed vorticity on the reference grid.
	///	</summary>
	DataMatrix3D<double> m_dataVorticity;

	///	<summary>
	///		Flag indicating whether divergence should be computed and output.
	///	</summary>
	bool m_fOutputDivergence;

	///	<summary>
	///		Divergence output variable.
	///	</summary>
	NcVar * m_varDivergence;

	///	<summary>
	///		Computed divergence on the reference grid.
	///	</summary>
	DataMatrix3D<double> m_dataDivergence;

};

///////////////////////////////////////////////////////////////////////////////

#endif

