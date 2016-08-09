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

#include "DataArray3D.h"

class Time;

///////////////////////////////////////////////////////////////////////////////

#ifndef TEMPEST_NETCDF
typedef int NcFile;
typedef int NcVar;
#endif

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
		const Time & timeOutputFrequency,
		std::string strOutputDir,
		std::string strOutputPrefix,
		int nOutputsPerFile,
		int nXReference,
		int nYReference,
		int nZReference,
		bool fOutputAllVarsOnNodes = true,
		bool fRemoveReferenceProfile = false
	);

	///	<summary>
	///		Destructor.
	///	</summary>
	virtual ~OutputManagerReference();

	///	<summary>
	///		Get the name of the OutputManager.
	///	</summary>
	virtual const char * GetName() const {
		return "Reference";
	}

	///	<summary>
	///		Modify the flag which indicates whether vorticity should be
	///		computed and output.
	///	</summary>
	void OutputVorticity(
		bool fOutputVorticity = true
	);

	///	<summary>
	///		Modify the flag which indicates whether divergence should be
	///		computed and output.
	///	</summary>
	void OutputDivergence(
		bool fOutputDivergence = true
	);

	///	<summary>
	///		Modify the flag which indicates whether temperature should be
	///		computed and output.
	///	</summary>
	void OutputTemperature(
		bool fOutputTemperature = true
	);

	///	<summary>
	///		Modify the flag which indicates whether surface pressure should be
	///		computed and output.
	///	</summary>
	void OutputSurfacePressure(
		bool fOutputSurfacePressure = true
	);

	///	<summary>
	///		Modify the flag which indicates whether Richardson number should be
	///		computed and output.
	///	</summary>
	void OutputRichardson(
		bool fOutputRichardson = true
	);

private:
	///	<summary>
	///		Calculate the patch coordinates of the reference points.
	///	</summary>
	bool CalculatePatchCoordinates();

protected:
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
	virtual void Output(
		const Time & time
	);

protected:
	///	<summary>
	///		Grid stamp.
	///	</summary>
	int m_iGridStamp;

	///	<summary>
	///		Flag indicating a fresh output file (no prior Outputs).
	///	</summary>
	bool m_fFreshOutputFile;

	///	<summary>
	///		Number of reference points in the X direction.
	///	</summary>
	int m_nXReference;

	///	<summary>
	///		Number of reference points in the Y direction.
	///	</summary>
	int m_nYReference;

	///	<summary>
	///		Number of reference points in the Z direction.
	///	</summary>
	int m_nZReference;

	///	<summary>
	///		Flag indicating that all state variables should be output on nodes.
	///	</summary>
	bool m_fOutputAllVarsOnNodes;

	///	<summary>
	///		Flag indicating that the reference profile should not be used.
	///	</summary>
	bool m_fRemoveReferenceProfile;

	///	<summary>
	///		Vector of longitudes.
	///	</summary>
	DataArray1D<double> m_dXCoord;

	///	<summary>
	///		Vector of latitudes.
	///	</summary>
	DataArray1D<double> m_dYCoord;

	///	<summary>
	///		Vector containing a single value (zero) corresponding to REta
	///		at the surface.
	///	</summary>
	DataArray1D<double> m_dREtaSurface;

	///	<summary>
	///		Vector of vertical reference points
	///	</summary>
	DataArray1D<double> m_dREtaCoord;

	///	<summary>
	///		Vector of alpha reference points.
	///	</summary>
	DataArray1D<double> m_dAlpha;

	///	<summary>
	///		Vector of beta reference points.
	///	</summary>
	DataArray1D<double> m_dBeta;

	///	<summary>
	///		Vector of patch indices.
	///	</summary>
	DataArray1D<int> m_iPatch;

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

	///	<summary>
	///		Vector of user data variables.
	///	</summary>
	std::vector<NcVar *> m_vecUserData2DVar;

private:
	///	<summary>
	///		Topography output variable.
	///	</summary>
	NcVar * m_varTopography;

	///	<summary>
	///		Interpolated topography on the reference grid.
	///	</summary>
	DataArray3D<double> m_dataTopography;

private:
	///	<summary>
	///		Interpolated state data on nodes on the reference grid.
	///	</summary>
	DataArray3D<double> m_dataStateNode;

	///	<summary>
	///		Interpolated state data on redges on the reference grid.
	///	</summary>
	DataArray3D<double> m_dataStateREdge;

	///	<summary>
	///		Interpolated tracers data on the reference grid.
	///	</summary>
	DataArray3D<double> m_dataTracers;

	///	<summary>
	///		Interpolated user data on the reference grid.
	///	</summary>
	DataArray3D<double> m_dataUserData2D;

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
	DataArray3D<double> m_dataVorticity;

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
	DataArray3D<double> m_dataDivergence;

	///	<summary>
	///		Flag indicating whether temperature should be computed and output.
	///	</summary>
	bool m_fOutputTemperature;

	///	<summary>
	///		Temperature output variable.
	///	</summary>
	NcVar * m_varTemperature;

	///	<summary>
	///		Computed temperature on the reference grid.
	///	</summary>
	DataArray3D<double> m_dataTemperature;

	///	<summary>
	///		Flag indicating whether surface pressure should be computed
	///		and output.
	///	</summary>
	bool m_fOutputSurfacePressure;

	///	<summary>
	///		Surface pressure output variable.
	///	</summary>
	NcVar * m_varSurfacePressure;

	///	<summary>
	///		Computed surface pressure on the reference grid.
	///	</summary>
	DataArray3D<double> m_dataSurfacePressure;

	///	<summary>
	///		Flag indicating whether Richardson should be computed and output.
	///	</summary>
	bool m_fOutputRichardson;

	///	<summary>
	///		Richardson number output variable.
	///	</summary>
	NcVar * m_varRichardson;

	///	<summary>
	///		Computed Richardson number on the reference grid.
	///	</summary>
	DataArray3D<double> m_dataRichardson;

};

///////////////////////////////////////////////////////////////////////////////

#endif

