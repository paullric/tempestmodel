///////////////////////////////////////////////////////////////////////////////
///
///	\file    ExtractSurface.cpp
///	\author  Paul Ullrich
///	\version July 8, 2014
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

#include "DataVector.h"
#include "DataMatrix.h"
#include "DataMatrix3D.h"
#include "CommandLine.h"
#include "Announce.h"

#include <netcdfcpp.h>

#include "mpi.h"

///////////////////////////////////////////////////////////////////////////////

void InterpolationWeightsLinear(
	double dP,
	const DataVector<double> & dataP,
	int & kBegin,
	int & kEnd,
	DataVector<double> & dW
) {
	const int nLev = dataP.GetRows();

	if (dP > dataP[0]) {
		kBegin = 0;
		kEnd = 1;
		dW[0] = 1.0;

	} else if (dP < dataP[nLev-1]) {
		kBegin = nLev-1;
		kEnd = nLev;
		dW[nLev-1] = 1.0;

	} else {
		for (int k = 0; k < nLev-1; k++) {
			if (dP >= dataP[k+1]) {
				kBegin = k;
				kEnd = k+2;

				dW[k] = (dataP[k+1] - dP)
				      / (dataP[k+1] - dataP[k]);

				dW[k+1] = 1.0 - dW[k];

				break;
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char ** argv) {

	// Initialize MPI
	MPI_Init(&argc, &argv);

try {

	// Input filename
	std::string strInputFile;

	// Output filename
	std::string strOutputFile;

	// Separate topography file
	std::string strTopographyFile;

	// List of variables to extract
	std::string strVariables;

	// Extract geopotential height
	bool fGeopotentialHeight;

	// Pressure level to extract
	double dPressureLevel;

	// Extract variables at the surface
	bool fExtractSurface;

	// Time index
	int iTime;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in", "");
		CommandLineString(strOutputFile, "out", "");
		CommandLineString(strVariables, "var", "");
		CommandLineBool(fGeopotentialHeight, "output_z");
		CommandLineDouble(dPressureLevel, "p", 0.0);
		CommandLineBool(fExtractSurface, "surf");
		CommandLineInt(iTime, "time", 0);

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Check command line arguments
	if (strInputFile == "") {
		_EXCEPTIONT("No input file specified");
	}
	if (strOutputFile == "") {
		_EXCEPTIONT("No output file specified");
	}
	if (strVariables == "") {
		_EXCEPTIONT("No variables specified");
	}
	if ((dPressureLevel == 0.0) && (!fExtractSurface)) {
		_EXCEPTIONT("No instruction specified");
	}
	if ((dPressureLevel != 0.0) && (fExtractSurface)) {
		_EXCEPTIONT("Only one instruction allowed");
	}

	// Open input file
	AnnounceStartBlock("Loading input file");
	NcFile ncdf_in(strInputFile.c_str(), NcFile::ReadOnly);
	if (!ncdf_in.is_valid()) {
		_EXCEPTION1("Unable to open file \"%s\" for reading",
			strInputFile.c_str());
	}

	// Load time array
	Announce("Time");
	NcVar * varTime = ncdf_in.get_var("time");
	int nTime = varTime->get_dim(0)->size();

	DataVector<double> dTime;
	dTime.Initialize(nTime);
	varTime->set_cur((long)0);
	varTime->get(&(dTime[0]), nTime);

	if (iTime >= nTime) {
		_EXCEPTION2("iTime index (%i) out of range [0,%i]",
			iTime, nTime-1);
	}

	// Load latitude array
	Announce("Latitude");
	NcVar * varLat = ncdf_in.get_var("lat");
	int nLat = varLat->get_dim(0)->size();

	DataVector<double> dLat;
	dLat.Initialize(nLat);
	varLat->set_cur((long)0);
	varLat->get(&(dLat[0]), nLat);

	// Load longitude array
	Announce("Longitude");
	NcVar * varLon = ncdf_in.get_var("lon");
	int nLon = varLon->get_dim(0)->size();

	DataVector<double> dLon;
	dLon.Initialize(nLon);
	varLon->set_cur((long)0);
	varLon->get(&(dLon[0]), nLon);

	// Load level array
	Announce("Level");
	NcVar * varLev = ncdf_in.get_var("lev");
	int nLev = varLev->get_dim(0)->size();

	DataVector<double> dLev;
	dLev.Initialize(nLev);
	varLev->set_cur((long)0);
	varLev->get(&(dLev[0]), nLev);

	AnnounceEndBlock("Done");

	// Load topography
	NcVar *varZs = ncdf_in.get_var("Zs");

	DataMatrix<double> dZs;
	dZs.Initialize(nLat, nLon);
	varZs->set_cur((long)0, (long)0);
	varZs->get(&(dZs[0][0]), nLat, nLon);

	// Open output file
	AnnounceStartBlock("Constructing output file");

	NcFile ncdf_out(strOutputFile.c_str(), NcFile::Replace);
	if (!ncdf_out.is_valid()) {
		_EXCEPTION1("Unable to open file \"%s\" for writing",
			strOutputFile.c_str());
	}

	// Output latitude and longitude array
	Announce("Latitude");
	NcDim * dimOutLat = ncdf_out.add_dim("lat", nLat);
	NcVar * varOutLat = ncdf_out.add_var("lat", ncDouble, dimOutLat);
	varOutLat->set_cur((long)0);
	varOutLat->put(&(dLat[0]), nLat);

	Announce("Longitude");
	NcDim * dimOutLon = ncdf_out.add_dim("lon", nLon);
	NcVar * varOutLon = ncdf_out.add_var("lon", ncDouble, dimOutLon);
	varOutLon->set_cur((long)0);
	varOutLon->put(&(dLon[0]), nLon);

	// Output pressure array
	Announce("Pressure");
	NcDim * dimOutP = ncdf_out.add_dim("p", 1);
	NcVar * varOutP = ncdf_out.add_var("p", ncDouble, dimOutP);
	varOutP->set_cur((long)0);
	varOutP->put(&dPressureLevel, 1);

	AnnounceEndBlock("Done");

	// Parse variable string
	Announce("Loading variables");
	int iVarBegin = 0;
	int iVarCurrent = 0;

	std::vector<std::string> vecVariableStrings;

	// Parse variable name
	for (;;) {
		if ((iVarCurrent >= strVariables.length()) ||
			(strVariables[iVarCurrent] == ',') ||
			(strVariables[iVarCurrent] == ' ')
		) {
			if (iVarCurrent == iVarBegin) {
				if (iVarCurrent >= strVariables.length()) {
					break;
				}

				continue;
			}

			vecVariableStrings.push_back(
				strVariables.substr(iVarBegin, iVarCurrent - iVarBegin));

			iVarBegin = iVarCurrent + 1;
		}

		iVarCurrent++;
	}

	// Check variables
	if (vecVariableStrings.size() == 0) {
		_EXCEPTIONT("No variables specified");
	}

	// Load all variables
	std::vector<NcVar *> vecNcVar;
	for (int v = 0; v < vecVariableStrings.size(); v++) {
		vecNcVar.push_back(ncdf_in.get_var(vecVariableStrings[v].c_str()));
		if (vecNcVar[v] == NULL) {
			_EXCEPTION1("Unable to load variable \"%s\" from file",
				vecVariableStrings[v].c_str());
		}
	}

	// Physical constants
	Announce("Initializing thermodynamic variables");
	NcAtt * attRd = ncdf_in.get_att("Rd");
	double dRd = attRd->as_double(0);

	NcAtt * attCp = ncdf_in.get_att("Cp");
	double dCp = attCp->as_double(0);

	double dGamma = dCp / (dCp - dRd);

	NcAtt * attP0 = ncdf_in.get_att("P0");
	double dP0 = attP0->as_double(0);

	double dPressureScaling = dP0 * pow(dRd / dP0, dGamma);

	NcAtt * attZtop = ncdf_in.get_att("Ztop");
	double dZtop = attZtop->as_double(0);

	// Rho
	DataMatrix3D<double> dataRho;
	dataRho.Initialize(nLev, nLat, nLon);

	NcVar * varRho = ncdf_in.get_var("Rho");
	varRho->set_cur(iTime, 0, 0, 0);
	varRho->get(&(dataRho[0][0][0]), 1, nLev, nLat, nLon);

	// Theta
	DataMatrix3D<double> dataTheta;
	dataTheta.Initialize(nLev, nLat, nLon);

	NcVar * varTheta = ncdf_in.get_var("Theta");
	varTheta->set_cur(iTime, 0, 0, 0);
	varTheta->get(&(dataTheta[0][0][0]), 1, nLev, nLat, nLon);

	// Pressure everywhere
	DataMatrix3D<double> dataP;
	dataP.Initialize(nLev, nLat, nLon);

	for (int k = 0; k < nLev; k++) {
	for (int i = 0; i < nLat; i++) {
	for (int j = 0; j < nLon; j++) {
		dataP[k][i][j] = dPressureScaling
			* exp(log(dataRho[k][i][j] * dataTheta[k][i][j]) * dGamma);
	}
	}
	}

	// Input data
	DataMatrix3D<double> dataIn;
	dataIn.Initialize(nLev, nLat, nLon);

	// Output data
	DataMatrix<double> dataOut;
	dataOut.Initialize(nLat, nLon);

	// Pressure in column
	DataVector<double> dataColumnP;
	dataColumnP.Initialize(nLev);

	// Column weights
	DataVector<double> dW;
	dW.Initialize(nLev);

	// Loop through all variables
	AnnounceStartBlock("Interpolating");

	for (int v = 0; v < vecNcVar.size(); v++) {

		Announce(vecVariableStrings[v].c_str());

		// Create new variable in output array
		NcVar * varOut = ncdf_out.add_var(
			vecVariableStrings[v].c_str(), ncDouble, dimOutLat, dimOutLon);

		// Load in the data array
		vecNcVar[v]->set_cur(iTime, 0, 0, 0);
		vecNcVar[v]->get(&(dataIn[0][0][0]), 1, nLev, nLat, nLon);

		// Loop thorugh all latlon indices
		for (int i = 0; i < nLat; i++) {
		for (int j = 0; j < nLon; j++) {

			// Find weights
			int kBegin = 0;
			int kEnd = 0;

			// On a pressure surface
			if (dPressureLevel != 0.0) {
				for (int k = 0; k < nLev; k++) {
					dataColumnP[k] = dataP[k][i][j];
				}

				InterpolationWeightsLinear(
					dPressureLevel,
					dataColumnP,
					kBegin,
					kEnd,
					dW);
			}

			// At the physical surface
			if (fExtractSurface) {
				kBegin = 0;
				kEnd = 2;

				dW[0] =   dLev[1] / (dLev[1] - dLev[0]);
				dW[1] = - dLev[0] / (dLev[1] - dLev[0]);
			}

/*
			// DEBUG
			double dTestRho = 0.0;
			double dTestTheta = 0.0;
			double dP;
			for (int k = kBegin; k < kEnd; k++) {
				dTestRho += dW[k] * dataRho[k][i][j];
				dTestTheta += dW[k] * dataTheta[k][i][j];
			}

			dP = dPressureScaling
				* exp(log(dTestRho * dTestTheta) * dGamma);

			printf("%1.5e %1.5e %1.5e\n", dataP[kBegin], dP, dataP[kEnd-1]);
*/
			// Interpolate in the vertical
			dataOut[i][j] = 0.0;
			for (int k = kBegin; k < kEnd; k++) {
				dataOut[i][j] += dW[k] * dataIn[k][i][j];
			}
		}
		}

		// Write variable
		varOut->set_cur(0, 0);
		varOut->put(&(dataOut[0][0]), nLat, nLon);
	}

	// Output geopotential height
	if (fGeopotentialHeight) {

		Announce("Geopotential height");

		// Loop thorugh all latlon indices
		for (int i = 0; i < nLat; i++) {
		for (int j = 0; j < nLon; j++) {

			int kBegin;
			int kEnd;

			// Interpolate onto pressure level
			if (dPressureLevel != 0.0) {
				for (int k = 0; k < nLev; k++) {
					dataColumnP[k] = dataP[k][i][j];
				}

				InterpolationWeightsLinear(
					dPressureLevel,
					dataColumnP,
					kBegin,
					kEnd,
					dW);
			}

			// At the physical surface
			if (fExtractSurface) {
				kBegin = 0;
				kEnd = 2;

				dW[0] =   dLev[1] / (dLev[1] - dLev[0]);
				dW[1] = - dLev[0] / (dLev[1] - dLev[0]);
			}

			dataOut[i][j] = 0.0;
			for (int k = kBegin; k < kEnd; k++) {
				dataOut[i][j] += dW[k] * dLev[k];
			}

			dataOut[i][j] = dZs[i][j] + dataOut[i][j] * (dZtop - dZs[i][j]);
		}
		}

		// Write variable
		NcVar * varOut = ncdf_out.add_var(
			"Z", ncDouble, dimOutLat, dimOutLon);

		varOut->set_cur(0, 0);
		varOut->put(&(dataOut[0][0]), nLat, nLon);
	}

	AnnounceEndBlock("Done");

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}

	// Finalize MPI
	MPI_Finalize();

}

///////////////////////////////////////////////////////////////////////////////

