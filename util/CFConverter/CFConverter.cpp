///////////////////////////////////////////////////////////////////////////////
///
///	\file    CFConverter.cpp
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

#include "DataArray1D.h"
#include "DataArray2D.h"
#include "DataArray3D.h"
#include "CommandLine.h"
#include "Announce.h"

#include <cstdlib>
#include <cmath>
#include <netcdfcpp.h>

#include "PolynomialInterp.h"

#include <mpi.h>

////////////////////////////////////////////////////////////////////////////////

void CopyNcFileAttributes(
	NcFile * fileIn,
	NcFile * fileOut
) {
	for (int a = 0; a < fileIn->num_atts(); a++) {
		NcAtt * att = fileIn->get_att(a);
		long num_vals = att->num_vals();

		NcValues * pValues = att->values();

		if (att->type() == ncByte) {
			fileOut->add_att(att->name(), num_vals,
				(const ncbyte*)(pValues->base()));

		} else if (att->type() == ncChar) {
			fileOut->add_att(att->name(), num_vals,
				(const char*)(pValues->base()));

		} else if (att->type() == ncShort) {
			fileOut->add_att(att->name(), num_vals,
				(const short*)(pValues->base()));

		} else if (att->type() == ncInt) {
			fileOut->add_att(att->name(), num_vals,
				(const int*)(pValues->base()));

		} else if (att->type() == ncFloat) {
			fileOut->add_att(att->name(), num_vals,
				(const float*)(pValues->base()));

		} else if (att->type() == ncDouble) {
			fileOut->add_att(att->name(), num_vals,
				(const double*)(pValues->base()));

		} else {
			_EXCEPTIONT("Invalid attribute type");
		}

		delete pValues;
	}
}

////////////////////////////////////////////////////////////////////////////////

void CopyNcVarAttributes(
	NcVar * varIn,
	NcVar * varOut
) {
	for (int a = 0; a < varIn->num_atts(); a++) {
		NcAtt * att = varIn->get_att(a);
		long num_vals = att->num_vals();

		NcValues * pValues = att->values();

		if (att->type() == ncByte) {
			varOut->add_att(att->name(), num_vals,
				(const ncbyte*)(pValues->base()));

		} else if (att->type() == ncChar) {
			varOut->add_att(att->name(), num_vals,
				(const char*)(pValues->base()));

		} else if (att->type() == ncShort) {
			varOut->add_att(att->name(), num_vals,
				(const short*)(pValues->base()));

		} else if (att->type() == ncInt) {
			varOut->add_att(att->name(), num_vals,
				(const int*)(pValues->base()));

		} else if (att->type() == ncFloat) {
			varOut->add_att(att->name(), num_vals,
				(const float*)(pValues->base()));

		} else if (att->type() == ncDouble) {
			varOut->add_att(att->name(), num_vals,
				(const double*)(pValues->base()));

		} else {
			_EXCEPTIONT("Invalid attribute type");
		}

		delete pValues;
	}
}

///////////////////////////////////////////////////////////////////////////////

void ParseLevelArray(
	const std::string & strLevels,
	std::vector<double> & vecLevels
) {
	int iLevelBegin = 0;
	int iLevelCurrent = 0;

	vecLevels.clear();

	if (strLevels == "") {
		return;
	}

	// Parse pressure levels
	bool fRangeMode = false;
	for (;;) {
		if ((iLevelCurrent >= strLevels.length()) ||
			(strLevels[iLevelCurrent] == ',') ||
			(strLevels[iLevelCurrent] == ' ') ||
			(strLevels[iLevelCurrent] == ':')
		) {
			// Range mode
			if ((!fRangeMode) &&
				(strLevels[iLevelCurrent] == ':')
			) {
				if (vecLevels.size() != 0) {
					_EXCEPTIONT("Invalid set of pressure levels");
				}
				fRangeMode = true;
			}
			if (fRangeMode) {
				if ((strLevels[iLevelCurrent] != ':') &&
					(iLevelCurrent < strLevels.length())
				) {
					_EXCEPTION1("Invalid character in pressure range (%c)",
						strLevels[iLevelCurrent]);
				}
			}

			if (iLevelCurrent == iLevelBegin) {
				if (iLevelCurrent >= strLevels.length()) {
					break;
				}

				continue;
			}

			std::string strPressureLevelSubStr = 
				strLevels.substr(
					iLevelBegin, iLevelCurrent - iLevelBegin);

			vecLevels.push_back(atof(strPressureLevelSubStr.c_str()));
			
			iLevelBegin = iLevelCurrent + 1;
		}

		iLevelCurrent++;
	}

	// Range mode -- repopulate array
	if (fRangeMode) {
		if (vecLevels.size() != 3) {
			_EXCEPTIONT("Exactly three pressure level entries required "
				"for range mode");
		}
		double dLevelBegin = vecLevels[0];
		double dLevelStep = vecLevels[1];
		double dLevelEnd = vecLevels[2];

		if (dLevelStep == 0.0) {
			_EXCEPTIONT("Level step size cannot be zero");
		}
		if ((dLevelEnd - dLevelBegin) / dLevelStep > 10000.0) {
			_EXCEPTIONT("Too many levels in range (limit 10000)");
		}
		if ((dLevelEnd - dLevelBegin) / dLevelStep < 0.0) {
			_EXCEPTIONT("Sign mismatch in level step");
		}

		vecLevels.clear();
		for (int i = 0 ;; i++) {
			double dLevel = dLevelBegin + static_cast<double>(i) * dLevelStep;

			if ((dLevelStep > 0.0) && (dLevel > dLevelEnd)) {
				break;
			}
			if ((dLevelStep < 0.0) && (dLevel < dLevelEnd)) {
				break;
			}

			vecLevels.push_back(dLevel);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void ParseVariableList(
	const std::string & strVariables,
	std::vector< std::string > & vecVariableStrings
) {
	int iVarBegin = 0;
	int iVarCurrent = 0;

	// Parse variable name
	for (;;) {
		if ((iVarCurrent >= strVariables.length()) ||
			(strVariables[iVarCurrent] == ',') ||
			(strVariables[iVarCurrent] == '.') ||
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
}

///////////////////////////////////////////////////////////////////////////////

void InterpolationWeightsLinear(
	double dP,
	const DataArray1D<double> & dataP,
	int & kBegin,
	int & kEnd,
	DataArray1D<double> & dW
) {
	const int nLev = dataP.GetRows();

	// Monotone increasing coordinate
	if (dataP[1] > dataP[0]) {
		if (dP < dataP[0]) {
			kBegin = 0;
			kEnd = 2;

		} else if (dP > dataP[nLev-1]) {
			kBegin = nLev-2;
			kEnd = nLev;

		} else {
			for (int k = 0; k < nLev-1; k++) {
				if (dP <= dataP[k+1]) {
					kBegin = k;
					kEnd = k+2;

					break;
				}
			}
		}

	// Monotone decreasing coordinate
	} else {
		if (dP > dataP[0]) {
			kBegin = 0;
			kEnd = 2;

		} else if (dP < dataP[nLev-1]) {
			kBegin = nLev-2;
			kEnd = nLev;

		} else {
			for (int k = 0; k < nLev-1; k++) {
				if (dP >= dataP[k+1]) {
					kBegin = k;
					kEnd = k+2;

					break;
				}
			}
		}
	}

	// Weights
	dW[kBegin] =
		  (dataP[kBegin+1] - dP)
	    / (dataP[kBegin+1] - dataP[kBegin]);

	dW[kBegin+1] = 1.0 - dW[kBegin];
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char ** argv) {

	MPI_Init(&argc, &argv);

	NcError error(NcError::silent_nonfatal);

try {

	// Input filename
	std::string strInputFile;

	// List of variables to extract
	std::string strVariables;

	// Project Id
	std::string strProjectId;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in", "");
		CommandLineString(strVariables, "var", "");
		CommandLineString(strProjectId, "project_id", "DCMIP2016");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Check command line arguments
	if (strInputFile == "") {
		_EXCEPTIONT("No input file specified");
	}
	if (strVariables == "") {
		_EXCEPTIONT("No variables specified");
	}

	// Parse input file string
	std::vector< std::string > vecConfigStrings;

	ParseVariableList(strInputFile, vecConfigStrings);

	if ((vecConfigStrings.size() < 7) ||
	    (vecConfigStrings.size() > 8) ||
		(vecConfigStrings[vecConfigStrings.size()-1] != "nc")
	) {
		_EXCEPTIONT("\nInvalid input filename.  Expected:\n"
			"model.experiment_id.resolution.levels.grid.equation[.description].nc\n"
			"with optional [.description]");
	}

	// Parse variable string
	std::vector< std::string > vecVariableStrings;

	ParseVariableList(strVariables, vecVariableStrings);

	// Check variables
	if (vecVariableStrings.size() == 0) {
		_EXCEPTIONT("No variables specified");
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
	if (varTime == NULL) {
		_EXCEPTION1("File \"%s\" does not contain variable \"time\"",
			strInputFile.c_str());
	}
	int nTime = varTime->get_dim(0)->size();

	DataArray1D<double> dTime(nTime);
	varTime->set_cur((long)0);
	varTime->get(&(dTime[0]), nTime);

	// Load latitude array
	Announce("Latitude");
	NcVar * varLat = ncdf_in.get_var("lat");
	if (varLat == NULL) {
		_EXCEPTION1("File \"%s\" does not contain variable \"lat\"",
			strInputFile.c_str());
	}
	int nLat = varLat->get_dim(0)->size();

	DataArray1D<double> dLat(nLat);
	varLat->set_cur((long)0);
	varLat->get(&(dLat[0]), nLat);

	// Load longitude array
	Announce("Longitude");
	NcVar * varLon = ncdf_in.get_var("lon");
	if (varLon == NULL) {
		_EXCEPTION1("File \"%s\" does not contain variable \"lon\"",
			strInputFile.c_str());
	}
	int nLon = varLon->get_dim(0)->size();

	DataArray1D<double> dLon(nLon);
	varLon->set_cur((long)0);
	varLon->get(&(dLon[0]), nLon);

	// Load Ztop parameter
	NcAtt * attZtop = ncdf_in.get_att("Ztop");
	if (attZtop == NULL) {
		_EXCEPTIONT("Attribute \"Ztop\" not found");
	}
	double dZtop = attZtop->as_double(0);

	// Load level array
	Announce("Level");
	NcDim * dimLev = ncdf_in.get_dim("lev");
	if (dimLev == NULL) {
		_EXCEPTIONT("Dimension \"lev\" not found");
	}
	NcVar * varLev = ncdf_in.get_var("lev");
	if (varLev == NULL) {
		_EXCEPTIONT("Variable \"lev\" not found");
	}
	int nLev = dimLev->size();

	DataArray1D<double> dLev(nLev);
	varLev->set_cur((long)0);
	varLev->get(&(dLev[0]), nLev);

	for (int k = 0; k < nLev; k++) {
		dLev[k] *= dZtop;
	}

	// Load level interface array
	Announce("Interface");
	NcDim * dimILev = ncdf_in.get_dim("ilev");
	if (dimILev == NULL) {
		_EXCEPTIONT("Dimension \"ilev\" not found");
	}
	NcVar * varILev = ncdf_in.get_var("ilev");
	int nILev = 0;
	DataArray1D<double> dILev;
	if (varILev == NULL) {
		_EXCEPTIONT("Variable \"ilev\" not found");
	} else {
		nILev = varILev->get_dim(0)->size();
		if (nILev != nLev + 1) {
			_EXCEPTIONT("Variable \"ilev\" must have size lev+1");
		}
		dILev.Allocate(nILev);
		varILev->set_cur((long)0);
		varILev->get(&(dILev[0]), nILev);
	}

	for (int k = 0; k < nILev; k++) {
		dILev[k] *= dZtop;
	}
/*
	// Load topography
	Announce("Topography");
	NcVar * varZs = ncdf_in.get_var("Zs");
	if (varZs == NULL) {
		_EXCEPTION1("File \"%s\" does not contain variable \"Zs\"",
			strInputFile.c_str());
	}

	DataArray2D<double> dZs(nLat, nLon);
	varZs->set_cur((long)0, (long)0);
	varZs->get(&(dZs[0][0]), nLat, nLon);

	AnnounceEndBlock("Done");
*/
	// Load density
	Announce("Density");
	bool fHasDensity = false;
	NcVar * varRho = ncdf_in.get_var("Rho");
	DataArray1D<double> dRho;
	if (varRho != NULL) {
		fHasDensity = true;
		if (varRho->get_dim(0)->size() != nTime) {
			_EXCEPTION2("Invalid dim 0 size on Rho: %i %i",
				varRho->get_dim(0)->size(), nTime);
		}
		if (varRho->get_dim(1)->size() != nLev) {
			_EXCEPTION2("Invalid dim 1 size on Rho: %i %i",
				varRho->get_dim(1)->size(), nLev);
		}
		if (varRho->get_dim(2)->size() != nLat) {
			_EXCEPTION2("Invalid dim 2 size on Rho: %i %i",
				varRho->get_dim(2)->size(), nLat);
		}
		if (varRho->get_dim(3)->size() != nLon) {
			_EXCEPTION2("Invalid dim 3 size on Rho: %i %i",
				varRho->get_dim(3)->size(), nLon);
		}

		dRho.Allocate(nTime * nLev * nLon * nLat);
		varRho->set_cur(0, 0, 0, 0);
		varRho->get(&(dRho[0]), nTime, nLev, nLat, nLon);

		for (int i = 0; i < dRho.GetRows(); i++) {
			if (dRho[i] == 0.0) {
				_EXCEPTIONT("Zero density detected");
			}
		}
	}

	// Open output file
	AnnounceStartBlock("Constructing output files");

	// Loop through all variables
	for (int v = 0; v < vecVariableStrings.size(); v++) {

		// Start block
		AnnounceStartBlock(vecVariableStrings[v].c_str());

		// Input variable
		NcVar * varIn = NULL;
		int nDim = 0;
		DataArray1D<double> data;

		bool f2D = false;
		bool fOnLevels = false;
		bool fOnInterfaces = false;

		// Tracers
		if (vecVariableStrings[v][0] == 'Q') {
			std::string strQ = vecVariableStrings[v];
			if (strQ == "Q") {
				strQ = "Qv";
			}
			varIn = ncdf_in.get_var(strQ.c_str());
			if (varIn == NULL) {
				std::string strRhoQ = "Rho" + strQ;
				if (!fHasDensity) {
					_EXCEPTION1("Cannot find variable \"%s\" and density "
						"not available", strQ.c_str());
				}
				varIn = ncdf_in.get_var(strRhoQ.c_str());
				if (varIn == NULL) {
					_EXCEPTION2("Cannot find variable \"%s\" or \"%s\"",
						strQ.c_str(), strRhoQ.c_str());
				}
			}
			nDim = varIn->num_dims();
			if (nDim != 4) {
				_EXCEPTION1("Invalid number of dimensions for variable \"%s\"",
					strQ.c_str());
			}
			data.Allocate(nTime * nLev * nLat * nLon);
			varIn->get(&(data[0]), nTime, nLev, nLat, nLon);

			for (int i = 0; i < data.GetRows(); i++) {
				data[i] /= dRho[i];
			}

			fOnLevels = true;

		// Non-tracers
		} else {
			varIn = ncdf_in.get_var(vecVariableStrings[v].c_str());
			if (varIn == NULL) {
				_EXCEPTION1("Unable to load variable \"%s\"",
					vecVariableStrings[v].c_str());
			}

			// Number of dimensions
			nDim = varIn->num_dims();
			if (nDim == 3) {
				f2D = true;
			} else if (nDim == 4) {
				NcDim * dimLev = varIn->get_dim(1);
				if (strcmp(dimLev->name(), "lev") == 0) {
					fOnLevels = true;
				}
				if (strcmp(dimLev->name(), "ilev") == 0) {
					fOnInterfaces = true;
				}
			} else {
				_EXCEPTION1("Invalid number of dimensions for variable \"%s\"",
					vecVariableStrings[v].c_str());
			}

			if (f2D) {
				data.Allocate(nTime * nLat * nLon);
				varIn->get(&(data[0]), nTime, nLat, nLon);
			}
			if (fOnLevels) {
				data.Allocate(nTime * nLev * nLat * nLon);
				varIn->get(&(data[0]), nTime, nLev, nLat, nLon);
			}
			if (fOnInterfaces) {
				data.Allocate(nTime * nILev * nLat * nLon);
				varIn->get(&(data[0]), nTime, nILev, nLat, nLon);
			}
		}

		if (varIn == NULL) {
			_EXCEPTIONT("Logic error");
		}

		// Open output file
		std::string strOutputFile;
		for (int k = 0; k < vecConfigStrings.size()-1; k++) {
			strOutputFile += vecConfigStrings[k];
			strOutputFile += ".";
		}
		strOutputFile += vecVariableStrings[v];
		strOutputFile += ".nc";

		NcFile ncdf_out(strOutputFile.c_str(), NcFile::Replace);
		if (!ncdf_out.is_valid()) {
			_EXCEPTION1("Unable to open file \"%s\" for writing",
				strOutputFile.c_str());
		}

		CopyNcFileAttributes(&ncdf_in, &ncdf_out);

		// Add ESMF indexing attributes
		{
			ncdf_out.add_att("Conventions", "CF-1.0");
			ncdf_out.add_att("project_id", strProjectId.c_str());
			ncdf_out.add_att("institute_id", "UCD/LBNL");
			ncdf_out.add_att("experiment_id", vecConfigStrings[1].c_str());
			ncdf_out.add_att("model_id", vecConfigStrings[0].c_str());
			ncdf_out.add_att("frequency", "unknown");
			ncdf_out.add_att("modeling_realm", "atmos");
			ncdf_out.add_att("horizontal_resolution", vecConfigStrings[2].c_str());
			ncdf_out.add_att("levels", vecConfigStrings[3].c_str());
			ncdf_out.add_att("grid", vecConfigStrings[4].c_str());
			ncdf_out.add_att("equation", vecConfigStrings[5].c_str());

			if (vecConfigStrings.size() == 8) {
				ncdf_out.add_att("description", vecConfigStrings[6].c_str());
			}
		}

		// Output time array
		Announce("Time");
		NcDim * dimOutTime = ncdf_out.add_dim("time");
		NcVar * varOutTime = ncdf_out.add_var("time", ncDouble, dimOutTime);
		varOutTime->set_cur((long)0);
		varOutTime->put(&(dTime[0]), nTime);

		CopyNcVarAttributes(varTime, varOutTime);

		// Output height array
		NcDim * dimOutZ = NULL;
		NcVar * varOutZ = NULL;
		if (fOnLevels) {
			Announce("Altitude");
			dimOutZ = ncdf_out.add_dim("z", nLev);
			varOutZ = ncdf_out.add_var("z", ncDouble, dimOutZ);
			varOutZ->set_cur((long)0);
			varOutZ->put(&(dLev[0]), nLev);
			varOutZ->add_att("long_name", "altitude");
			varOutZ->add_att("units", "m");
		}
		if (fOnInterfaces) {
			Announce("Altitude");
			dimOutZ = ncdf_out.add_dim("z", nILev);
			varOutZ = ncdf_out.add_var("z", ncDouble, dimOutZ);
			varOutZ->set_cur((long)0);
			varOutZ->put(&(dILev[0]), nILev);
			varOutZ->add_att("long_name", "altitude");
			varOutZ->add_att("units", "m");
		}

		// Output latitude and longitude array
		Announce("Latitude");
		NcDim * dimOutLat = ncdf_out.add_dim("lat", nLat);
		NcVar * varOutLat = ncdf_out.add_var("lat", ncDouble, dimOutLat);
		varOutLat->set_cur((long)0);
		varOutLat->put(&(dLat[0]), nLat);

		CopyNcVarAttributes(varLat, varOutLat);

		Announce("Longitude");
		NcDim * dimOutLon = ncdf_out.add_dim("lon", nLon);
		NcVar * varOutLon = ncdf_out.add_var("lon", ncDouble, dimOutLon);
		varOutLon->set_cur((long)0);
		varOutLon->put(&(dLon[0]), nLon);

		CopyNcVarAttributes(varLon, varOutLon);

		// Copy over data
		Announce("Variable data");
		NcVar * varOut = NULL;
		if (f2D) {
	   		varOut =
				ncdf_out.add_var(
					vecVariableStrings[v].c_str(),
					ncDouble,
					dimOutTime,
					dimOutLat,
					dimOutLon);

			varOut->put(&(data[0]), nTime, nLat, nLon);
		}
		if (fOnLevels) {
	   		varOut =
				ncdf_out.add_var(
					vecVariableStrings[v].c_str(),
					ncDouble,
					dimOutTime,
					dimOutZ,
					dimOutLat,
					dimOutLon);

			varOut->put(&(data[0]), nTime, nLev, nLat, nLon);
		}
		if (fOnInterfaces) {
	   		varOut =
				ncdf_out.add_var(
					vecVariableStrings[v].c_str(),
					ncDouble,
					dimOutTime,
					dimOutZ,
					dimOutLat,
					dimOutLon);

			varOut->put(&(data[0]), nTime, nILev, nLat, nLon);
		}
		if (varOut == NULL) {
			_EXCEPTIONT("Logic error");
		}

		// Add long names
		if (vecVariableStrings[v] == "U") {
			varOut->add_att("units", "m/s");
			varOut->add_att("long_name", "eastward_wind");

		} else if (vecVariableStrings[v] == "V") {
			varOut->add_att("units", "m/s");
			varOut->add_att("long_name", "northward_wind");

		} else if (vecVariableStrings[v] == "W") {
			varOut->add_att("units", "m/s");
			varOut->add_att("long_name", "upward_air_velocity");

		} else if (vecVariableStrings[v] == "Rho") {
			varOut->add_att("units", "kg/m^3");
			varOut->add_att("long_name", "air_density");

		} else if (vecVariableStrings[v] == "Theta") {
			varOut->add_att("units", "K");
			varOut->add_att("long_name", "potential_temperature");

		} else if (vecVariableStrings[v] == "ThetaV") {
			varOut->add_att("units", "K");
			varOut->add_att("long_name", "virtual_potential_temperature");

		} else if (vecVariableStrings[v] == "T") {
			varOut->add_att("units", "K");
			varOut->add_att("long_name", "temperature");

		} else if (vecVariableStrings[v] == "TV") {
			varOut->add_att("units", "K");
			varOut->add_att("long_name", "virtual_temperature");

		} else if (vecVariableStrings[v] == "DELTA") {
			varOut->add_att("units", "1/s");
			varOut->add_att("long_name", "divergence");

		} else if (vecVariableStrings[v] == "ZETA") {
			varOut->add_att("units", "1/s");
			varOut->add_att("long_name", "relative_vorticity");

		} else if (vecVariableStrings[v] == "Q") {
			varOut->add_att("units", "kg/m^3");
			varOut->add_att("long_name", "specific_humidity");

		} else if (vecVariableStrings[v] == "Qc") {
			varOut->add_att("units", "kg/m^3");
			varOut->add_att("long_name", "cloud_water_mixing_ratio");

		} else if (vecVariableStrings[v] == "Qr") {
			varOut->add_att("units", "kg/m^3");
			varOut->add_att("long_name", "rain_water_mixing_ratio");

		} else if (vecVariableStrings[v] == "QCl") {
			varOut->add_att("units", "kg/m^3");
			varOut->add_att("long_name", "singlet_chlorine_mixing_ratio");

		} else if (vecVariableStrings[v] == "QCl2") {
			varOut->add_att("units", "kg/m^3");
			varOut->add_att("long_name", "chlorine_gas_mixing_ratio");

		} else if (vecVariableStrings[v] == "PS") {
			varOut->add_att("units", "Pa");
			varOut->add_att("long_name", "surface_pressure");

		} else {
			printf("WARNING: No CF-Compliant version of variable \"%s\" known",
				vecVariableStrings[v].c_str());
		}

		// Done with variable
		AnnounceEndBlock("Done");
	}

	AnnounceEndBlock("Done");

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}

	// Finalize MPI
	MPI_Finalize();
}

///////////////////////////////////////////////////////////////////////////////

