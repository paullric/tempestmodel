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

	// Output filename
	std::string strOutputFile;

	// Separate topography file
	std::string strTopographyFile;

	// List of variables to extract
	std::string strVariables;

	// Extract geopotential height
	bool fGeopotentialHeight;

	// Pressure levels to extract
	std::string strPressureLevels;

	// Height levels to extract
	std::string strHeightLevels;

	// Extract variables at the surface
	bool fExtractSurface;

	// Extract total energy
	bool fExtractTotalEnergy;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in", "");
		CommandLineString(strOutputFile, "out", "");
		CommandLineString(strVariables, "var", "");
		CommandLineBool(fGeopotentialHeight, "output_z");
		CommandLineBool(fExtractTotalEnergy, "output_energy");
		CommandLineString(strPressureLevels, "p", "");
		CommandLineString(strHeightLevels, "z", "");
		CommandLineBool(fExtractSurface, "surf");

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

	// Parse variable string
	std::vector< std::string > vecVariableStrings;

	ParseVariableList(strVariables, vecVariableStrings);

	// Check variables
	if (vecVariableStrings.size() == 0) {
		_EXCEPTIONT("No variables specified");
	}

	// Parse pressure level string
	std::vector<double> vecPressureLevels;

	ParseLevelArray(strPressureLevels, vecPressureLevels);

	int nPressureLevels = (int)(vecPressureLevels.size());

	for (int k = 0; k < nPressureLevels; k++) {
		if (vecPressureLevels[k] <= 0.0) {
			_EXCEPTIONT("Non-positive pressure values not allowed");
		}
	}

	// Parse height level string
	std::vector<double> vecHeightLevels;

	ParseLevelArray(strHeightLevels, vecHeightLevels);

	int nHeightLevels = (int)(vecHeightLevels.size());

	// Check pressure levels
	if ((nPressureLevels == 0) &&
		(nHeightLevels == 0) &&
		(!fExtractSurface)
	) {
		_EXCEPTIONT("No pressure / height levels to process");
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

	// Load level array
	Announce("Level");
	NcVar * varLev = ncdf_in.get_var("lev");
	if (varLev == NULL) {
		_EXCEPTION1("File \"%s\" does not contain variable \"lev\"",
			strInputFile.c_str());
	}
	int nLev = varLev->get_dim(0)->size();

	DataArray1D<double> dLev(nLev);
	varLev->set_cur((long)0);
	varLev->get(&(dLev[0]), nLev);

	// Load level interface array
	Announce("Interface");
	NcVar * varILev = ncdf_in.get_var("ilev");
	int nILev = 0;
	DataArray1D<double> dILev;
	if (varILev == NULL) {
		Announce("Warning: Variable \"ilev\" not found");
	} else {
		nILev = varILev->get_dim(0)->size();
		if (nILev != nLev + 1) {
			_EXCEPTIONT("Variable \"ilev\" must have size lev+1");
		}
		dILev.Allocate(nILev);
		varILev->set_cur((long)0);
		varILev->get(&(dILev[0]), nILev);
	}

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

	// Open output file
	AnnounceStartBlock("Constructing output file");

	NcFile ncdf_out(strOutputFile.c_str(), NcFile::Replace);
	if (!ncdf_out.is_valid()) {
		_EXCEPTION1("Unable to open file \"%s\" for writing",
			strOutputFile.c_str());
	}

	CopyNcFileAttributes(&ncdf_in, &ncdf_out);

	// Output time array
	Announce("Time");
	NcDim * dimOutTime = ncdf_out.add_dim("time");
	NcVar * varOutTime = ncdf_out.add_var("time", ncDouble, dimOutTime);
	varOutTime->set_cur((long)0);
	varOutTime->put(&(dTime[0]), nTime);

	CopyNcVarAttributes(varTime, varOutTime);

	// Output pressure array
	NcDim * dimOutP = NULL;
	NcVar * varOutP = NULL;
	if (nPressureLevels > 0) {
		Announce("Pressure");
		dimOutP = ncdf_out.add_dim("p", nPressureLevels);
		varOutP = ncdf_out.add_var("p", ncDouble, dimOutP);
		varOutP->set_cur((long)0);
		varOutP->put(&(vecPressureLevels[0]), nPressureLevels);
	}

	// Output height array
	NcDim * dimOutZ = NULL;
	NcVar * varOutZ = NULL;
	if (nHeightLevels > 0) {
		Announce("Height");
		dimOutZ = ncdf_out.add_dim("z", nHeightLevels);
		varOutZ = ncdf_out.add_var("z", ncDouble, dimOutZ);
		varOutZ->set_cur((long)0);
		varOutZ->put(&(vecHeightLevels[0]), nHeightLevels);
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

	// Output topography
	Announce("Topography");
	NcVar * varOutZs = ncdf_out.add_var(
		"Zs", ncDouble, dimOutLat, dimOutLon);

	varOutZs->set_cur((long)0, (long)0);
	varOutZs->put(&(dZs[0][0]), nLat, nLon);

	AnnounceEndBlock("Done");

	// Done
	AnnounceEndBlock("Done");

	// Load all variables
	Announce("Loading variables");

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

	NcAtt * attEarthRadius = ncdf_in.get_att("earth_radius");
	double dEarthRadius = attEarthRadius->as_double(0);

	NcAtt * attRd = ncdf_in.get_att("Rd");
	double dRd = attRd->as_double(0);

	NcAtt * attCp = ncdf_in.get_att("Cp");
	double dCp = attCp->as_double(0);

	double dGamma = dCp / (dCp - dRd);

	NcAtt * attP0 = ncdf_in.get_att("P0");
	double dP0 = attP0->as_double(0);

	double dPressureScaling = dP0 * std::pow(dRd / dP0, dGamma);

	NcAtt * attZtop = ncdf_in.get_att("Ztop");
	double dZtop = attZtop->as_double(0);

	// Input data
	DataArray3D<double> dataIn(nLev, nLat, nLon);
	DataArray3D<double> dataInt(nILev, nLat, nLon);

	// Output data
	DataArray2D<double> dataOut(nLat, nLon);

	// Pressure in column
	DataArray1D<double> dataColumnP(nLev);

	// Height in column
	DataArray1D<double> dataColumnZ(nLev);
	DataArray1D<double> dataColumnIZ(nILev);

	// Column weights
	DataArray1D<double> dW(nLev);
	DataArray1D<double> dIW(nILev);

	// Loop through all times, pressure levels and variables
	AnnounceStartBlock("Interpolating");

	// Add energy variable
	NcVar * varEnergy;
	if (fExtractTotalEnergy) {
		varEnergy = ncdf_out.add_var("TE", ncDouble, dimOutTime);
	}

	// Create output pressure variables
	std::vector<NcVar *> vecOutNcVarP;
	if (nPressureLevels > 0) {
		for (int v = 0; v < vecVariableStrings.size(); v++) {
			vecOutNcVarP.push_back(
				ncdf_out.add_var(
					vecVariableStrings[v].c_str(), ncDouble,
						dimOutTime, dimOutP, dimOutLat, dimOutLon));

			// Copy attributes
			CopyNcVarAttributes(vecNcVar[v], vecOutNcVarP[v]);
		}
	}

	// Create output height variables
	std::vector<NcVar *> vecOutNcVarZ;
	if (nHeightLevels > 0) {
		for (int v = 0; v < vecVariableStrings.size(); v++) {
			std::string strVarName = vecVariableStrings[v];
			if (nPressureLevels > 0) {
				strVarName += "z";
			}
			vecOutNcVarZ.push_back(
				ncdf_out.add_var(
					strVarName.c_str(), ncDouble,
						dimOutTime, dimOutZ, dimOutLat, dimOutLon));

			// Copy attributes
			CopyNcVarAttributes(vecNcVar[v], vecOutNcVarZ[v]);
		}
	}

	// Create output surface variable
	std::vector<NcVar *> vecOutNcVarS;
	if (fExtractSurface) {
		for (int v = 0; v < vecVariableStrings.size(); v++) {
			std::string strVarName = vecVariableStrings[v];
			strVarName += "S";

			vecOutNcVarS.push_back(
				ncdf_out.add_var(
					strVarName.c_str(), ncDouble,
						dimOutTime, dimOutLat, dimOutLon));

			// Copy attributes
			CopyNcVarAttributes(vecNcVar[v], vecOutNcVarS[v]);
		}
	}

	// Loop over all times
	for (int t = 0; t < nTime; t++) {

		char szAnnounce[256];
		sprintf(szAnnounce, "Time %i", t); 
		AnnounceStartBlock(szAnnounce);

		// Rho
		DataArray3D<double> dataRho(nLev, nLat, nLon);

		NcVar * varRho = ncdf_in.get_var("Rho");
		if (varRho == NULL) {
			_EXCEPTIONT("Unable to load variable \"Rho\" from file");
		}
		varRho->set_cur(t, 0, 0, 0);
		varRho->get(&(dataRho[0][0][0]), 1, nLev, nLat, nLon);

		// Pressure
		DataArray3D<double> dataP(nLev, nLat, nLon);

		if (nPressureLevels != 0) {
			NcVar * varP = ncdf_in.get_var("P");
			if (varP == NULL) {
				_EXCEPTIONT("Unable to load variable \"P\" from file");
			}
			varP->set_cur(t, 0, 0, 0);
			varP->get(&(dataP[0][0][0]), 1, nLev, nLat, nLon);
		}
/*
		// Populate pressure array
		if (nPressureLevels > 0) {

			// Calculate pointwise pressure
			for (int k = 0; k < nLev; k++) {
			for (int i = 0; i < nLat; i++) {
			for (int j = 0; j < nLon; j++) {
				dataP[k][i][j] = dPressureScaling
					* exp(log(dataRho[k][i][j] * dataP[k][i][j]) * dGamma);
			}
			}
			}
		}
*/
		// Height everywhere
		DataArray3D<double> dataZ(nLev, nLat, nLon);
		DataArray3D<double> dataIZ;
		if (nILev != 0) {
			dataIZ.Allocate(nILev, nLat, nLon);
		}

		// Populate height array
		if ((nHeightLevels > 0) || (fGeopotentialHeight)) {
			for (int k = 0; k < nLev; k++) {
			for (int i = 0; i < nLat; i++) {
			for (int j = 0; j < nLon; j++) {
				dataZ[k][i][j] = dZs[i][j] + dLev[k] * (dZtop - dZs[i][j]);
			}
			}
			}

			for (int k = 0; k < nILev; k++) {
			for (int i = 0; i < nLat; i++) {
			for (int j = 0; j < nLon; j++) {
				dataIZ[k][i][j] = dZs[i][j] + dILev[k] * (dZtop - dZs[i][j]);
			}
			}
			}
		}

		// Loop through all pressure levels and variables
		for (int v = 0; v < vecNcVar.size(); v++) {

			bool fOnInterfaces = false;

			// Load in the data array
			vecNcVar[v]->set_cur(t, 0, 0, 0);

			if (vecNcVar[v]->get_dim(1)->size() == nLev) {
				vecNcVar[v]->get(&(dataIn[0][0][0]), 1, nLev, nLat, nLon);

				Announce("%s (n)", vecVariableStrings[v].c_str());

			} else if (vecNcVar[v]->get_dim(1)->size() == nILev) {
				vecNcVar[v]->get(&(dataInt[0][0][0]), 1, nILev, nLat, nLon);
				fOnInterfaces = true;

				Announce("%s (i)", vecVariableStrings[v].c_str());
			} else {
				_EXCEPTION1("Variable \"%s\" has invalid level dimension",
					vecVariableStrings[v].c_str());
			}

			// At the physical surface
			if (fExtractSurface) {

				if (fOnInterfaces) {
					for (int i = 0; i < nLat; i++) {
					for (int j = 0; j < nLon; j++) {
						dataOut[i][j] = dataInt[0][i][j];
					}
					}

				} else {

					int kBegin = 0;
					int kEnd = 3;

					PolynomialInterp::LagrangianPolynomialCoeffs(
						3, dLev, dW, 0.0);

					// Loop thorugh all latlon indices
					for (int i = 0; i < nLat; i++) {
					for (int j = 0; j < nLon; j++) {

						// Interpolate in the vertical
						dataOut[i][j] = 0.0;
						for (int k = kBegin; k < kEnd; k++) {
							dataOut[i][j] += dW[k] * dataIn[k][i][j];
						}
					}
					}
				}

				// Write variable
				vecOutNcVarS[v]->set_cur(t, 0, 0);
				vecOutNcVarS[v]->put(&(dataOut[0][0]), 1, nLat, nLon);

			}

			// Loop through all pressure levels
			for (int p = 0; p < nPressureLevels; p++) {

				// Loop thorugh all latlon indices
				for (int i = 0; i < nLat; i++) {
				for (int j = 0; j < nLon; j++) {

					// Store column pressure
					for (int k = 0; k < nLev; k++) {
						dataColumnP[k] = dataP[k][i][j];
					}

					// Find weights
					int kBegin = 0;
					int kEnd = 0;

					// On a pressure surface
					InterpolationWeightsLinear(
						vecPressureLevels[p],
						dataColumnP,
						kBegin,
						kEnd,
						dW);

					// Interpolate in the vertical
					dataOut[i][j] = 0.0;
					for (int k = kBegin; k < kEnd; k++) {
						dataOut[i][j] += dW[k] * dataIn[k][i][j];
					}

				}
				}

				// Write variable
				vecOutNcVarP[v]->set_cur(t, p, 0, 0);
				vecOutNcVarP[v]->put(&(dataOut[0][0]), 1, 1, nLat, nLon);
			}

			// Loop through all height levels
			for (int z = 0; z < nHeightLevels; z++) {

				// Loop thorugh all latlon indices
				for (int i = 0; i < nLat; i++) {
				for (int j = 0; j < nLon; j++) {

					// Find weights
					int kBegin = 0;
					int kEnd = 0;

					// Interpolate from levels to z surfaces
					if (!fOnInterfaces) {
						for (int k = 0; k < nLev; k++) {
							dataColumnZ[k] = dataZ[k][i][j];
						}

						InterpolationWeightsLinear(
							vecHeightLevels[z],
							dataColumnZ,
							kBegin,
							kEnd,
							dW);

						dataOut[i][j] = 0.0;
						for (int k = kBegin; k < kEnd; k++) {
							dataOut[i][j] += dW[k] * dataIn[k][i][j];
						}

					// Interpolate from interfaces to z surfaces
					} else {
						for (int k = 0; k < nILev; k++) {
							dataColumnIZ[k] = dataIZ[k][i][j];
						}

						InterpolationWeightsLinear(
							vecHeightLevels[z],
							dataColumnIZ,
							kBegin,
							kEnd,
							dIW);

						dataOut[i][j] = 0.0;
						for (int k = kBegin; k < kEnd; k++) {
							dataOut[i][j] += dIW[k] * dataInt[k][i][j];
						}
					}
				}
				}

				// Write variable
				vecOutNcVarZ[v]->set_cur(t, z, 0, 0);
				vecOutNcVarZ[v]->put(&(dataOut[0][0]), 1, 1, nLat, nLon);
			}
		}

		// Output geopotential height
		if (fGeopotentialHeight) {

			Announce("Geopotential height");

			// Output variables
			NcVar * varOutZ;
			NcVar * varOutZs;

			if (nPressureLevels > 0) {
				varOutZ = ncdf_out.add_var(
					"PHIZ", ncDouble, dimOutTime, dimOutP, dimOutLat, dimOutLon);
			}
			if (fExtractSurface) {
				varOutZs = ncdf_out.add_var(
					"PHIZS", ncDouble, dimOutTime, dimOutLat, dimOutLon);
			}

			// Interpolate onto pressure levels
			for (int p = 0; p < nPressureLevels; p++) {

				// Loop thorugh all latlon indices
				for (int i = 0; i < nLat; i++) {
				for (int j = 0; j < nLon; j++) {

					int kBegin = 0;
					int kEnd = 0;

					for (int k = 0; k < nLev; k++) {
						dataColumnP[k] = dataP[k][i][j];
					}

					InterpolationWeightsLinear(
						vecPressureLevels[p],
						dataColumnP,
						kBegin,
						kEnd,
						dW);

					// Interpolate in the vertical
					dataOut[i][j] = 0.0;
					for (int k = kBegin; k < kEnd; k++) {
						dataOut[i][j] += dW[k] * dataZ[k][i][j];
					}
				}
				}

				// Write variable
				varOutZ->set_cur(t, p, 0, 0);
				varOutZ->put(&(dataOut[0][0]), 1, 1, nLat, nLon);

			}

			// Interpolate onto the physical surface
			if (fExtractSurface) {

				int kBegin = 0;
				int kEnd = 3;

				PolynomialInterp::LagrangianPolynomialCoeffs(
					3, dLev, dW, 0.0);

				// Loop thorugh all latlon indices
				for (int i = 0; i < nLat; i++) {
				for (int j = 0; j < nLon; j++) {

					// Interpolate in the vertical
					dataOut[i][j] = 0.0;
					for (int k = kBegin; k < kEnd; k++) {
						dataOut[i][j] += dW[k] * dataZ[k][i][j];
					}
				}
				}

				// Write variable
				varOutZs->set_cur(t, 0, 0);
				varOutZs->put(&(dataOut[0][0]), 1, nLat, nLon);

			}
		}

		// Extract total energy
		if (fExtractTotalEnergy) {
			Announce("Total Energy");

			// Zonal velocity
			DataArray3D<double> dataU(nLev, nLat, nLon);

			NcVar * varU = ncdf_in.get_var("U");
			varU->set_cur(t, 0, 0, 0);
			varU->get(&(dataU[0][0][0]), 1, nLev, nLat, nLon);

			// Meridional velocity
			DataArray3D<double> dataV(nLev, nLat, nLon);

			NcVar * varV = ncdf_in.get_var("V");
			varV->set_cur(t, 0, 0, 0);
			varV->get(&(dataV[0][0][0]), 1, nLev, nLat, nLon);

			// Vertical velocity
			DataArray3D<double> dataW(nLev, nLat, nLon);

			NcVar * varW = ncdf_in.get_var("W");
			varW->set_cur(t, 0, 0, 0);
			varW->get(&(dataW[0][0][0]), 1, nLev, nLat, nLon);

			// Calculate total energy
			double dTotalEnergy = 0.0;

			double dElementRefArea =
				dEarthRadius * dEarthRadius
				* M_PI / static_cast<double>(nLat)
				* 2.0 * M_PI / static_cast<double>(nLon);

			for (int k = 0; k < nLev; k++) {
			for (int i = 0; i < nLat; i++) {
			for (int j = 0; j < nLon; j++) {
				double dKineticEnergy =
					0.5 * dataRho[k][i][j] *
						( dataU[k][i][j] * dataU[k][i][j]
						+ dataV[k][i][j] * dataV[k][i][j]
						+ dataW[k][i][j] * dataW[k][i][j]);

				double dInternalEnergy =
					dataP[k][i][j] / (dGamma - 1.0);

				dTotalEnergy +=
					(dKineticEnergy + dInternalEnergy)
						* std::cos(M_PI * dLat[i] / 180.0) * dElementRefArea
						* (dZtop - dZs[i][j]) / static_cast<double>(nLev);
			}
			}
			}

			// Put total energy into file
			varEnergy->set_cur(t);
			varEnergy->put(&dTotalEnergy, 1);
		}

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

