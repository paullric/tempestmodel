///////////////////////////////////////////////////////////////////////////////
///
///	\file    ZonalTemporalAverage.cpp
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

#include <mpi.h>

#include <cstdlib>

#include <netcdfcpp.h>

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

void SetupArraysFromFile(
	const std::string & strSourceFile,
	const std::vector<std::string> & vecVariableStrings, 
	NcFile & ncdf_out,
	std::vector<NcVar *> & vecOutputVars,
	NcDim ** pDimOutLev,
	NcDim ** pDimOutLat
) {
	AnnounceStartBlock("Setup arrays from file");

	// Open source file
	NcFile ncdf_in(strSourceFile.c_str(), NcFile::ReadOnly);
	if (!ncdf_in.is_valid()) {
		_EXCEPTION1("Unable to open file \"%s\" for reading",
			strSourceFile.c_str());
	}

	// Copy file attributes
	CopyNcFileAttributes(&ncdf_in, &ncdf_out);

	// Load first variable
	NcVar * varIn = ncdf_in.get_var(vecVariableStrings[0].c_str());

	// Dimension names
	NcDim * dimLev = varIn->get_dim(1);
	NcDim * dimLat = varIn->get_dim(2);
	NcDim * dimLon = varIn->get_dim(3);

	NcToken szDimLevName = dimLev->name();
	NcToken szDimLatName = dimLat->name();
	NcToken szDimLonName = dimLon->name();

	int nLev = dimLev->size();
	int nLat = dimLat->size();
	int nLon = dimLon->size();

	NcVar * varLev = ncdf_in.get_var(szDimLevName);
	NcVar * varLat = ncdf_in.get_var(szDimLatName);
	NcVar * varLon = ncdf_in.get_var(szDimLonName);

	// Load level array
	DataVector<double> dLev;
	dLev.Initialize(nLev);
	varLev->set_cur((long)0);
	varLev->get(&(dLev[0]), nLev);

	// Load latitude array
	DataVector<double> dLat;
	dLat.Initialize(nLat);
	varLat->set_cur((long)0);
	varLat->get(&(dLat[0]), nLat);

	// Load longitude array
	DataVector<double> dLon;
	dLon.Initialize(nLon);
	varLon->set_cur((long)0);
	varLon->get(&(dLon[0]), nLon);

	// Add level array
	NcDim * dimOutLev = ncdf_out.add_dim(szDimLevName, nLev);
	NcVar * varOutLev = ncdf_out.add_var(szDimLevName, ncDouble, dimOutLev);
	varOutLev->set_cur((long)0);
	varOutLev->put(&(dLev[0]), nLev);

	// Add latitude array
	NcDim * dimOutLat = ncdf_out.add_dim(szDimLatName, nLat);
	NcVar * varOutLat = ncdf_out.add_var(szDimLatName, ncDouble, dimOutLat);
	varOutLat->set_cur((long)0);
	varOutLat->put(&(dLat[0]), nLat);

	// Store output dimensions
	(*pDimOutLev) = dimOutLev;
	(*pDimOutLat) = dimOutLat;

/*
	// Add longitude array
	NcDim * dimOutLon = ncdf_out.add_dim(szDimLonName, nLon);
	NcVar * varOutLon = ncdf_out.add_var(szDimLonName, ncDouble, dimOutLon);
	varOutLon->set_cur((long)0);
	varOutLon->put(&(dLon[0]), nLon);
*/
	// Add variable array
	for (int v = 0; v < vecVariableStrings.size(); v++) {
		NcVar * varIn =
			ncdf_in.get_var(vecVariableStrings[v].c_str());

		NcVar * varOut =
			ncdf_out.add_var(
				vecVariableStrings[v].c_str(),
				ncDouble,
				dimOutLev,
				dimOutLat);

		CopyNcVarAttributes(varIn, varOut);

		vecOutputVars.push_back(varOut);
	}

	AnnounceEndBlock("Done");
}

///////////////////////////////////////////////////////////////////////////////

void AddToAverage(
	const DataMatrix3D<double> & dataIn,
	DataMatrix3D<double> & dataRunningAvg
) {
	for (int k = 0; k < dataIn.GetRows(); k++) {
	for (int i = 0; i < dataIn.GetColumns(); i++) {
	for (int j = 0; j < dataIn.GetSubColumns(); j++) {
		dataRunningAvg[k][i][j] += dataIn[k][i][j];
	}
	}
	}
}

///////////////////////////////////////////////////////////////////////////////

void TimeAverage(
	const std::vector<std::string> & strFileList,
	const std::vector<std::string> & vecVariableStrings,
	std::vector< DataMatrix3D<double> > & vecDataAverage
) {
	AnnounceStartBlock("Temporal Average");

	// Check input parameters
	if (vecVariableStrings.size() != vecDataAverage.size()) {
		_EXCEPTIONT("Array size mismatch");
	}
	if (vecDataAverage.size() == 0) {
		_EXCEPTIONT("Zero size DataAverage array received");
	}

	// Total time 
	int nTotalTime = 0;

	// Reference number of elements
	int nRefLev;
	int nRefLat;
	int nRefLon;

	// Temporary data storage for input data
	DataMatrix3D<double> dataInput;

	// Loop through all files
	for (int f = 0; f < strFileList.size(); f++) {

		// Announce file
		Announce("%s", strFileList[f].c_str());

		// Open the file
		NcFile ncdf_in(strFileList[f].c_str(), NcFile::ReadOnly);
		if (!ncdf_in.is_valid()) {
			_EXCEPTION1("Unable to open file \"%s\" for reading",
				strFileList[f].c_str());
		}

		// Get the time dimension
		NcDim * dimTime = ncdf_in.get_dim("time");
		int nTime = dimTime->size();

		// Loop through all times
		for (int t = 0; t < nTime; t++) {

			// Loop through all variables
			for (int v = 0; v < vecVariableStrings.size(); v++) {

				NcVar * varIn =
					ncdf_in.get_var(vecVariableStrings[v].c_str());

				// Load level array
				NcDim * dimLev = varIn->get_dim(1);
				int nLev = dimLev->size();

				// Load latitude array
				NcDim * dimLat = varIn->get_dim(2);
				int nLat = dimLat->size();

				// Load longitude array
				NcDim * dimLon = varIn->get_dim(3);
				int nLon = dimLon->size();

				// Allocate space for data
				if ((f == 0) && (t == 0) && (v == 0)) {

					// Store reference sizes
					nRefLev = nLev;
					nRefLat = nLat;
					nRefLon = nLon;

					dataInput.Initialize(nRefLev, nRefLat, nRefLon);

					for (int w = 0; w < vecDataAverage.size(); w++) {
						vecDataAverage[w].Initialize(nLev, nLat, nLon);
					}
				}

				// Check array sizes
				if ((nRefLev != nLev) ||
					(nRefLat != nLat) ||
					(nRefLon != nLon)
				) {
					_EXCEPTIONT("Dimension size mismatch");
				}

				// Load data
				varIn->set_cur(t, 0, 0, 0);
				varIn->get(&(dataInput[0][0][0]), 1, nLev, nLat, nLon);

				// Add to average
				for (int k = 0; k < nLev; k++) {
				for (int i = 0; i < nLat; i++) {
				for (int j = 0; j < nLon; j++) {
					vecDataAverage[v][k][i][j] += dataInput[k][i][j];
				}
				}
				}
			}
		}

		// Advance the total time
		nTotalTime += nTime;
	}

	// Average by total time increments
	for (int v = 0; v < vecDataAverage.size(); v++) {
		for (int k = 0; k < nRefLev; k++) {
		for (int i = 0; i < nRefLat; i++) {
		for (int j = 0; j < nRefLon; j++) {
			vecDataAverage[v][k][i][j] /= 
				static_cast<double>(nTotalTime);
		}
		}
		}
	}

	AnnounceEndBlock("Done");
}

///////////////////////////////////////////////////////////////////////////////

void EddyStatistics(
	const std::vector<std::string> & strFileList,
	const std::vector<std::string> & vecVariableStrings,
	const std::vector< DataMatrix<double> > & dataDifference,
	std::vector< DataMatrix3D<double> > & vecDataAverage
) {
	// Find eddy variables
	int ixVarU = (-1);
	int ixVarV = (-1);
	int ixVarT = (-1);

	for (int v = 0; v < vecVariableStrings.size(); v++) {
		if (vecVariableStrings[v] == "U") {
			ixVarU = v;
		} else if (vecVariableStrings[v] == "V") {
			ixVarV = v;
		} else if (vecVariableStrings[v] == "T") {
			ixVarT = v;
		}
	}

	if (ixVarU == (-1)) {
		_EXCEPTIONT("Variable U must be included in variable list"
			" for eddy statistics");
	}
	if (ixVarV == (-1)) {
		_EXCEPTIONT("Variable V must be included in variable list"
			" for eddy statistics");
	}
	if (ixVarT == (-1)) {
		_EXCEPTIONT("Variable T must be included in variable list"
			" for eddy statistics");
	}

	// Total time 
	int nTotalTime = 0;

	// Reference number of elements
	int nRefLev;
	int nRefLat;
	int nRefLon;

	// Temporary data storage for input data
	DataMatrix3D<double> dataInputU;
	DataMatrix3D<double> dataInputV;
	DataMatrix3D<double> dataInputT;

	// Loop through all files
	for (int f = 0; f < strFileList.size(); f++) {

		// Announce file
		Announce("%s", strFileList[f].c_str());

		// Open the file
		NcFile ncdf_in(strFileList[f].c_str(), NcFile::ReadOnly);
		if (!ncdf_in.is_valid()) {
			_EXCEPTION1("Unable to open file \"%s\" for reading",
				strFileList[f].c_str());
		}

		// Get the time dimension
		NcDim * dimTime = ncdf_in.get_dim("time");
		int nTime = dimTime->size();

		// Loop through all times
		for (int t = 0; t < nTime; t++) {

			// Primitive variables
			NcVar * varInU = ncdf_in.get_var("U");
			NcVar * varInV = ncdf_in.get_var("V");
			NcVar * varInT = ncdf_in.get_var("T");

			// Load level array
			NcDim * dimLev = varInU->get_dim(1);
			int nLev = dimLev->size();

			// Load latitude array
			NcDim * dimLat = varInU->get_dim(2);
			int nLat = dimLat->size();

			// Load longitude array
			NcDim * dimLon = varInU->get_dim(3);
			int nLon = dimLon->size();

			// Allocate space for data
			if ((f == 0) && (t == 0)) {

				// Store reference sizes
				nRefLev = nLev;
				nRefLat = nLat;
				nRefLon = nLon;

				dataInputU.Initialize(nRefLev, nRefLat, nRefLon);
				dataInputV.Initialize(nRefLev, nRefLat, nRefLon);
				dataInputT.Initialize(nRefLev, nRefLat, nRefLon);

				// Allocate input data
				vecDataAverage.resize(5);
				vecDataAverage[0].Initialize(nLev, nLat, nLon);
				vecDataAverage[1].Initialize(nLev, nLat, nLon);
				vecDataAverage[2].Initialize(nLev, nLat, nLon);
				vecDataAverage[3].Initialize(nLev, nLat, nLon);
				vecDataAverage[4].Initialize(nLev, nLat, nLon);
			}

			// Check array sizes
			if ((nRefLev != nLev) ||
				(nRefLat != nLat) ||
				(nRefLon != nLon)
			) {
				_EXCEPTIONT("Dimension size mismatch");
			}

			// Load eddy statistics variables
			varInU->set_cur(t, 0, 0, 0);
			varInU->get(&(dataInputU[0][0][0]), 1, nLev, nLat, nLon);

			varInV->set_cur(t, 0, 0, 0);
			varInV->get(&(dataInputV[0][0][0]), 1, nLev, nLat, nLon);

			varInT->set_cur(t, 0, 0, 0);
			varInT->get(&(dataInputT[0][0][0]), 1, nLev, nLat, nLon);

			// Add to average
			for (int k = 0; k < nLev; k++) {
			for (int i = 0; i < nLat; i++) {
			for (int j = 0; j < nLon; j++) {
				double dUprime =
					dataInputU[k][i][j] - dataDifference[ixVarU][k][i];
				double dVprime =
					dataInputV[k][i][j] - dataDifference[ixVarV][k][i];
				double dTprime =
					dataInputT[k][i][j] - dataDifference[ixVarT][k][i];

				vecDataAverage[0][k][i][j] += dUprime * dUprime;
				vecDataAverage[1][k][i][j] += dUprime * dVprime;
				vecDataAverage[2][k][i][j] += dVprime * dVprime;
				vecDataAverage[3][k][i][j] += dVprime * dTprime;
				vecDataAverage[4][k][i][j] += dTprime * dTprime;
			}
			}
			}
		}

		// Advance the total time
		nTotalTime += nTime;
	}

	// Average by total time increments
	for (int v = 0; v < vecDataAverage.size(); v++) {
		for (int k = 0; k < nRefLev; k++) {
		for (int i = 0; i < nRefLat; i++) {
		for (int j = 0; j < nRefLon; j++) {
			vecDataAverage[v][k][i][j] /= 
				static_cast<double>(nTotalTime);
		}
		}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void ZonalAverage(
	const DataMatrix3D<double> & dataIn,
	DataMatrix<double> & dataOut
) {
	dataOut.Initialize(
		dataIn.GetRows(),
		dataIn.GetColumns());

	dataOut.Zero();

	for (int k = 0; k < dataIn.GetRows(); k++) {
	for (int i = 0; i < dataIn.GetColumns(); i++) {
	for (int j = 0; j < dataIn.GetSubColumns(); j++) {
		dataOut[k][i] += dataIn[k][i][j];
	}
	}
	}

	for (int k = 0; k < dataIn.GetRows(); k++) {
	for (int i = 0; i < dataIn.GetColumns(); i++) {
		dataOut[k][i] /= static_cast<double>(dataIn.GetSubColumns());
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

	// Output filename
	std::string strVariables;

	// Eddy statistics
	bool fCalculateEddyStatistics;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "inlist", "");
		CommandLineString(strOutputFile, "out", "");
		CommandLineString(strVariables, "var", "");

		CommandLineBool(fCalculateEddyStatistics, "eddystats");

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

	// Open input file list
	AnnounceStartBlock("Loading input file list");
	FILE * fp = fopen(strInputFile.c_str(), "r");

	std::vector<std::string> strFileList;
	for (;;) {
		char szFilename[1024];
		fgets(szFilename, 1024, fp);

		if (feof(fp)) {
			break;
		}

		std::string strFilename(szFilename);
		strFilename[strFilename.length()-1] = '\0';

		strFileList.push_back(strFilename);
	}
	fclose(fp);

	if (strFileList.size() == 0) {
		_EXCEPTIONT("No input files specified");
	}
	AnnounceEndBlock("Done");

	// Output file
	NcFile ncdf_out(strOutputFile.c_str(), NcFile::Replace);

	// Setup output arrays from file
	std::vector<NcVar *> vecOutVariables;

	NcDim * dimOutLev;
	NcDim * dimOutLat;

	SetupArraysFromFile(
		strFileList[0],
		vecVariableStrings,
		ncdf_out,
		vecOutVariables,
		&dimOutLev,
		&dimOutLat);

	// Data arrays for storing time average
	std::vector< DataMatrix3D<double> > vecTimeAverage;
	vecTimeAverage.resize(vecVariableStrings.size());

	// Temporal average
	TimeAverage(
		strFileList,
		vecVariableStrings,
		vecTimeAverage);

	// Zonal average
	std::vector< DataMatrix<double> > vecZonalTimeAverage;

	vecZonalTimeAverage.resize(vecTimeAverage.size());

	AnnounceStartBlock("Zonal Average");
	for (int v = 0; v < vecTimeAverage.size(); v++) {
		ZonalAverage(
			vecTimeAverage[v],
			vecZonalTimeAverage[v]);
	}
	AnnounceEndBlock("Done");

	// Output
	AnnounceStartBlock("Output zonal averages");
	for (int v = 0; v < vecTimeAverage.size(); v++) {
		vecOutVariables[v]->put(
			&(vecZonalTimeAverage[v][0][0]),
			vecZonalTimeAverage[v].GetRows(),
			vecZonalTimeAverage[v].GetColumns());
	}
	AnnounceEndBlock("Done");

	// Calculate eddy statistics
	if (fCalculateEddyStatistics) {
		AnnounceStartBlock("Eddy statistics");

		std::vector< DataMatrix3D<double> > vecEddyTimeAverages;
		vecEddyTimeAverages.resize(5);

		EddyStatistics(
			strFileList,
			vecVariableStrings,
			vecZonalTimeAverage,
			vecEddyTimeAverages);

		std::vector< DataMatrix<double> > vecZonalTimeEddyAverage;
		vecZonalTimeEddyAverage.resize(5);
		for (int v = 0; v < vecEddyTimeAverages.size(); v++) {
			ZonalAverage(
				vecEddyTimeAverages[v],
				vecZonalTimeEddyAverage[v]);
		}

		int nLev = vecZonalTimeEddyAverage[0].GetRows();
		int nLat = vecZonalTimeEddyAverage[0].GetColumns();

		NcVar * varUU =
			ncdf_out.add_var("UU", ncDouble, dimOutLev, dimOutLat); 
		NcVar * varUV =
			ncdf_out.add_var("UV", ncDouble, dimOutLev, dimOutLat); 
		NcVar * varVV =
			ncdf_out.add_var("VV", ncDouble, dimOutLev, dimOutLat); 
		NcVar * varVT =
			ncdf_out.add_var("VT", ncDouble, dimOutLev, dimOutLat); 
		NcVar * varTT =
			ncdf_out.add_var("TT", ncDouble, dimOutLev, dimOutLat);

		varUU->put(&(vecZonalTimeEddyAverage[0][0][0]), nLev, nLat);
		varUV->put(&(vecZonalTimeEddyAverage[1][0][0]), nLev, nLat);
		varVV->put(&(vecZonalTimeEddyAverage[2][0][0]), nLev, nLat);
		varVT->put(&(vecZonalTimeEddyAverage[3][0][0]), nLev, nLat);
		varTT->put(&(vecZonalTimeEddyAverage[4][0][0]), nLev, nLat);

		AnnounceEndBlock("Done");
	}
	
} catch(Exception & e) {
	Announce(e.ToString().c_str());
}

	// Finalize MPI
	MPI_Finalize();

}

///////////////////////////////////////////////////////////////////////////////

