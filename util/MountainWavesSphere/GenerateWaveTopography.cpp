///////////////////////////////////////////////////////////////////////////////
///
///	\file    GenerateWaveTopography.cpp
///	\author  Paul Ullrich
///	\version April 11, 2014
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

#include "Exception.h"
#include "CommandLine.h"
#include "DataVector.h"
#include "DataMatrix.h"
#include "Announce.h"
#include "Parameters.h"

#include <cmath>

#include <netcdfcpp.h>

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char ** argv) {

try {

	// Number of latitudes
	int nLat;

	// Number of longitudes
	int nLon;

	// Zonal wavenumber
	int nK;

	// Meridional power
	int nLpow;

	// Output file
	std::string strOutputFile;

	// Parse the command line
	BeginCommandLine()
		CommandLineInt(nLat, "lat", 40);
		CommandLineInt(nLon, "lon", 80);
		CommandLineInt(nK, "k", 6);
		CommandLineInt(nLpow, "lpow", 2);
		CommandLineString(strOutputFile, "out", "topo.nc");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	// Generate longitude and latitude arrays
	AnnounceBanner();
	AnnounceStartBlock("Generating longitude and latitude arrays");

	DataVector<double> dLon;
	dLon.Initialize(nLon);

	Parameters param;
	param.GenerateLatituteArray(nLat);

	double dDeltaLon = 2.0 * M_PI / static_cast<double>(nLon);
	for (int i = 0; i < nLon; i++) {
		dLon[i] = (static_cast<double>(i) + 0.5) * dDeltaLon;
	}

	AnnounceEndBlock("Done");

	// Open NetCDF output file
	AnnounceStartBlock("Writing to file");

	NcFile ncdf_out(strOutputFile.c_str(), NcFile::Replace);

	// Output coordinates
	NcDim * dimLat = ncdf_out.add_dim("lat", nLat);
	NcDim * dimLon = ncdf_out.add_dim("lon", nLon);

	NcVar * varLat = ncdf_out.add_var("lat", ncDouble, dimLat);
	varLat->set_cur((long)0);
	varLat->put(&(param.vecNode[0]), nLat);

	NcVar * varLon = ncdf_out.add_var("lon", ncDouble, dimLon);
	varLon->set_cur((long)0);
	varLon->put(&(dLon[0]), nLon);

	// Generate topography
	DataMatrix<double> dTopo;
	dTopo.Initialize(nLat, nLon);

	double dK = static_cast<double>(nK);
	double dLpow = static_cast<double>(nLpow);

	for (int j = 0; j < nLat; j++) {
	for (int i = 0; i < nLon; i++) {
		dTopo[j][i] = sin(dK * dLon[i]) * pow(cos(param.vecNode[j]), dLpow);
	}
	}

	// Write topography
	NcVar * varZs = ncdf_out.add_var("Zs", ncDouble, dimLat, dimLon);
	varZs->set_cur(0, 0);
	varZs->put(&(dTopo[0][0]), nLat, nLon);

	AnnounceEndBlock("Done");
	Announce("Completed successfully!");
	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}
}

///////////////////////////////////////////////////////////////////////////////

