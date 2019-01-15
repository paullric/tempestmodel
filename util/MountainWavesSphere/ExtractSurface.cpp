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

#include <cstdlib>
#include <netcdfcpp.h>

#include "PolynomialInterp.h"

#include <mpi.h>

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char ** argv) {

	// Initialize MPI
	MPI_Init(&argc, &argv);

try {

	// Input filename
	std::string strInputFile;

	// Output filename
	std::string strOutputFile;

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


	AnnounceEndBlock("Done");

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}

	// Finalize MPI
	MPI_Finalize();

}

///////////////////////////////////////////////////////////////////////////////

