///////////////////////////////////////////////////////////////////////////////
///
///	\file    Parameters.h
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

#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

#include <vector>
#include <cmath>

#include <netcdfcpp.h>

///////////////////////////////////////////////////////////////////////////////

class Parameters {
public:

	///	<summary>
	///		Number of elements in the phi direction
	///	</summary>
	int nPhiElements;

	///	<summary>
	///		Earth scale factor
	///	</summary>
	double dXscale;
	
	///	<summary>
	///		Background temperature
	///	</summary>
	double dT0;
	
	///	<summary>
	///		Background flow velocity
	///	</summary>
	double dU0;
	
	///	<summary>
	///		Gravity
	///	</summary>
	double dG;
	
	///	<summary>
	///		Base rotation rate
	///	</summary>
	double dOmega;
	
	///	<summary>
	///		Gamma
	///	</summary>
	double dGamma;
	
	///	<summary>
	///		Nodal points
	///	</summary>
	std::vector<double> vecNode;
	
	///	<summary>
	///		Edge points
	///	</summary>
	std::vector<double> vecEdge;

public:
	///	<summary>
	///		Generate latitude arrays.
	///	</summary>
	void GenerateLatituteArray(
		int a_nPhiElements
	) {
		nPhiElements = a_nPhiElements;

		// Generate the node and edge coordinate array
		vecNode.resize(nPhiElements);
		vecEdge.resize(nPhiElements+1);

		double dDeltaPhi = M_PI / static_cast<double>(nPhiElements);

		//PhiX = -pi/2 * cos((2*[1:Nphi]-1)*pi/(2*Nphi));
		double dPhiElements = static_cast<double>(nPhiElements);

		for (int i = 0; i < nPhiElements; i++) {
			double dI = static_cast<double>(i);

			vecNode[i] =
				- 0.5 * M_PI + (static_cast<double>(i) + 0.5) * dDeltaPhi;
				//- 0.5 * M_PI
				//	* cos((2.0 * dI + 1.0) * M_PI / (2.0 * dPhiElements));
		}
		for (int i = 0; i <= nPhiElements; i++) {
			vecEdge[i] =
				- 0.5 * M_PI + static_cast<double>(i) * dDeltaPhi;
		}
	}

public:
	///	<summary>
	///		Write parameters to a NetCDF file.
	///	</summary>
	void WriteToNcFile(
		NcFile & ncdf_out,
		NcDim * dimLat,
		NcDim * dimLatS = NULL
	) {
		// Scalar parameters
		ncdf_out.add_att("X", dXscale);
		ncdf_out.add_att("T0", dT0);
		ncdf_out.add_att("u0", dU0);
		ncdf_out.add_att("g", dG);
		ncdf_out.add_att("omega", dOmega);
		ncdf_out.add_att("gamma", dGamma);

		// Latitudes
		NcVar *varLat = ncdf_out.add_var("lat", ncDouble, dimLat);
	
		varLat->set_cur((long)0);
		varLat->put(&(vecNode[0]), nPhiElements);

		if (dimLatS != NULL) {
			NcVar *varLatS = ncdf_out.add_var("lats", ncDouble, dimLatS);

			varLatS->set_cur((long)0);
			varLatS->put(&(vecEdge[1]), nPhiElements - 1);
		}
	}

	///	<summary>
	///		Read parameters from a NetCDF file.
	///	</summary>
	void ReadFromNcFile(
		NcFile & ncdf_in
	) {
		// Scalar parameters
		dXscale = ncdf_in.get_att("X")->as_double(0);
		dT0 = ncdf_in.get_att("T0")->as_double(0);
		dU0 = ncdf_in.get_att("u0")->as_double(0);
		dG  = ncdf_in.get_att("g")->as_double(0);
		dOmega = ncdf_in.get_att("omega")->as_double(0);
		dGamma = ncdf_in.get_att("gamma")->as_double(0);

		// Latitudes
		NcVar *varLat = ncdf_in.get_var("lat");
		nPhiElements = varLat->get_dim(0)->size();

		vecNode.resize(nPhiElements);
		varLat->set_cur((long)0);
		varLat->get(&(vecNode[0]), nPhiElements);
/*
		vecEdge.resize(nPhiElements+1);
		vecEdge[0] = - 0.5 * M_PI;
		vecEdge[1] = + 0.5 * M_PI;

		NcVar *varLatS = ncdf_in.get_var("lats");
		varLatS->set_cur((long)0);
		varLatS->get(&(vecEdge[1]), nPhiElements-1);
*/
	}

};

///////////////////////////////////////////////////////////////////////////////

#endif

