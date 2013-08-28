///////////////////////////////////////////////////////////////////////////////
///
///	\file    PhysicalConstants.h
///	\author  Paul Ullrich
///	\version February 24, 2013
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

#include "PhysicalConstants.h"

#ifdef NETCDFENABLED

///////////////////////////////////////////////////////////////////////////////

void PhysicalConstants::NcOutputPhysicalConstants(
	NcFile * ncOut
) {
	// Add the const dimension to the file
	NcDim * dimConst = ncOut->add_dim("const", 1);

	// Output the Earth radius
	NcVar * varEarthRadius =
		ncOut->add_var("earth_radius", ncDouble, dimConst);
	varEarthRadius->put(&m_dEarthRadius, 1);

	// Output the constant of surface gravity
	NcVar * varGravity =
		ncOut->add_var("gravity", ncDouble, dimConst);
	varEarthRadius->put(&m_dG, 1);

	// Output the rotation rate of the Earth
	NcVar * varOmega =
		ncOut->add_var("omega", ncDouble, dimConst);
	varEarthRadius->put(&m_dOmega, 1);

	// Output the inclination of the grid
	NcVar * varAlpha =
		ncOut->add_var("alpha", ncDouble, dimConst);
	varEarthRadius->put(&m_dAlpha, 1);

	// Output the ideal gas constant of dry air
	NcVar * varR =
		ncOut->add_var("R_d", ncDouble, dimConst);
	varEarthRadius->put(&m_dR, 1);

	// Output the specific heat capacity of dry air at constant pressure
	NcVar * varCp =
		ncOut->add_var("c_p", ncDouble, dimConst);
	varEarthRadius->put(&m_dCp, 1);

	// Output the reference temperature
	NcVar * varT0 =
		ncOut->add_var("T0", ncDouble, dimConst);
	varEarthRadius->put(&m_dT0, 1);

	// Output the reference pressure
	NcVar * varP0 =
		ncOut->add_var("p0", ncDouble, dimConst);
	varEarthRadius->put(&m_dP0, 1);

	// Output the reference density of water
	NcVar * varRhoWater =
		ncOut->add_var("rho_water", ncDouble, dimConst);
	varRhoWater->put(&m_dRhoWater, 1);

	// Output the ideal gas constant for water vapor
	NcVar * varRvap =
		ncOut->add_var("R_v", ncDouble, dimConst);
	varRvap->put(&m_dRvap, 1);

	// Output the ratio of the molar mass of dry air to molar mass of water
	NcVar * varMvap =
		ncOut->add_var("M_v", ncDouble, dimConst);
	varRvap->put(&m_dMvap, 1);

	// Output the Latent heat of vaporization of water vapor
	NcVar * varLvap =
		ncOut->add_var("L_v", ncDouble, dimConst);
	varRvap->put(&m_dLvap, 1);
}

///////////////////////////////////////////////////////////////////////////////

void PhysicalConstants::NcInputPhysicalConstants(
	NcFile * ncIn
) {
}

///////////////////////////////////////////////////////////////////////////////

#else

///////////////////////////////////////////////////////////////////////////////

void PhysicalConstants::Stub() {
	_EXCEPTION();
}

///////////////////////////////////////////////////////////////////////////////

#endif
