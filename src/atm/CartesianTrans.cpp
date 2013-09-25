////////////////////////////////////////////////////////////////////////////////
///
///	\file    CartesianTrans.cpp
///	\author  Paul Ullrich, Jorge Guerra
///	\version September 24, 2013
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

#include "CartesianTrans.h"
#include "Exception.h"

#include "MathHelper.h"
#include <cfloat>

////////////////////////////////////////////////////////////////////////////////

void CartesianTrans::XYPFromXYP(
	double dX_in,
	double dY_in,
	int nP,
	double &dX_out,
	double &dY_out
) {
	// Local variables
	double xx, yy, zz;
	double sx, sy, sz;

	// Convert to relative Cartesian coordinates
	sz = 1.0 / sqrt(1.0 + dX_in * dX_in + dY_in * dY_in);
	sx = sz * dX_in;
	sy = sz * dY_in;

	// Convert to full Cartesian coordinates
    yy = sx;
    zz = sy;
    xx = sz;

	// (DBG) Check the free z coordinate
	if (sz <= 0.0) {
		_EXCEPTIONT("Invalid relative Z coordinate");
	}

	// Use panel information to calculate (alpha, beta) coords
	dX_out = sx / sz;
	dY_out  = sy / sz;

	// Invalid panel
    if (nP != 0)
    {
      _EXCEPTION1(
      "Invalid Panel ID: nP. Given: %d, Expected: 0 for Cartesian Grid.\n", nP);
    }
}

////////////////////////////////////////////////////////////////////////////////

void CartesianTrans::RLLFromXYP(
	double dX,
	double dY,
	int nP,
	double &lon,
	double &lat
) {

	// Cartesian panel 0
    lon = atan(dX);
    lat = atan(dY / sqrt(1.0 + dX * dX));

    // Invalid panel
    if (nP != 0)
    {
      _EXCEPTION1(
      "Invalid Panel ID: nP. Given: %d, Expected: 0 for Cartesian Grid.\n", nP);
    }
}

////////////////////////////////////////////////////////////////////////////////

void CartesianTrans::RLLFromABP(
	double dA,
	double dB,
	int nP,
	double & lon,
	double & lat
) {
	return RLLFromXYP(tan(dA), tan(dB), nP, lon, lat);
}

////////////////////////////////////////////////////////////////////////////////

void CartesianTrans::XYPFromRLL(
	double lon,
	double lat,
	double &dX,
	double &dY,
	int &nP
) {
	// Default panel to unattainable value
	nP = 6;

	// Translate from RLL coordinates to XYZ space
	double xx, yy, zz, pm;

	xx = cos(lon) * cos(lat);
	yy = sin(lon) * cos(lat);
	zz = sin(lat);

	pm = Max(fabs(xx), Max(fabs(yy), fabs(zz)));

	// Check maxmality of the x coordinate
	if (pm == fabs(xx)) {
		if (xx > 0) {
			nP = 0;
		} else {
			nP = 2;
		}
	}

	// Check maximality of the y coordinate
	if (pm == fabs(yy)) {
		if (yy > 0) {
			nP = 1;
		} else {
			nP = 3;
		}
	}

	// Check maximality of the z coordinate
	if (pm == fabs(zz)) {
		if (zz > 0) {
			nP = 4;
		} else {
			nP = 5;
		}
	}

	// Panel assignments
	double sx, sy, sz;
	if (nP == 0) {
		sx = yy;
		sy = zz;
		sz = xx;

	} else if (nP == 1) {
		sx = -xx;
		sy = zz;
		sz = yy;

	} else if (nP == 2) {
		sx = -yy;
		sy = zz;
		sz = -xx;

	} else if (nP == 3) {
		sx = xx;
		sy = zz;
		sz = -yy;

	} else if (nP == 4) {
		sx = yy;
		sy = -xx;
		sz = zz;

	} else if (nP == 5) {
		sx = yy;
		sy = xx;
		sz = -zz;

	} else {
		_EXCEPTIONT("Logic error.");
	}

	// Convert to gnomonic coordinates
	dX = sx / sz;
	dY = sy / sz;
}

////////////////////////////////////////////////////////////////////////////////

void CartesianTrans::ABPFromRLL(
	double lon,
	double lat,
	double & dA,
	double & dB,
	int & nP
) {
	XYPFromRLL(lon, lat, dA, dB, nP);
	dA = atan(dA);
	dB = atan(dB);
}

////////////////////////////////////////////////////////////////////////////////
