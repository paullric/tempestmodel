///////////////////////////////////////////////////////////////////////////////
///
///	\file    CubedSphereTrans.cpp
///	\author  Paul Ullrich
///	\version March 6, 2013
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

#include "CubedSphereTrans.h"
#include "Exception.h"

#include "MathHelper.h"
#include <cfloat>

////////////////////////////////////////////////////////////////////////////////

void CubedSphereTrans::XYZFromXYP(
	double dXX,
	double dYY,
	int np,
	double &xx,
	double &yy,
	double &zz
) {
	// X, Y, Z coordinates in the local panel-centric coordinate system
	double sx, sy, sz;

	// Convert to Cartesian coordinates
	sz = 1.0 / sqrt(1.0 + dXX * dXX + dYY * dYY);
	sx = sz * dXX;
	sy = sz * dYY;

	// Rotate panel-centroid Cartesian coordinates to global
	// Cartesian coordinates.
	switch(np) {
		case 0:
			yy = sx;
			zz = sy;
			xx = sz;
			break;

		case 1:
			xx = -sx;
			zz = sy;
			yy = sz;
			break;

		case 2:
			yy = -sx;
			zz = sy;
			xx = -sz;
			break;

		case 3:
			xx = sx;
			zz = sy;
			yy = -sz;
			break;

		case 4:
			yy = sx;
			xx = -sy;
			zz = sz;
			break;

		case 5:
			yy = sx;
			xx = sy;
			zz = -sz;
			break;

		default:
			_EXCEPTIONT("Invalid panel id");
	}
}

////////////////////////////////////////////////////////////////////////////////

void CubedSphereTrans::XYPFromXYP(
	double dX_in,
	double dY_in,
	int isource,
	int idest,
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
	switch(isource) {
		case 0:
			yy = sx;
			zz = sy;
			xx = sz;
			break;

		case 1:
			xx = -sx;
			zz = sy;
			yy = sz;
			break;

		case 2:
			yy = -sx;
			zz = sy;
			xx = -sz;
			break;

		case 3:
			xx = sx;
			zz = sy;
			yy = -sz;
			break;

		case 4:
			yy = sx;
			xx = -sy;
			zz = sz;
			break;

		case 5:
			yy = sx;
			xx = sy;
			zz = -sz;
			break;

		default:
			_EXCEPTIONT("Invalid source panel id");
	}

	// Convert to relative Cartesian coordinates on destination panel
	switch (idest) {
		case 0:
			sx = yy;
			sy = zz;
			sz = xx;
			break;

		case 1:
			sx = -xx;
			sy = zz;
			sz = yy;
			break;

		case 2:
			sx = -yy;
			sy = zz;
			sz = -xx;
			break;

		case 3:
			sx = xx;
			sy = zz;
			sz = -yy;
			break;

		case 4:
			sx = yy;
			sy = -xx;
			sz = zz;
			break;

		case 5:
			sx = yy;
			sy = xx;
			sz = -zz;
			break;

		default:
			_EXCEPTIONT("Invalid destination panel id");
	}

	// (DBG) Check the free z coordinate
	if (sz <= 0.0) {
		_EXCEPTIONT("Invalid relative Z coordinate");
	}

	// Use panel information to calculate (alpha, beta) coords
	dX_out = sx / sz;
	dY_out  = sy / sz;
}

////////////////////////////////////////////////////////////////////////////////

void CubedSphereTrans::RLLFromXYP(
	double dX,
	double dY,
	int nP,
	double &lon,
	double &lat
) {
	switch (nP) {
		// Equatorial panel 1
		case 0:
			lon = atan(dX);
			lat = atan(dY / sqrt(1.0 + dX * dX));
			break;

		// Equatorial panel 2
		case 1:
			lon = atan(dX) + 0.5 * M_PI;
			lat = atan(dY / sqrt(1.0 + dX * dX));
			break;

		// Equatorial panel 3
		case 2:
			lon = atan(dX) + M_PI;
			lat = atan(dY / sqrt(1.0 + dX * dX));
			break;

		// Equatorial panel 4
		case 3:
			lon = atan(dX) + 1.5 * M_PI;
			lat = atan(dY / sqrt(1.0 + dX * dX));
			break;

		// North polar panel
		case 4:
			if (fabs(dX) > DBL_EPSILON) {
				lon = atan2(dX, -dY);
			} else if (dY <= 0.0) {
				lon = 0.0;
			} else {
				lon = M_PI;
			}
			lat = 0.5 * M_PI - atan(sqrt(dX * dX + dY * dY));
			break;

		// South polar panel
		case 5:
			if (fabs(dX) > DBL_EPSILON) {
				lon = atan2(dX, dY);
			} else if (dY > 0.0) {
				lon = 0.0;
			} else {
				lon = M_PI;
			}
			lat = -0.5 * M_PI + atan(sqrt(dX * dX + dY * dY));
			break;

		// Invalid panel
		default:
			_EXCEPTION1(
				"Invalid nP coordinate.  Given: %d, Expected: [0-5].\n", nP);
	}

	// Map to the interval [0, 2 pi]
	if (lon < 0.0) {
		lon += 2.0 * M_PI;
	}
}

////////////////////////////////////////////////////////////////////////////////

void CubedSphereTrans::RLLFromABP(
	double dA,
	double dB,
	int nP,
	double & lon,
	double & lat
) {
	return RLLFromXYP(tan(dA), tan(dB), nP, lon, lat);
}

////////////////////////////////////////////////////////////////////////////////

void CubedSphereTrans::XYPFromRLL(
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

void CubedSphereTrans::ABPFromRLL(
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

void CubedSphereTrans::VecTransABPFromRLL(
	double dX,
	double dY,
	int nP,
	double dUlon,
	double dUlat,
	double & dUalpha,
	double & dUbeta
) {
	double dDelta2 = 1.0 + dX * dX + dY * dY;
	double dRadius;

	double lat;

	if ((nP > 3) && (fabs(dX) < 1.0e-13) && (fabs(dY) < 1.0e-13)) {
		if (nP == 4) {
			dUalpha = dUlon;
		} else {
			dUalpha = - dUlon;
		}
		dUbeta = dUlat;
		return;
	}

	switch (nP) {
		// Equatorial panels
		case 0:
		case 1:
		case 2:
		case 3:
			// Convert spherical coords to geometric basis
			lat = atan(dY / sqrt(1.0 + dX * dX));
			dUlon = dUlon / cos(lat);

			// Calculate new vector components
			dUalpha = dUlon;
			dUbeta =
				dX * dY / (1.0 + dY * dY) * dUlon
				+ dDelta2 / ((1.0 + dY * dY) * sqrt(1.0 + dX * dX)) * dUlat;
			break;

		// North polar panel
		case 4:
			// Convert spherical coords to geometric basis
			lat = 0.5 * M_PI - atan(sqrt(dX * dX + dY * dY));
			dUlon = dUlon / cos(lat);

			// Calculate new vector components
			dRadius = sqrt(dX * dX + dY * dY);

			dUalpha =
				- dY / (1.0 + dX * dX) * dUlon
				- dDelta2 * dX / ((1.0 + dX * dX) * dRadius) * dUlat;

			dUbeta =
				dX / (1.0 + dY * dY) * dUlon
				- dDelta2 * dY / ((1.0 + dY * dY) * dRadius) * dUlat;
			break;

		// South polar panel
		case 5:
			// Convert spherical coords to geometric basis
			lat = -0.5 * M_PI + atan(sqrt(dX * dX + dY * dY));
			dUlon = dUlon / cos(lat);

			// Calculate new vector components
			dRadius = sqrt(dX * dX + dY * dY);

			dUalpha =
				dY / (1.0 + dX * dX) * dUlon
				+ dDelta2 * dX / ((1.0 + dX * dX) * dRadius) * dUlat;

			dUbeta =
				- dX / (1.0 + dY * dY) * dUlon
				+ dDelta2 * dY / ((1.0 + dY * dY) * dRadius) * dUlat;
			break;

		// Invalid panel
		default:
			_EXCEPTION1(
				"Invalid nP coordinate.  Given: %d, Expected: [0-5].\n", nP);
	}
}

////////////////////////////////////////////////////////////////////////////////

void CubedSphereTrans::VecTransRLLFromABP(
	double dX,
	double dY,
	int nP,
	double dUalpha,
	double dUbeta,
	double & dUlon,
	double & dUlat
) {
	double dDelta2 = 1.0 + dX * dX + dY * dY;
	double dRadius;

	double lat;

	switch (nP) {
		// Equatorial panels
		case 0:
		case 1:
		case 2:
		case 3:
			// Calculate new vector components
			dUlon = dUalpha;
			dUlat = 
				- dX * dY * sqrt(1.0 + dX * dX) / dDelta2 * dUalpha
				+ (1.0 + dY * dY) * sqrt(1.0 + dX * dX) / dDelta2 * dUbeta;

			// Convert spherical coords to unit basis
			lat = atan(dY / sqrt(1.0 + dX * dX));
			dUlon = dUlon * cos(lat);
			break;

		// North polar panel
		case 4:
			// Calculate new vector components
			dRadius = sqrt(dX * dX + dY * dY);

			dUlon = 
				- dY * (1.0 + dX * dX) / (dRadius * dRadius) * dUalpha
				+ dX * (1.0 + dY * dY) / (dRadius * dRadius) * dUbeta;

			dUlat =
				- dX * (1.0 + dX * dX) / (dDelta2 * dRadius) * dUalpha
				- dY * (1.0 + dY * dY) / (dDelta2 * dRadius) * dUbeta;

			// Convert spherical coords to unit basis
			lat = 0.5 * M_PI - atan(sqrt(dX * dX + dY * dY));
			dUlon = dUlon * cos(lat);

			break;

		// South polar panel
		case 5:
			// Calculate new vector components
			dRadius = sqrt(dX * dX + dY * dY);

			dUlon =
				dY * (1.0 + dX * dX) / (dRadius * dRadius) * dUalpha
				- dX * (1.0 + dY * dY) / (dRadius * dRadius) * dUbeta;

			dUlat =
				dX * (1.0 + dX * dX) / (dDelta2 * dRadius) * dUalpha
				+ dY * (1.0 + dY * dY) / (dDelta2 * dRadius) * dUbeta;
		
			// Convert spherical coords to unit basis
			lat = -0.5 * M_PI + atan(sqrt(dX * dX + dY * dY));
			dUlon = dUlon * cos(lat);

			break;

		// Invalid panel
		default:
			_EXCEPTION1(
				"Invalid nP coordinate.  Given: %d, Expected: [0-5].\n", nP);
	}
}

////////////////////////////////////////////////////////////////////////////////

void CubedSphereTrans::CoVecTransABPFromRLL(
	double dX,
	double dY,
	int nP,
	double dUlon,
	double dUlat,
	double & dUalpha,
	double & dUbeta
) {
	double dDelta2 = 1.0 + dX * dX + dY * dY;
	double dRadius;

	double lat;

	if ((nP > 3) && (fabs(dX) < 1.0e-13) && (fabs(dY) < 1.0e-13)) {
		if (nP == 4) {
			dUalpha = dUlon;
		} else {
			dUalpha = - dUlon;
		}
		dUbeta = dUlat;
		return;
	}

	switch (nP) {
		// Equatorial panels
		case 0:
		case 1:
		case 2:
		case 3:
			// Convert Ulon from geometric basis
			lat = atan(dY / sqrt(1.0 + dX * dX));
			dUlon = dUlon / cos(lat);

			// Calculate vector component
			dUalpha =
				  (1.0 + dX * dX) / dDelta2 * dUlon
				- dX * dY * sqrt(1.0 + dX * dX) / dDelta2 * dUlat;
			dUbeta = 
				sqrt(1.0 + dX * dX) * (1.0 + dY * dY) / dDelta2 * dUlat;

			break;

		// North polar panel
		case 4:
			// Convert Ulon from geometric basis
			lat = 0.5 * M_PI - atan(sqrt(dX * dX + dY * dY));
			dUlon = dUlon / cos(lat);

			// Calculate vector component
			dRadius = sqrt(dX * dX + dY * dY);

			dUalpha = 
				- dY * (1.0 + dX * dX) / dDelta2 * dUlon
				- dX * (1.0 + dX * dX) / (dDelta2 * dRadius) * dUlat;

			dUbeta =
				+ dX * (1.0 + dY * dY) / dDelta2 * dUlon
				- dY * (1.0 + dY * dY) / (dDelta2 * dRadius) * dUlat;

			break;

		// South polar panel
		case 5:

			// Convert Ulon from geometric basis
			lat = -0.5 * M_PI + atan(sqrt(dX * dX + dY * dY));
			dUlon = dUlon / cos(lat);

			// Calculate vector component
			dRadius = sqrt(dX * dX + dY * dY);

			dUalpha = 
				+ dY * (1.0 + dX * dX) / dDelta2 * dUlon
				+ dX * (1.0 + dX * dX) / (dDelta2 * dRadius) * dUlat;

			dUbeta =
				- dX * (1.0 + dY * dY) / dDelta2 * dUlon
				+ dY * (1.0 + dY * dY) / (dDelta2 * dRadius) * dUlat;

			break;

		// Invalid panel
		default:
			_EXCEPTION1(
				"Invalid nP coordinate.  Given: %d, Expected: [0-5].\n", nP);
	}
}

////////////////////////////////////////////////////////////////////////////////

void CubedSphereTrans::CoVecTransRLLFromABP(
	double dX,
	double dY,
	int nP,
	double dUalpha,
	double dUbeta,
	double & dUlon,
	double & dUlat
) {
	double dDelta2 = 1.0 + dX * dX + dY * dY;
	double dRadius;
	double dRadius2;

	double lat;

	if ((nP > 3) && (fabs(dX) < 1.0e-13) && (fabs(dY) < 1.0e-13)) {
		if (nP == 4) {
			dUlon = dUalpha;
		} else {
			dUlon = - dUalpha;
		}
		dUlat = dUbeta;
		return;
	}

	switch (nP) {
		// Equatorial panels
		case 0:
		case 1:
		case 2:
		case 3:
			// Calculate new vector components
			dUlon =
				  dDelta2 / (1.0 + dX * dX) * dUalpha
				+ dDelta2 * dX * dY / (1.0 + dX * dX) / (1.0 + dY * dY) * dUbeta;
			dUlat =
				  dDelta2 / sqrt(1.0 + dX * dX) / (1.0 + dY * dY) * dUbeta;

			// Convert spherical coords to geometric basis
			lat = atan(dY / sqrt(1.0 + dX * dX));
			dUlon *= cos(lat);

			break;

		// North polar panel
		case 4:
			// Calculate new vector components
			dRadius2 = (dX * dX + dY * dY);
			dRadius = sqrt(dRadius2);

			dUlon =
				- dDelta2 * dY / (1.0 + dX * dX) / dRadius2 * dUalpha
				+ dDelta2 * dX / (1.0 + dY * dY) / dRadius2 * dUbeta;

			dUlat =
				- dDelta2 * dX / (1.0 + dX * dX) / dRadius * dUalpha
				- dDelta2 * dY / (1.0 + dY * dY) / dRadius * dUbeta;

			// Convert spherical coords to geometric basis
			lat = 0.5 * M_PI - atan(sqrt(dX * dX + dY * dY));
			dUlon *= cos(lat);

			break;

		// South polar panel
		case 5:
			// Calculate new vector components
			dRadius2 = (dX * dX + dY * dY);
			dRadius = sqrt(dRadius2);

			dUlon =
				+ dDelta2 * dY / (1.0 + dX * dX) / dRadius2 * dUalpha
				- dDelta2 * dX / (1.0 + dY * dY) / dRadius2 * dUbeta;

			dUlat =
				+ dDelta2 * dX / (1.0 + dX * dX) / dRadius * dUalpha
				+ dDelta2 * dY / (1.0 + dY * dY) / dRadius * dUbeta;

			// Convert spherical coords to geometric basis
			lat = 0.5 * M_PI - atan(sqrt(dX * dX + dY * dY));
			dUlon *= cos(lat);

			break;

		// Invalid panel
		default:
			_EXCEPTION1(
				"Invalid nP coordinate.  Given: %d, Expected: [0-5].\n", nP);
	}
}

////////////////////////////////////////////////////////////////////////////////

