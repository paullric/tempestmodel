///////////////////////////////////////////////////////////////////////////////
///
///	\file    MathHelper.cpp
///	\author  Paul Ullrich
///	\version July 26, 2010
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

#include "MathHelper.h"

///////////////////////////////////////////////////////////////////////////////

double LinearInterpolate(
	double dXout,
	double dX1,
	double dX2,
	double dY1,
	double dY2
) {
	double dW = (dXout - dX1) / (dX2 - dX1);

	return ((1.0 - dW) * dY1 + dW * dY2);
}

///////////////////////////////////////////////////////////////////////////////

double QuadraticInterpolate(
	double dXout,
	double dDX,
	double dY1,
	double dY2,
	double dY3
) {
	double dAdj = 1.0 / (2.0 * dDX * dDX);

	return dAdj * (
		+ (dY1 * (dXout - dDX) * (dXout - 2.0 * dDX))
		- (2.0 * dY2 * (dXout) * (dXout - 2.0 * dDX))
		+ (dY3 * (dXout) * (dXout - dDX)));
}

///////////////////////////////////////////////////////////////////////////////

double CubicInterpolate(
	double dXout,
	double dDX,
	double dY1,
	double dY2,
	double dY3,
	double dY4
) {
	double dAdj = 1.0 / (6.0 * dDX * dDX * dDX);

	double dXDX = dXout - dDX;
	double dX2DX = dXout - 2.0 * dDX;
	double dX3DX = dXout - 3.0 * dDX;

	return dAdj * (
		- (dY1 * dXDX * dX2DX * dX3DX)
		+ (dY2 * 3.0 * dXout * dX2DX * dX3DX)
		- (dY3 * 3.0 * dXout * dXDX * dX3DX)
		+ (dY4 * dXout * dXDX * dX2DX));
}

///////////////////////////////////////////////////////////////////////////////

