///////////////////////////////////////////////////////////////////////////////
///
///	\file    GaussianQuadrature.cpp
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

#include "GaussianQuadrature.h"

#include "Exception.h"

///////////////////////////////////////////////////////////////////////////////

int GaussianQuadrature::GetGaussPointCount(
	int nOrder
) {
	return ((nOrder - 1) / 2 + 1);
}

///////////////////////////////////////////////////////////////////////////////

void GaussianQuadrature::GetGaussPointsWeights(
	int nOrder,
	double * dG,
	double * dW
) {

	// Second-order scheme
	if (nOrder <= 2) {
		dG[0] = 0.0;

		dW[0] = 2.0;

	// Fourth-order scheme
	} else if (nOrder <= 4) {
		dG[0] = - 0.57735026918962576451;
		dG[1] = + 0.57735026918962576451;

		dW[0] = 1.0;
		dW[1] = 1.0;

	// Sixth-order scheme
	} else if (nOrder <= 6) {
		dG[0] = - 0.77459666924148337704;
		dG[1] =   0.0;
		dG[2] = + 0.77459666924148337704;

		dW[0] = 0.55555555555555555556;
		dW[1] = 0.88888888888888888889;
		dW[2] = 0.55555555555555555556;
 
 	// Eigth-order scheme
	} else if (nOrder <= 8) {
		dG[0] = - 0.86113631159405257524;
		dG[1] = - 0.33998104358485626481;
		dG[2] = + 0.33998104358485626481;
		dG[3] = + 0.86113631159405257524;

		dW[0] = 0.34785484513745385737;
		dW[1] = 0.65214515486254614263;
		dW[2] = 0.65214515486254614263;
		dW[3] = 0.34785484513745385737;

	// Tenth-order scheme
	} else if (nOrder <= 10) {
		dG[0] = - 0.90617984593866399282;
		dG[1] = - 0.53846931010568309105;
		dG[2] =   0.0;
		dG[3] = + 0.53846931010568309105;
		dG[4] = + 0.90617984593866399282;

		dW[0] = 0.23692688505618908752;
		dW[1] = 0.47862867049936646804;
		dW[2] = 0.56888888888888888889;
		dW[3] = 0.47862867049936646804;
		dW[4] = 0.23692688505618908752;

	// Unimplemented
	} else {
		_EXCEPTIONT("Unimplemented");
	}
}

///////////////////////////////////////////////////////////////////////////////

