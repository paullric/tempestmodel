///////////////////////////////////////////////////////////////////////////////
///
///	\file    FluxCorrectionFunction.cpp
///	\author  Paul Ullrich
///	\version September 11, 2013
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

#include "FluxCorrectionFunction.h"

#include "DataArray1D.h"
#include "DataArray2D.h"
#include "LinearAlgebra.h"

///////////////////////////////////////////////////////////////////////////////

void FluxCorrectionFunction::GetDerivatives(
	int iType,
	int nOrder,
	const DataArray1D<double> & dNodes,
	DataArray1D<double> & dDeriv
) {
	// Verify parameters
	if (iType < 1) {
		_EXCEPTIONT("iType must be at least 1");
	}
	if (nOrder < 1) {
		_EXCEPTIONT("Order must be at least 1");
	}

	// Compute edge reconstruction function (DG type)
	DataArray1D<double> dB(nOrder+1);

	DataArray2D<double> dVan(nOrder+1, nOrder+1);

	// Left value of reconstruction function is 1
	double dSign = 1.0;
	for (int i = nOrder; i >= 0; i--) {
		dVan[i][0] = dSign;
		dSign *= (-1.0);
	}
	dB[0] = 1.0;

	// Right-value has a zero of multiplicity = type
	DataArray1D<double> dCoeff(nOrder+1);
	for (int i = 0; i <= nOrder; i++) {
		dCoeff[i] = 1.0;
	}
	for (int n = 0; n < iType; n++) {
		for (int i = 0; i <= nOrder; i++) {
			dVan[i][n+1] = dCoeff[i];
		}
		for (int i = 0; i < nOrder - n; i++) {
			dCoeff[i] =
				static_cast<double>(nOrder - n - i) * dCoeff[i];
		}
		for (int i = nOrder - n; i <= nOrder; i++) {
			dCoeff[i] = 0.0;
		}
	}

	// Polynomial is perpendicular to all spaces of polynomials up to P_{type-2}
	for (int n = 0; n < nOrder - iType; n++) {
	for (int m = 0; m <= nOrder; m++) {
		int nSumCoeff = (nOrder - m + n);
		if ((nSumCoeff % 2) == 0) {
			dVan[m][iType+1+n] = 2.0 / static_cast<double>(nSumCoeff + 1.0);
		}
	}
	}

	// Solve Vandermonde system for polynomial coefficients
	DataArray1D<int> iPIV(nOrder+1);
	LAPACK::DGESV(dVan, dB, iPIV);

	// Compute derivative of reconstruction polynomial on [-1,1], with P(1)=1
	dSign = 1.0;
	for (int i = nOrder; i >= 0; i--) {
		dB[i] *= dSign;
		dSign *= - 1.0;
	}
	for (int i = 0; i < nOrder; i++) {
		dB[nOrder-i] =
			static_cast<double>(i+1) * dB[nOrder-i-1];
	}
	dB[0] = 0.0;

	// Compute derivatives of reconstruction polynomial on nodes for the
	// reference element [0,1]
	dDeriv.Allocate(dNodes.GetRows());
	for (int n = 0; n < dNodes.GetRows(); n++) {
		double dX = 1.0;
		for (int i = 0; i < nOrder; i++) {
			dDeriv[n] += dB[nOrder-i] * dX;
			dX *= (2.0 * dNodes[n] - 1.0);
		}
		dDeriv[n] *= 2.0;
	}
}

///////////////////////////////////////////////////////////////////////////////

