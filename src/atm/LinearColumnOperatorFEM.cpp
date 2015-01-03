///////////////////////////////////////////////////////////////////////////////
///
///	\file    LinearColumnOperatorFEM.cpp
///	\author  Paul Ullrich
///	\version July 29, 2014
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

#include "LinearColumnOperatorFEM.h"

#include "FluxCorrectionFunction.h"

#include "PolynomialInterp.h"

#include <cmath>

///////////////////////////////////////////////////////////////////////////////
/// LinearColumnInterpFEM
///////////////////////////////////////////////////////////////////////////////

void LinearColumnInterpFEM::Initialize(
	InterpSource eInterpSource,
	int nVerticalOrder,
	const DataVector<double> & dREtaNode,
	const DataVector<double> & dREtaREdge,
	const DataVector<double> & dREtaOut
) {
	const double ParamEpsilon = 1.0e-12;

	const int nRElementsIn  = dREtaNode.GetRows();
	const int nRElementsOut = dREtaOut.GetRows();

	const int nFiniteElementsIn = nRElementsIn / nVerticalOrder;

	// Verify input parameters
	if (nRElementsIn == 0) {
		_EXCEPTIONT("At least one row required for dREtaNode");
	}
	if (dREtaREdge.GetRows() != nRElementsIn + 1) {
		_EXCEPTIONT("REtaNode / REtaREdge mismatch");
	}
	if (nRElementsIn % nVerticalOrder != 0) {
		_EXCEPTIONT("Column RElements / VerticalOrder mismatch");
	}

	// Initialize LinearColumnOperator
	if (eInterpSource == InterpSource_Interfaces) {
		LinearColumnOperator::Initialize(nRElementsIn+1, nRElementsOut);
	} else {
		LinearColumnOperator::Initialize(nRElementsIn, nRElementsOut);
	}

	// Loop through all output elements
	for (int l = 0; l < nRElementsOut; l++) {

		// Determine input element index
		bool fOnREdge = false;
		int a;

		for (a = 0; a < nFiniteElementsIn - 1; a++) {
			double dNextREtaREdgeIn =
				dREtaREdge[(a+1) * nVerticalOrder] - ParamEpsilon;

			if (dREtaOut[l] < dNextREtaREdgeIn) {
				break;
			}
			if (dREtaOut[l] < dNextREtaREdgeIn + 2.0 * ParamEpsilon) {
				fOnREdge = true;
				break;
			}
		}

		// Interpolation coefficients for a continuous basis
		if (eInterpSource == InterpSource_Interfaces) {

			if (fOnREdge) {
				m_dCoeff[l][(a+1) * nVerticalOrder] = 1.0;

			} else {
				PolynomialInterp::LagrangianPolynomialCoeffs(
					nVerticalOrder + 1,
					&(dREtaREdge[a * nVerticalOrder]),
					&(m_dCoeff[l][a * nVerticalOrder]),
					dREtaOut[l]);

				m_iBegin[l] =  a    * nVerticalOrder;
				m_iEnd[l]   = (a+1) * nVerticalOrder + 1;
			}

		// Interpolation coefficients for a discontinuous basis
		} else {

			PolynomialInterp::LagrangianPolynomialCoeffs(
				nVerticalOrder,
				&(dREtaNode[a * nVerticalOrder]),
				&(m_dCoeff[l][a * nVerticalOrder]),
				dREtaOut[l]);

			m_iBegin[l] =  a    * nVerticalOrder;
			m_iEnd[l]   = (a+1) * nVerticalOrder;

			// Interpolating from nodes to edges, weight left and right
			// interpolant to minimize error.  Note that this approach does not
			// guarantee minimum error for sub-element stretching of REta.
			if (fOnREdge) {

				// Calculate one-sided errors and interpolant weights
				double dDeltaREtaL =
					  dREtaREdge[(a+1) * nVerticalOrder]
					- dREtaREdge[ a    * nVerticalOrder];

				double dDeltaREtaR =
					  dREtaREdge[(a+2) * nVerticalOrder]
					- dREtaREdge[(a+1) * nVerticalOrder];

				double dErrorL =
					pow(dDeltaREtaL, static_cast<double>(nVerticalOrder));

				double dErrorR =
					pow(dDeltaREtaR, static_cast<double>(nVerticalOrder));

				double dWeightL = dErrorR / (dErrorL + dErrorR);
				double dWeightR = dErrorL / (dErrorL + dErrorR);

				// Right interpolant coefficients
				PolynomialInterp::LagrangianPolynomialCoeffs(
					nVerticalOrder,
					&(dREtaNode[(a+1) * nVerticalOrder]),
					&(m_dCoeff[l][m_iEnd[l]]),
					dREtaOut[l]);

				// Weight interpolants by error
				for (int k = m_iBegin[l]; k < m_iEnd[l]; k++) {
					m_dCoeff[l][k] *= dWeightL;
				}

				int iNewEnd = m_iEnd[l] + nVerticalOrder;
				for (int k = m_iEnd[l]; k < iNewEnd; k++) {
					m_dCoeff[l][k] *= dWeightR;
				}

				m_iEnd[l] = iNewEnd;
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
/// LinearColumnDiffFEM
///////////////////////////////////////////////////////////////////////////////

void LinearColumnDiffFEM::Initialize(
	InterpSource eInterpSource,
	int nVerticalOrder,
	const DataVector<double> & dREtaNode,
	const DataVector<double> & dREtaREdge,
	const DataVector<double> & dREtaOut,
	bool fZeroBoundaries
) {
	const int ParamFluxCorrectionType = 2;

	const double ParamEpsilon = 1.0e-12;

	const int nRElementsIn  = dREtaNode.GetRows();
	const int nRElementsOut = dREtaOut.GetRows();

	const int nFiniteElements = nRElementsIn / nVerticalOrder;

	// Verify input parameters
	if (nRElementsIn == 0) {
		_EXCEPTIONT("At least one row required for dREtaNode");
	}
	if (dREtaREdge.GetRows() != nRElementsIn + 1) {
		_EXCEPTIONT("REtaNode / REtaREdge mismatch");
	}
	if (nRElementsIn % nVerticalOrder != 0) {
		_EXCEPTIONT("Column RElements / VerticalOrder mismatch");
	}

	// Initialize LinearColumnOperator
	if (eInterpSource == InterpSource_Interfaces) {
		LinearColumnOperator::Initialize(nRElementsIn+1, nRElementsOut);
	} else {
		LinearColumnOperator::Initialize(nRElementsIn, nRElementsOut);
	}

	// Loop through all output elements
	for (int l = 0; l < nRElementsOut; l++) {

		// Determine input element index and whether we are on a finite
		// element edge (excluding top and bottom boundary)
		bool fOnREdge = false;
		int a;

		for (a = 0; a < nFiniteElements - 1; a++) {
			double dNextREtaREdgeIn =
				dREtaREdge[(a+1) * nVerticalOrder] - ParamEpsilon;

			if (dREtaOut[l] < dNextREtaREdgeIn) {
				break;
			}
			if (dREtaOut[l] < dNextREtaREdgeIn + 2.0 * ParamEpsilon) {
				fOnREdge = true;
				break;
			}
		}

		// Differentiation coefficients for a continuous basis
		if (eInterpSource == InterpSource_Interfaces) {

			PolynomialInterp::DiffLagrangianPolynomialCoeffs(
				nVerticalOrder + 1,
				&(dREtaREdge[a * nVerticalOrder]),
				&(m_dCoeff[l][a * nVerticalOrder]),
				dREtaOut[l]);

			// Set bounds on coefficients
			if (!fOnREdge) {
				m_iBegin[l] =  a * nVerticalOrder;
				m_iEnd[l]   = (a+1) * nVerticalOrder + 1;
			}

			// Special treatment of derivatives at interfaces
			if (fOnREdge) {

				// Temporary coefficients
				DataVector<double> dTempCoeff;
				dTempCoeff.Initialize(nVerticalOrder + 1);
				
				// Calculate one-sided errors and derivative weights
				double dDeltaREtaL =
					  dREtaREdge[(a+1) * nVerticalOrder]
					- dREtaREdge[ a    * nVerticalOrder];

				double dDeltaREtaR =
					  dREtaREdge[(a+2) * nVerticalOrder]
					- dREtaREdge[(a+1) * nVerticalOrder];

				double dErrorL =
					pow(dDeltaREtaL, static_cast<double>(nVerticalOrder));

				double dErrorR =
					pow(dDeltaREtaR, static_cast<double>(nVerticalOrder));

				double dWeightL = dErrorR / (dErrorL + dErrorR);
				double dWeightR = dErrorL / (dErrorL + dErrorR);

				// Calculate right-side derivative coefficients
				PolynomialInterp::DiffLagrangianPolynomialCoeffs(
					nVerticalOrder + 1,
					&(dREtaREdge[(a+1) * nVerticalOrder]),
					&(dTempCoeff[0]),
					dREtaOut[l]);

				for (int k = 0; k <= nVerticalOrder; k++) {
					m_dCoeff[l][a * nVerticalOrder + k] *= dWeightL;
				}

				for (int k = 0; k <= nVerticalOrder; k++) {
					m_dCoeff[l][(a+1) * nVerticalOrder + k] +=
						dWeightR * dTempCoeff[k];
				}

				// Set bounds on coefficients
				m_iBegin[l] =  a * nVerticalOrder;
				m_iEnd[l]   = (a+2) * nVerticalOrder + 1;
			}

		// Differentiation coefficients for a discontinuous basis
		} else {

			if ((dREtaOut[l] < dREtaREdge[0]) ||
				(dREtaOut[l] > dREtaREdge[nRElementsIn])
			) {
				_EXCEPTIONT("REta out of range");
			}

			// Element spacing
			double dDeltaREta = 
				dREtaREdge[(a+1) * nVerticalOrder]
				- dREtaREdge[a * nVerticalOrder];

			// Contribution due to local derivative
			PolynomialInterp::DiffLagrangianPolynomialCoeffs(
				nVerticalOrder,
				&(dREtaNode[a * nVerticalOrder]),
				&(m_dCoeff[l][a * nVerticalOrder]),
				dREtaOut[l]);

			if (fOnREdge) {
				PolynomialInterp::DiffLagrangianPolynomialCoeffs(
					nVerticalOrder,
					&(dREtaNode[(a+1) * nVerticalOrder]),
					&(m_dCoeff[l][(a+1) * nVerticalOrder]),
					dREtaOut[l]);

				for (int k = 0; k < nRElementsIn; k++) {
					m_dCoeff[l][k] *= 0.5 * dDeltaREta;
				}

			} else {
				for (int k = 0; k < nRElementsIn; k++) {
					m_dCoeff[l][k] *= dDeltaREta;
				}
			}

			// Interpolation coefficients to finite element edges
			DataVector<double> dTempCoeffLL;
			dTempCoeffLL.Initialize(nVerticalOrder);

			DataVector<double> dTempCoeffLR;
			dTempCoeffLR.Initialize(nVerticalOrder);

			DataVector<double> dTempCoeffRL;
			dTempCoeffRL.Initialize(nVerticalOrder);

			DataVector<double> dTempCoeffRR;
			dTempCoeffRR.Initialize(nVerticalOrder);

			// Derivatives of the flux reconstruction function at this point
			DataVector<double> dNodesR;
			dNodesR.Initialize(1);
			dNodesR[0] =
				(dREtaOut[l] - dREtaREdge[a * nVerticalOrder]) / dDeltaREta;

			DataVector<double> dNodesL;
			dNodesL.Initialize(1);
			dNodesL[0] = 1.0 - dNodesR[0];

			DataVector<double> dDerivR;
			dDerivR.Initialize(1);

			DataVector<double> dDerivL;
			dDerivL.Initialize(1);

			FluxCorrectionFunction::GetDerivatives(
				ParamFluxCorrectionType,
				nVerticalOrder + 1,
				dNodesR,
				dDerivR);

			FluxCorrectionFunction::GetDerivatives(
				ParamFluxCorrectionType,
				nVerticalOrder + 1,
				dNodesL,
				dDerivL);

			dDerivL[0] *= -1.0;

			// Interpolation coefficients to element interfaces
			PolynomialInterp::LagrangianPolynomialCoeffs(
				nVerticalOrder,
				&(dREtaNode[a * nVerticalOrder]),
				&(dTempCoeffLR[0]),
				dREtaREdge[a * nVerticalOrder]);

			PolynomialInterp::LagrangianPolynomialCoeffs(
				nVerticalOrder,
				&(dREtaNode[a * nVerticalOrder]),
				&(dTempCoeffRL[0]),
				dREtaREdge[(a+1) * nVerticalOrder]);

			if (a != 0) {
				PolynomialInterp::LagrangianPolynomialCoeffs(
					nVerticalOrder,
					&(dREtaNode[(a-1) * nVerticalOrder]),
					&(dTempCoeffLL[0]),
					dREtaREdge[a * nVerticalOrder]);
			}

			if (a != nFiniteElements - 1) {
				PolynomialInterp::LagrangianPolynomialCoeffs(
					nVerticalOrder,
					&(dREtaNode[(a+1) * nVerticalOrder]),
					&(dTempCoeffRR[0]),
					dREtaREdge[(a+1) * nVerticalOrder]);
			}

			// Apply edge interpolation coefficients
			if (a != 0) {
				if (!fOnREdge) {
					for (int k = 0; k < nVerticalOrder; k++) {
						m_dCoeff[l][(a-1) * nVerticalOrder + k] +=
							0.5 * dDerivL[0] * dTempCoeffLL[k];
					}
				}

				for (int k = 0; k < nVerticalOrder; k++) {
					m_dCoeff[l][a * nVerticalOrder + k] -=
						0.5 * dDerivL[0] * dTempCoeffLR[k];
				}

			} else {
				if (!fZeroBoundaries) {
					for (int k = 0; k < nVerticalOrder; k++) {
						m_dCoeff[l][ a    * nVerticalOrder + k] +=
							0.5 * dDerivL[0] * dTempCoeffRL[k];
						m_dCoeff[l][(a+1) * nVerticalOrder + k] -=
							0.5 * dDerivL[0] * dTempCoeffRR[k];
					}
				}
			}

			if (a != nFiniteElements - 1) {
				for (int k = 0; k < nVerticalOrder; k++) {
					m_dCoeff[l][(a+1) * nVerticalOrder + k] +=
						0.5 * dDerivR[0] * dTempCoeffRR[k];
				}
				for (int k = 0; k < nVerticalOrder; k++) {
					m_dCoeff[l][a * nVerticalOrder + k] -=
						0.5 * dDerivR[0] * dTempCoeffRL[k];
				}

			} else {
				if (!fZeroBoundaries) {
					for (int k = 0; k < nVerticalOrder; k++) {
						m_dCoeff[l][ a    * nVerticalOrder + k] +=
							0.5 * dDerivR[0] * dTempCoeffLR[k];
						m_dCoeff[l][(a-1) * nVerticalOrder + k] -=
							0.5 * dDerivR[0] * dTempCoeffLL[k];
					}
				}
			}

			for (int k = 0; k < nRElementsIn; k++) {
				m_dCoeff[l][k] /= dDeltaREta;
			}

			// Set boundaries
			if (a != 0) {
				m_iBegin[l] = (a-1) * nVerticalOrder;
			} else {
				m_iBegin[l] =  a    * nVerticalOrder;
			}

			if (a != nFiniteElements-1) {
				m_iEnd[l]   = (a+2) * nVerticalOrder;
			} else {
				m_iEnd[l]   = (a+1) * nVerticalOrder;
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

