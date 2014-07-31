///////////////////////////////////////////////////////////////////////////////
///
///	\file    LinearColumnOperator.cpp
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

#include "LinearColumnOperator.h"

///////////////////////////////////////////////////////////////////////////////

LinearColumnOperator::LinearColumnOperator() {
}

///////////////////////////////////////////////////////////////////////////////

LinearColumnOperator::LinearColumnOperator(
	int nRElementsIn,
	int nRElementsOut
) {
	Initialize(nRElementsIn, nRElementsOut);
}

///////////////////////////////////////////////////////////////////////////////

void LinearColumnOperator::Initialize(
	int nRElementsIn,
	int nRElementsOut
) {
	m_dCoeff.Initialize(nRElementsOut, nRElementsIn);
	m_iBegin.Initialize(nRElementsOut);
	m_iEnd  .Initialize(nRElementsOut);
}

///////////////////////////////////////////////////////////////////////////////

void LinearColumnOperator::Apply(
	const double * dColumnIn,
	double * dColumnOut
) const {
	const int nRElementsOut = m_dCoeff.GetRows();

	for (int k = 0; k < nRElementsOut; k++) {
		dColumnOut[k] = 0.0;

		for (int l = m_iBegin[k]; l < m_iEnd[k]; l++) {
			dColumnOut[k] += m_dCoeff[k][l] * dColumnIn[l];
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void LinearColumnOperator::Apply(
	const double * dColumnIn,
	double * dColumnOut,
	int nStrideIn,
	int nStrideOut
) const {
	const int nRElementsOut = m_dCoeff.GetRows();

	for (int k = 0; k < nRElementsOut; k++) {
		int kx = k * nStrideOut;

		dColumnOut[kx] = 0.0;

		for (int l = m_iBegin[k]; l < m_iEnd[k]; l++) {
			int lx = l * nStrideIn;

			dColumnOut[kx] += m_dCoeff[k][l] * dColumnIn[lx];
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void LinearColumnOperator::Apply(
	const double * dColumnIn,
	const double * dColumnRefIn,
	double * dColumnOut,
	const double * dColumnRefOut
) const {
	const int nRElementsOut = m_dCoeff.GetRows();

	for (int k = 0; k < nRElementsOut; k++) {
		dColumnOut[k] = dColumnRefOut[k];

		for (int l = m_iBegin[k]; l < m_iEnd[k]; l++) {
			dColumnOut[k] +=
				m_dCoeff[k][l] * (dColumnIn[l] - dColumnRefIn[l]);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void LinearColumnOperator::Apply(
	const double * dColumnIn,
	const double * dColumnRefIn,
	double * dColumnOut,
	const double * dColumnRefOut,
	int nStrideIn,
	int nStrideOut
) const {
	const int nRElementsOut = m_dCoeff.GetRows();

	for (int k = 0; k < nRElementsOut; k++) {
		int kx = k * nStrideOut;

		dColumnOut[kx] = dColumnRefOut[kx];

		for (int l = m_iBegin[k]; l < m_iEnd[k]; l++) {
			int lx = l * nStrideIn;

			dColumnOut[kx] +=
				m_dCoeff[k][l] * (dColumnIn[lx] - dColumnRefIn[lx]);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////


