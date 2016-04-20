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

LinearColumnOperator::LinearColumnOperator() :
	m_fInitialized(false)
{
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
	m_fInitialized = true;
	m_dCoeff.Allocate(nRElementsOut, nRElementsIn);
	m_iBegin.Allocate(nRElementsOut);
	m_iEnd  .Allocate(nRElementsOut);
}

///////////////////////////////////////////////////////////////////////////////

void LinearColumnOperator::InitializeIdentity(
	int nRElements
) {
	m_fInitialized = true;
	m_dCoeff.Allocate(nRElements, nRElements);
	m_iBegin.Allocate(nRElements);
	m_iEnd  .Allocate(nRElements);

	for (int i = 0; i < nRElements; i++) {
		m_dCoeff[i][i] = 1.0;
		m_iBegin[i] = i;
		m_iEnd[i] = i+1;
	}
}

///////////////////////////////////////////////////////////////////////////////

void LinearColumnOperator::ComposeWith(
	const LinearColumnOperator & op
) {
	// Verify dimensions
	if (m_dCoeff.GetColumns() != op.m_dCoeff.GetRows()) {
		_EXCEPTION2("Composed operator dimension mismatch (%i, %i)",
			m_dCoeff.GetColumns(), op.m_dCoeff.GetRows());
	}

	// Intermediate size
	int nRElementsInter = m_dCoeff.GetColumns();

	// Construct composed operator
	int nRElementsOut = m_dCoeff.GetRows();
	int nRElementsIn  = op.m_dCoeff.GetColumns();

	if (nRElementsOut == 0) {
		_EXCEPTIONT("Invalid number of rows in source operator");
	}
	if (nRElementsIn == 0) {
		_EXCEPTIONT("Invalid number of columns in target operator");
	}

	DataArray2D<double> dNewCoeff(nRElementsOut, nRElementsIn);
	DataArray1D<int> iNewBegin(nRElementsOut);
	DataArray1D<int> iNewEnd(nRElementsOut);

	// Perform matrix product
	for (int i = 0; i < nRElementsOut; i++) {
	for (int j = 0; j < nRElementsIn; j++) {
		for (int k = 0; k < nRElementsInter; k++) {
			dNewCoeff[i][j] += m_dCoeff[i][k] * op.m_dCoeff[k][j];
		}
	}
	}

	// Determine non-zero intervals
	for (int i = 0; i < nRElementsOut; i++) {
		bool fFoundBegin = false;
		for (int j = 0; j < nRElementsIn; j++) {
			if (dNewCoeff[i][j] != 0.0) {
				iNewEnd[i] = j + 1;
				if (!fFoundBegin) {
					fFoundBegin = true;
					iNewBegin[i] = j;
				}
			}
		}
	}

	// Overwrite current op
	m_dCoeff = dNewCoeff;
	m_iBegin = iNewBegin;
	m_iEnd   = iNewEnd;
}

///////////////////////////////////////////////////////////////////////////////

void LinearColumnOperator::DebugOutput(
	const DataArray1D<double> * pREtaNode,
	const DataArray1D<double> * pREtaREdge,
	const std::string & strTag,
	bool fExcept
) {
	std::string strOpFile = "op" + strTag + ".txt";
	std::string strRnFile = "rn" + strTag + ".txt";
	std::string strRiFile = "ri" + strTag + ".txt";

	FILE * fp = fopen(strOpFile.c_str(), "w");
	for (int i = 0; i < m_dCoeff.GetRows(); i++) {
		for (int j = 0; j < m_dCoeff.GetColumns(); j++) {
			fprintf(fp, "%1.15e\t", m_dCoeff[i][j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	if (pREtaNode != NULL) {
		FILE * fpn = fopen(strRnFile.c_str(), "w");
		for (int i = 0; i < pREtaNode->GetRows(); i++) {
			fprintf(fpn, "%1.15e\n", (*pREtaNode)[i]);
		}
		fclose(fpn);
	}

	if (pREtaREdge != NULL) {
		FILE * fpi = fopen(strRiFile.c_str(), "w");
		for (int i = 0; i < pREtaREdge->GetRows(); i++) {
			fprintf(fpi, "%1.15e\n", (*pREtaREdge)[i]);
		}
		fclose(fpi);
	}

	for (int i = 0; i < m_iBegin.GetRows(); i++) {
		printf("%i %i\n", m_iBegin[i], m_iEnd[i]);
	}

	if (fExcept) {
		_EXCEPTION();
	}
}

///////////////////////////////////////////////////////////////////////////////

double LinearColumnOperator::Apply(
	const double * dColumnIn,
	int iRout,
	int nStride
) const {
	if (!m_fInitialized) {
		_EXCEPTIONT("Attempting to Apply uninitialized LinearColumnOperator");
	}

	double dOut = 0.0;
	for (int l = m_iBegin[iRout]; l < m_iEnd[iRout]; l++) {
		int lx = l * nStride;

		dOut += m_dCoeff[iRout][l] * dColumnIn[lx];
	}
	return dOut;
}

///////////////////////////////////////////////////////////////////////////////

double LinearColumnOperator::Apply(
	const double * dColumnIn,
	const double * dColumnRefIn,
	double dColumnRefOut,
	int iRout,
	int nStride
) const {
	if (!m_fInitialized) {
		_EXCEPTIONT("Attempting to Apply uninitialized LinearColumnOperator");
	}

	double dOut = dColumnRefOut;
	for (int l = m_iBegin[iRout]; l < m_iEnd[iRout]; l++) {
		int lx = l * nStride;

		dOut += m_dCoeff[iRout][l] * (dColumnIn[lx] - dColumnRefIn[lx]);
	}
	return dOut;
}

///////////////////////////////////////////////////////////////////////////////

void LinearColumnOperator::Apply(
	const double * dColumnIn,
	double * dColumnOut
) const {
	if (!m_fInitialized) {
		_EXCEPTIONT("Attempting to Apply uninitialized LinearColumnOperator");
	}

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
	if (!m_fInitialized) {
		_EXCEPTIONT("Attempting to Apply uninitialized LinearColumnOperator");
	}

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

void LinearColumnOperator::ApplyWithRef(
	const double * dColumnIn,
	const double * dColumnRefIn,
	double * dColumnOut,
	const double * dColumnRefOut
) const {
	if (!m_fInitialized) {
		_EXCEPTIONT("Attempting to Apply uninitialized LinearColumnOperator");
	}

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

void LinearColumnOperator::ApplyWithRef(
	const double * dColumnIn,
	const double * dColumnRefIn,
	double * dColumnOut,
	const double * dColumnRefOut,
	int nStrideIn,
	int nStrideOut
) const {
	if (!m_fInitialized) {
		_EXCEPTIONT("Attempting to Apply uninitialized LinearColumnOperator");
	}

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


