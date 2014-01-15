///////////////////////////////////////////////////////////////////////////////
///
///	\file    LinearAlgebra.cpp
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

#include "LinearAlgebra.h"

#include "Exception.h"

#include <iostream>

//////////////////////////////////////////////////////////////////////////////

void LAPACK::DGEMM(
	DataMatrix<double> & dC,
	DataMatrix<double> & dA,
	DataMatrix<double> & dB,
	double dAlpha,
	double dBeta
) {
	// Check dimensions
	if ((dA.GetRows() != dB.GetColumns()) ||
		(dA.GetColumns() != dC.GetColumns()) ||
		(dB.GetRows() != dC.GetRows())
	) {
		_EXCEPTIONT("Invalid matrix / matrix dimensions");
	}

	// Store CLAPACK parameters
	char cTransA = 'N';
	char cTransB = 'N';

	int m     = dA.GetColumns();
	int n     = dB.GetRows();
	int k     = dA.GetRows();

	int nLDA  = m;
	int nLDB  = k;
	int nLDC  = m;

#ifdef USEACML
	// Call the matrix solve
	dgemm(cTransA, cTransB, m, n, k,
		dAlpha, &(dA[0][0]), nLDA, &(dB[0][0]), nLDB,
		dBeta, &(dC[0][0]), nLDC);
#endif
#ifdef USEESSL
	// Call the matrix solve
	dgemm(cTransA, cTransB, m, n, k,
		dAlpha, &(dA[0][0]), nLDA, &(dB[0][0]), nLDB,
		dBeta, &(dC[0][0]), nLDC);
#endif
#if defined USEVECLIB || defined USEMKL
	// Call the matrix solve
	dgemm_(&cTransA, &cTransB, &m, &n, &k, &dAlpha,
		&(dA[0][0]), &nLDA, &(dB[0][0]), &nLDB, &dBeta, &(dC[0][0]), &nLDC);
#endif
}

//////////////////////////////////////////////////////////////////////////////

int LAPACK::DGESV(
	DataMatrix<double> & dA,
	DataVector<double> & dBX,
	DataVector<int> & iPIV
) {
	// Check dimensions
	if ((dA.GetRows() != dA.GetColumns()) ||
		(dA.GetColumns() != dBX.GetRows())
	) {
		_EXCEPTIONT("Invalid matrix / vector dimensions");
	}

	// Store CLAPACK parameters
	int n     = dA.GetRows();
	int nRHS  = 1;
	int nLDA  = n;
	int nLDB  = n;

	int nInfo;

#ifdef USEACML
	// Call the matrix solve
	dgesv(n, nRHS, &(dA[0][0]), nLDA, &(iPIV[0]), &(dBX[0]), nLDB, &nInfo);
#endif
#ifdef USEESSL
	// Call the matrix solve
	dgesv(n, nRHS, &(dA[0][0]), nLDA, &(iPIV[0]), &(dBX[0]), nLDB, nInfo);
#endif
#if defined USEVECLIB || defined USEMKL
	// Call the matrix solve
	dgesv_(
		&n, &nRHS,
		&(dA[0][0]), &nLDA, &(iPIV[0]), &(dBX[0]), &nLDB, &nInfo);
#endif

	return nInfo;
}

//////////////////////////////////////////////////////////////////////////////

int LAPACK::DGESV(
	DataMatrix<double> & dA,
	DataMatrix<double> & dBX,
	DataVector<int> & iPIV
) {
	// Check dimensions
	if ((dA.GetRows() != dA.GetColumns()) ||
		(dA.GetColumns() != dBX.GetColumns())
	) {
		_EXCEPTIONT("Invalid matrix / vector dimensions");
	}

	// Store CLAPACK parameters
	int n     = dA.GetRows();
	int nRHS  = dBX.GetRows();
	int nLDA  = n;
	int nLDB  = n;

	int nInfo;

#ifdef USEACML
	// Call the matrix solve
	dgesv(n, nRHS,
		&(dA[0][0]), nLDA, &(iPIV[0]), &(dBX[0][0]), nLDB, &nInfo);
#endif
#ifdef USEESSL
	// Call the matrix solve
	dgesv(n, nRHS,
		&(dA[0][0]), nLDA, &(iPIV[0]), &(dBX[0][0]), nLDB, nInfo);
#endif
#if defined USEVECLIB || defined USEMKL
	// Call the matrix solve
	dgesv_(
		&n, &nRHS,
		&(dA[0][0]), &nLDA, &(iPIV[0]), &(dBX[0][0]), &nLDB, &nInfo);
#endif

	return nInfo;
}

///////////////////////////////////////////////////////////////////////////////

int LAPACK::DGBSV(
	DataMatrix<double> & dA,
	DataVector<double> & dBX,
	DataVector<int> & iPIV,
	int nKL,
	int nKU
) {
	// Check dimensions
	if ((dA.GetRows() != dA.GetColumns()) ||
		(dA.GetColumns() != dBX.GetRows())
	) {
		_EXCEPTIONT("Invalid matrix / vector dimensions");
	}

	// Store CLAPACK parameters
	int n     = dA.GetRows();
	int nRHS  = 1;
	int nLDAB = n;
	int nLDB  = n;

	int nInfo;

#ifdef USEACML
	// Call the banded diagonal matrix solve
	dgbsv(
		n, nKL, nKU, nRHS,
		&(dA[0][0]), nLDAB, &(iPIV[0]), &(dBX[0]), nLDB, &nInfo);
#endif
#ifdef USEESSL
	// Call the banded diagonal matrix solve
	dgbf(&(dA[0][0]), nLDAB, n, nKL, nKU, &(iPIV[0]));
	dgbs(&(dA[0][0]), nLDAB, n, nKL, nKU, &(iPIV[0]), &(dBX[0]));
	nInfo = 0;
#endif
#if defined USEVECLIB || defined USEMKL
	// Call the banded diagonal matrix solve
	dgbsv_(
		&n, &nKL, &nKU, &nRHS,
		&(dA[0][0]), &nLDAB, &(iPIV[0]), &(dBX[0]), &nLDB, &nInfo);
#endif

	return nInfo;
}

///////////////////////////////////////////////////////////////////////////////

int LAPACK::DTPSV(
	char chUpperLower,
	char chTrans,
	char chDiagonal,
	int nN,
	DataVector<double> & dA,
	DataVector<double> & dX
) {
	int iIncX = 1;

#if defined USEACML || defined USEESSL
	_EXCEPTION();
#endif
#if defined USEVECLIB || defined USEMKL
	// Call the triangular matrix solve
	dtpsv_(
		&chUpperLower, &chTrans, &chDiagonal,
		&nN, &(dA[0]), &(dX[0]), &iIncX);
#endif

	return 0;
}

///////////////////////////////////////////////////////////////////////////////

int LAPACK::DGETRF(
	DataMatrix<double> & dA,
	DataVector<int> & iPIV
) {
	int m = dA.GetRows();
	int n = dA.GetColumns();

	int lda = m;

	int nInfo;

#ifdef USEACML
	dgetrf(m, n, &(dA[0][0]), lda, &(iPIV[0]), &nInfo);
#endif
#ifdef USEESSL
	dgetrf(m, n, &(dA[0][0]), lda, &(iPIV[0]), nInfo);
#endif
#if defined USEVECLIB || defined USEMKL
	dgetrf_(&m, &n, &(dA[0][0]), &lda, &(iPIV[0]), &nInfo);
#endif

	return nInfo;
}

///////////////////////////////////////////////////////////////////////////////

int LAPACK::DGETRI(
	DataMatrix<double> & dA,
	DataVector<int> & iPIV,
	DataVector<double> & dWork
) {
	if (dA.GetRows() != dA.GetColumns()) {
		_EXCEPTIONT("Matrix A must be square.");
	}
	if (iPIV.GetRows() != dA.GetRows()) {
		_EXCEPTIONT("Dimension of PIV must be same as matrix dimension.");
	}

	int n = dA.GetRows();

	int nWork = dWork.GetRows();

	if (nWork < n) {
		_EXCEPTION2("Work matrix too small:  Expected %i, found %i.", n, nWork);
	}

	int nInfo;

#if defined USEACML || defined USEESSL
	//dgetri(n, &(dA[0][0]), n, &(iPIV[0]), &(dWork[0]), nWork, &nInfo);
	nInfo = 0;
	_EXCEPTIONT("Unimplemented.");
#endif
#if defined USEVECLIB || defined USEMKL
	dgetri_(&n, &(dA[0][0]), &n, &(iPIV[0]), &(dWork[0]), &nWork, &nInfo);
#endif

	return nInfo;
}

///////////////////////////////////////////////////////////////////////////////

int LAPACK::DTRTRI(
	char chUpperLower,
	char chDiagonal,
	DataMatrix<double> & dA
) {
	if (dA.GetRows() != dA.GetColumns()) {
		_EXCEPTIONT("Matrix A must be square.");
	}

	int n = dA.GetRows();

	int nInfo;

#ifdef USEACML
	dtrtri(chUpperLower, chDiagonal, n, &(dA[0][0]), n, &nInfo);
#endif
#ifdef USEESSL
	dtrtri(&chUpperLower, &chDiagonal, n, &(dA[0][0]), n, nInfo);
#endif
#if defined USEVECLIB || defined USEMKL
	dtrtri_(&chUpperLower, &chDiagonal, &n, &(dA[0][0]), &n, &nInfo);
#endif

	return nInfo;
}

///////////////////////////////////////////////////////////////////////////////

int LAPACK::DGEQRF(
	DataMatrix<double> & dA,
	DataVector<double> & dTau,
	DataVector<double> & dWork
) {
	int nRows = dA.GetColumns();
	int nCols = dA.GetRows();
	int nLDA = nRows;

	int nWork = dWork.GetRows();

	int nInfo;

#if defined USEACML || defined USEESSL
	//dgeqrf(nRows, nCols, &(dA[0][0]), nLDA, &(dTau[0]),
	//	&(dWork[0]), nWork, &nInfo);
	nInfo = 0;
	_EXCEPTION();
#endif
#if defined USEVECLIB || defined USEMKL
	dgeqrf_(&nRows, &nCols, &(dA[0][0]), &nLDA, &(dTau[0]),
		&(dWork[0]), &nWork, &nInfo);
#endif

	return nInfo;
}

///////////////////////////////////////////////////////////////////////////////

int LAPACK::DORGQR(
	DataMatrix<double> & dA,
	DataVector<double> & dTau,
	DataVector<double> & dWork
) {
	int nRows = dA.GetColumns();
	int nCols = dA.GetRows();

	int nWork = dWork.GetRows();

	int nInfo;

	if (nWork < nCols) {
		_EXCEPTION2("Work matrix too small:  Expected %i, found %i.",
			nCols, nWork);
	}

#if defined USEACML || defined USEESSL
	//dorgqr(nRows, nCols, nCols, &(dA[0][0]), nRows, &(dTau[0]),
	//	&(dWork[0]), nWork, &nInfo);
	nInfo = 0;
	_EXCEPTION();
#endif
#if defined USEVECLIB || defined USEMKL
	dorgqr_(&nRows, &nCols, &nCols, &(dA[0][0]), &nRows, &(dTau[0]),
		&(dWork[0]), &nWork, &nInfo);
#endif

	return nInfo;

}

///////////////////////////////////////////////////////////////////////////////

void LAPACK::GeneralizedInverse(
	DataMatrix<double> & dA,
	DataMatrix<double> & dOut
) {
	unsigned int i;
	unsigned int j;
	unsigned int k;

	int nInfo;

	// Check that A is row-dominant
	if (dA.GetRows() > dA.GetColumns()) {
		_EXCEPTION2("Number of columns must be at least equal to the number of "
			"rows to compute a generalized inverse (in fortran ordering).  "
			"Found: (%i,%i)", dA.GetRows(), dA.GetColumns());
	}

	// Initialize matrices
	DataVector<double> dTau;
	dTau.Initialize(dA.GetRows());

	DataVector<double> dWork;
	dWork.Initialize(2 * dA.GetColumns());

	DataVector<int> iPIV;
	iPIV.Initialize(dA.GetRows());

	DataMatrix<double> dR;
	dR.Initialize(dA.GetRows(), dA.GetRows());
	dR.Zero();

	// QR decomposition
	nInfo = DGEQRF(dA, dTau, dWork);
	if (nInfo < 0) {
		_EXCEPTION1("Illegal value in LAPACK_DGEQRF: %i", nInfo);
	}

	// Copy matrix R
	for (i = 0; i < dA.GetRows(); i++) {
	for (j = i; j < dA.GetRows(); j++) {
		dR[j][i] = dA[j][i];
	}
	}

	// Calculate the inverse of R
	nInfo = DTRTRI('U', 'N', dR);
	if (nInfo < 0) {
		_EXCEPTION1("Illegal value in LAPACK_DTRTRI: %i", nInfo);
	}
	if (nInfo > 0) {
		_EXCEPTION1("Matrix R is singular to working precision: %i", nInfo);
	}

	// Calculate matrix Q
	nInfo = DORGQR(dA, dTau, dWork);
	if (nInfo < 0) {
		_EXCEPTION1("Illegal value in LAPACK_DORGQR: %i", nInfo);
	}

	// Initialize the output matrix
	dOut.Initialize(dA.GetColumns(), dA.GetRows());
	dOut.Zero();

	// Multiply
	for (i = 0; i < dA.GetColumns(); i++) {
	for (j = 0; j < dA.GetRows(); j++) {
	for (k = 0; k < dA.GetRows(); k++) {
		dOut[i][j] += dR[k][j] * dA[k][i];
	}
	}
	}
}

///////////////////////////////////////////////////////////////////////////////

