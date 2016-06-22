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
#include <cmath>

//////////////////////////////////////////////////////////////////////////////

void LAPACK::DGEMM(
	DataArray2D<double> & dC,
	DataArray2D<double> & dA,
	DataArray2D<double> & dB,
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

#ifdef TEMPEST_LAPACK_ACML_INTERFACE
	// Call the matrix solve
	dgemm(cTransA, cTransB, m, n, k,
		dAlpha, &(dA[0][0]), nLDA, &(dB[0][0]), nLDB,
		dBeta, &(dC[0][0]), nLDC);
#endif
#ifdef TEMPEST_LAPACK_ESSL_INTERFACE
	// Call the matrix solve
	dgemm(cTransA, cTransB, m, n, k,
		dAlpha, &(dA[0][0]), nLDA, &(dB[0][0]), nLDB,
		dBeta, &(dC[0][0]), nLDC);
#endif
#ifdef TEMPEST_LAPACK_FORTRAN_INTERFACE
	// Call the matrix solve
	dgemm_(&cTransA, &cTransB, &m, &n, &k, &dAlpha,
		&(dA[0][0]), &nLDA, &(dB[0][0]), &nLDB, &dBeta, &(dC[0][0]), &nLDC);
#endif
}

//////////////////////////////////////////////////////////////////////////////

int LAPACK::DGESV(
	DataArray2D<double> & dA,
	DataArray1D<double> & dBX,
	DataArray1D<int> & iPIV
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

#ifdef TEMPEST_LAPACK_ACML_INTERFACE
	// Call the matrix solve
	dgesv(n, nRHS, &(dA[0][0]), nLDA, &(iPIV[0]), &(dBX[0]), nLDB, &nInfo);
#endif
#ifdef TEMPEST_LAPACK_ESSL_INTERFACE
	// Call the matrix solve
	dgesv(n, nRHS, &(dA[0][0]), nLDA, &(iPIV[0]), &(dBX[0]), nLDB, nInfo);
#endif
#ifdef TEMPEST_LAPACK_FORTRAN_INTERFACE
	// Call the matrix solve
	dgesv_(
		&n, &nRHS,
		&(dA[0][0]), &nLDA, &(iPIV[0]), &(dBX[0]), &nLDB, &nInfo);
#endif

	return nInfo;
}

//////////////////////////////////////////////////////////////////////////////

int LAPACK::DGESV(
	DataArray2D<double> & dA,
	DataArray2D<double> & dBX,
	DataArray1D<int> & iPIV
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

#ifdef TEMPEST_LAPACK_ACML_INTERFACE
	// Call the matrix solve
	dgesv(n, nRHS,
		&(dA[0][0]), nLDA, &(iPIV[0]), &(dBX[0][0]), nLDB, &nInfo);
#endif
#ifdef TEMPEST_LAPACK_ESSL_INTERFACE
	// Call the matrix solve
	dgesv(n, nRHS,
		&(dA[0][0]), nLDA, &(iPIV[0]), &(dBX[0][0]), nLDB, nInfo);
#endif
#ifdef TEMPEST_LAPACK_FORTRAN_INTERFACE
	// Call the matrix solve
	dgesv_(
		&n, &nRHS,
		&(dA[0][0]), &nLDA, &(iPIV[0]), &(dBX[0][0]), &nLDB, &nInfo);
#endif

	return nInfo;
}

///////////////////////////////////////////////////////////////////////////////

int LAPACK::DGBSV(
	DataArray2D<double> & dA,
	DataArray1D<double> & dBX,
	DataArray1D<int> & iPIV,
	int nKL,
	int nKU
) {
	// Check dimensions
	if (dA.GetRows() < 2 * nKL + nKU + 1) {
		_EXCEPTIONT("Matrix A has insufficient rows for DGBSV");
	}
	if (dBX.GetRows() < dA.GetRows()) {
		_EXCEPTIONT("Matrix A / B dimension mismatch in DGBSV");
	}
	if (iPIV.GetRows() < dA.GetRows()) {
		_EXCEPTIONT("Matrix A / IPIV dimension mismatch in DGBSV");
	}

	// Store CLAPACK parameters
	int n     = dA.GetRows();
	int nRHS  = 1;
	int nLDAB = dA.GetColumns();
	int nLDB  = n;

	int nInfo;

#ifdef TEMPEST_LAPACK_ACML_INTERFACE
	// Call the banded diagonal matrix solve
	dgbsv(
		n, nKL, nKU, nRHS,
		&(dA[0][0]), nLDAB, &(iPIV[0]), &(dBX[0]), nLDB, &nInfo);
#endif
#ifdef TEMPEST_LAPACK_ESSL_INTERFACE
	// Call the banded diagonal matrix solve
	dgbf(&(dA[0][0]), nLDAB, n, nKL, nKU, &(iPIV[0]));
	dgbs(&(dA[0][0]), nLDAB, n, nKL, nKU, &(iPIV[0]), &(dBX[0]));
	nInfo = 0;
#endif
#ifdef TEMPEST_LAPACK_FORTRAN_INTERFACE
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
	DataArray1D<double> & dA,
	DataArray1D<double> & dX
) {
	int iIncX = 1;

#if defined TEMPEST_LAPACK_ACML_INTERFACE \
 || defined TEMPEST_LAPACK_ESSL_INTERFACE
	_EXCEPTION();
#endif
#ifdef TEMPEST_LAPACK_FORTRAN_INTERFACE
	// Call the triangular matrix solve
	dtpsv_(
		&chUpperLower, &chTrans, &chDiagonal,
		&nN, &(dA[0]), &(dX[0]), &iIncX);
#endif

	return 0;
}

///////////////////////////////////////////////////////////////////////////////

int LAPACK::DGETRF(
	DataArray2D<double> & dA,
	DataArray1D<int> & iPIV
) {
	int m = dA.GetRows();
	int n = dA.GetColumns();

	int lda = m;

	int nInfo;

#ifdef TEMPEST_LAPACK_ACML_INTERFACE
	dgetrf(m, n, &(dA[0][0]), lda, &(iPIV[0]), &nInfo);
#endif
#ifdef TEMPEST_LAPACK_ESSL_INTERFACE
	dgetrf(m, n, &(dA[0][0]), lda, &(iPIV[0]), nInfo);
#endif
#ifdef TEMPEST_LAPACK_FORTRAN_INTERFACE
	dgetrf_(&m, &n, &(dA[0][0]), &lda, &(iPIV[0]), &nInfo);
#endif

	return nInfo;
}

///////////////////////////////////////////////////////////////////////////////

int LAPACK::DGBTRF(
	DataArray2D<double> & dA,
	DataArray1D<int> & iPIV,
	int iKL,
	int iKU
) {
	int m = dA.GetRows();
	int n = dA.GetColumns();

	int lda = m;

	int nInfo;

#ifdef TEMPEST_LAPACK_ACML_INTERFACE
	dgbtrf(m, n, &iKL, &iKU, &(dA[0][0]), lda, &(iPIV[0]), &nInfo);
#endif
#ifdef TEMPEST_LAPACK_ESSL_INTERFACE
	dgbtrf(m, n, iKL, iKU, &(dA[0][0]), lda, &(iPIV[0]), nInfo);
#endif
#ifdef TEMPEST_LAPACK_FORTRAN_INTERFACE
	dgbtrf_(&m, &n, &iKL, &iKU, &(dA[0][0]), &lda, &(iPIV[0]), &nInfo);
#endif

	return nInfo;
}

///////////////////////////////////////////////////////////////////////////////

int LAPACK::DGETRS(
	char chTrans,
	DataArray2D<double> & dA,
	DataArray1D<double> & dB,
	DataArray1D<int> & iPIV
) {
	int m = dA.GetRows();
	int n = dA.GetColumns();

	int nrhs = 1;
	int lda = m;
	int ldb = m;

	int nInfo = 0;

#ifdef TEMPEST_LAPACK_ACML_INTERFACE
	dgetrs(chTrans, n, nrhs, &(dA[0][0]), lda, &(iPIV[0]), &(dB[0]), ldb, &nInfo);
#endif
#ifdef TEMPEST_LAPACK_ESSL_INTERFACE
	dgetrs(chTrans, n, nrhs, &(dA[0][0]), lda, &(iPIV[0]), &(dB[0]), ldb, nInfo);
#endif
#ifdef TEMPEST_LAPACK_FORTRAN_INTERFACE
	dgetrs_(&chTrans, &n, &nrhs, &(dA[0][0]), &lda, &(iPIV[0]), &(dB[0]), &ldb, &nInfo);
#endif

	return nInfo;
}

///////////////////////////////////////////////////////////////////////////////

int LAPACK::DGBTRS(
	char chTrans,
	DataArray2D<double> & dA,
	DataArray1D<double> & dB,
	DataArray1D<int> & iPIV,
	int iKL,
	int iKU
) {
	int m = dA.GetRows();
	int n = dA.GetColumns();

	int lda = m;
	int ldb = dB.GetRows();
	int nRHS = 1;

	int nInfo;

#ifdef TEMPEST_LAPACK_ACML_INTERFACE
	dgbtrs(chTrans, n, iKL, iKU, nRHS, &(dA[0][0]), lda, &(iPIV[0]), &(dB[0]), ldb, &nInfo);
#endif
#ifdef TEMPEST_LAPACK_ESSL_INTERFACE
	dgbtrs(chTrans, n, iKL, iKU, nRHS, &(dA[0][0]), lda, &(iPIV[0]), &(dB[0]), ldb, nInfo);
#endif
#ifdef TEMPEST_LAPACK_FORTRAN_INTERFACE
	dgbtrs_(&chTrans, &n, &iKL, &iKU, &nRHS, &(dA[0][0]), &lda, &(iPIV[0]), &(dB[0]), &ldb, &nInfo);
#endif

	return nInfo;
}

///////////////////////////////////////////////////////////////////////////////

int LAPACK::DGETRI(
	DataArray2D<double> & dA,
	DataArray1D<int> & iPIV,
	DataArray1D<double> & dWork
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

#if defined TEMPEST_LAPACK_ACML_INTERFACE \
 || defined TEMPEST_LAPACK_ESSL_INTERFACE
	//dgetri(n, &(dA[0][0]), n, &(iPIV[0]), &(dWork[0]), nWork, &nInfo);
	nInfo = 0;
	_EXCEPTIONT("Unimplemented.");
#endif
#ifdef TEMPEST_LAPACK_FORTRAN_INTERFACE
	dgetri_(&n, &(dA[0][0]), &n, &(iPIV[0]), &(dWork[0]), &nWork, &nInfo);
#endif

	return nInfo;
}

///////////////////////////////////////////////////////////////////////////////

int LAPACK::DTRTRI(
	char chUpperLower,
	char chDiagonal,
	DataArray2D<double> & dA
) {
	if (dA.GetRows() != dA.GetColumns()) {
		_EXCEPTIONT("Matrix A must be square.");
	}

	int n = dA.GetRows();

	int nInfo;

#ifdef TEMPEST_LAPACK_ACML_INTERFACE
	dtrtri(chUpperLower, chDiagonal, n, &(dA[0][0]), n, &nInfo);
#endif
#ifdef TEMPEST_LAPACK_ESSL_INTERFACE
	dtrtri(&chUpperLower, &chDiagonal, n, &(dA[0][0]), n, nInfo);
#endif
#ifdef TEMPEST_LAPACK_FORTRAN_INTERFACE
	dtrtri_(&chUpperLower, &chDiagonal, &n, &(dA[0][0]), &n, &nInfo);
#endif

	return nInfo;
}

///////////////////////////////////////////////////////////////////////////////

int LAPACK::DGEQRF(
	DataArray2D<double> & dA,
	DataArray1D<double> & dTau,
	DataArray1D<double> & dWork
) {
	int nRows = dA.GetColumns();
	int nCols = dA.GetRows();
	int nLDA = nRows;

	int nWork = dWork.GetRows();

	int nInfo;

#if defined TEMPEST_LAPACK_ACML_INTERFACE \
 || defined TEMPEST_LAPACK_ESSL_INTERFACE
	//dgeqrf(nRows, nCols, &(dA[0][0]), nLDA, &(dTau[0]),
	//	&(dWork[0]), nWork, &nInfo);
	nInfo = 0;
	_EXCEPTION();
#endif
#ifdef TEMPEST_LAPACK_FORTRAN_INTERFACE
	dgeqrf_(&nRows, &nCols, &(dA[0][0]), &nLDA, &(dTau[0]),
		&(dWork[0]), &nWork, &nInfo);
#endif

	return nInfo;
}

///////////////////////////////////////////////////////////////////////////////

int LAPACK::DORGQR(
	DataArray2D<double> & dA,
	DataArray1D<double> & dTau,
	DataArray1D<double> & dWork
) {
	int nRows = dA.GetColumns();
	int nCols = dA.GetRows();

	int nWork = dWork.GetRows();

	int nInfo;

	if (nWork < nCols) {
		_EXCEPTION2("Work matrix too small:  Expected %i, found %i.",
			nCols, nWork);
	}

#if defined TEMPEST_LAPACK_ACML_INTERFACE \
 || defined TEMPEST_LAPACK_ESSL_INTERFACE
	//dorgqr(nRows, nCols, nCols, &(dA[0][0]), nRows, &(dTau[0]),
	//	&(dWork[0]), nWork, &nInfo);
	nInfo = 0;
	_EXCEPTION();
#endif
#ifdef TEMPEST_LAPACK_FORTRAN_INTERFACE
	dorgqr_(&nRows, &nCols, &nCols, &(dA[0][0]), &nRows, &(dTau[0]),
		&(dWork[0]), &nWork, &nInfo);
#endif

	return nInfo;

}

///////////////////////////////////////////////////////////////////////////////

void LAPACK::GeneralizedInverse(
	DataArray2D<double> & dA,
	DataArray2D<double> & dOut
) {
	int nInfo;

	// Check that A is row-dominant
	if (dA.GetRows() > dA.GetColumns()) {
		_EXCEPTION2("Number of columns must be at least equal to the number of "
			"rows to compute a generalized inverse (in fortran ordering).  "
			"Found: (%i,%i)", dA.GetRows(), dA.GetColumns());
	}

	// Initialize matrices
	DataArray1D<double> dTau(dA.GetRows());

	DataArray1D<double> dWork(2 * dA.GetColumns());

	DataArray1D<int> iPIV(dA.GetRows());

	DataArray2D<double> dR(dA.GetRows(), dA.GetRows());

	// QR decomposition
	nInfo = DGEQRF(dA, dTau, dWork);
	if (nInfo < 0) {
		_EXCEPTION1("Illegal value in LAPACK_DGEQRF: %i", nInfo);
	}

	// Copy matrix R
	for (size_t i = 0; i < dA.GetRows(); i++) {
	for (size_t j = i; j < dA.GetRows(); j++) {
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
	dOut.Allocate(dA.GetColumns(), dA.GetRows());

	// Multiply
	for (size_t i = 0; i < dA.GetColumns(); i++) {
	for (size_t j = 0; j < dA.GetRows(); j++) {
	for (size_t k = 0; k < dA.GetRows(); k++) {
		dOut[i][j] += dR[k][j] * dA[k][i];
	}
	}
	}
}

///////////////////////////////////////////////////////////////////////////////

void LAPACK::GeneralizedInverseSVD(
	DataArray2D<double> & dA,
	DataArray2D<double> & dOut
) {
	dOut.Allocate(dA.GetRows(), dA.GetColumns());

	// Compute SVD
	DataArray2D<double> dU(dA.GetColumns(), dA.GetColumns());
	DataArray2D<double> dVT(dA.GetRows(), dA.GetRows());

	int dimbig;
	int dimsmall;

	if (dA.GetRows() > dA.GetColumns()) {
		dimbig = dA.GetRows();
		dimsmall = dA.GetColumns();
	} else {
		dimbig = dA.GetColumns();
		dimsmall = dA.GetRows();
	}

	int lwork = 5 * dimbig;

	DataArray1D<double> dS(dimsmall);

	DataArray1D<double> dWork(lwork);

	char jobu = 'S';
	char jobvt = 'S';

	int m = dA.GetColumns();
	int n = dA.GetRows();

	int nInfo = 0;

#ifdef TEMPEST_LAPACK_FORTRAN_INTERFACE
	dgesvd_(&jobu, &jobvt, &m, &n, &(dA[0][0]), &m, &(dS[0]), &(dU[0][0]), &m, &(dVT[0][0]), &n, &(dWork[0]), &lwork, &nInfo);
#else
	_EXCEPTIONT("Unimplemented");
#endif

	if (nInfo > 0) {
		_EXCEPTIONT("Convergence failure in LAPACK_DGESVD");
	}
	if (nInfo < 0) {
		_EXCEPTION1("Illegal value in LAPACK_DGESVD: %i", nInfo);
	}

	// Invert the singular values
	for (int i = 0; i < dimsmall; i++) {
		if (fabs(dS[i]) < 1.0e-14) {
			dS[i] = 0.0;
		} else {
			dS[i] = 1.0 / dS[i];
		}
	}

	// Calculate pseudo-inverse
	if (dA.GetColumns() > dA.GetRows()) {

		// Apply inverted singular values to U
		for (int i = 0; i < dA.GetColumns(); i++) {
		for (int j = 0; j < dA.GetRows(); j++) {
			dU[i][j] *= dS[j];
		}
		}

		// Calculate matrix product
		for (int i = 0; i < dA.GetRows(); i++) {
		for (int j = 0; j < dA.GetColumns(); j++) {
		for (int k = 0; k < dA.GetRows(); k++) {
			dOut[i][j] += dVT[i][k] * dU[j][k];
		}
		}
		}

	} else {
		// Apply inverted singular values to VT
		for (int i = 0; i < dA.GetRows(); i++) {
		for (int j = 0; j < dA.GetColumns(); j++) {
			dVT[i][j] *= dS[j];
		}
		}

		// Calculate matrix product
		for (int i = 0; i < dA.GetRows(); i++) {
		for (int j = 0; j < dA.GetColumns(); j++) {
		for (int k = 0; k < dA.GetColumns(); k++) {
			dOut[i][j] += dVT[i][k] * dU[k][j];
		}
		}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

