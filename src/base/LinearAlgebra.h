///////////////////////////////////////////////////////////////////////////////
///
///	\file    LinearAlgebra.h
///	\author  Paul Ullrich
///	\version July 26, 2010
///
///	<summary>
///		This file provides wrappers for performing linear algebra operations
///		on DataMatrix and DataVector objects.
///	</summary>
///	<remarks>
///		Copyright 2000-2010 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _LINEARALGEBRA_H_
#define _LINEARALGEBRA_H_

///////////////////////////////////////////////////////////////////////////////

#include "Exception.h"

#include "DataVector.h"
#include "DataMatrix.h"
#include "MathHelper.h"

#include <cstdlib>
#include <iostream>

///////////////////////////////////////////////////////////////////////////////

#ifdef USEACML
#include <acml.h>
#endif

#ifdef USEESSL
#include <essl.h>
#endif

#if defined USEVECLIB || defined USEMKL

extern "C" {

///	Swap two vectors
void dswap_(int *n, double *x, int *incx, double *y, int *incy);

///	Scale a vector
void dscal_(int *n, double *alpha, double *x, int *incx);

///	Copy a vector
void dcopy_(int *n, double *x, int *incx, double *y, int *incy);

/// Compute a * X + Y
void daxpy_(int *n, double *alpha, double *x, int *incx, double *y, int *incy);

/// Dot product of two vectors
double ddot_(int *n, double *x, int *incx, double *y, int *incy);

///	Compute the norm of a vector
double dnrm2_(int *n, double *x, int *incx);

///	General matrix solver from CLAPACK
int dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);

/// Triangular matrix solver from CLAPACK
int dtpsv_(char *uplo, char *trace, char *diag, int *n, double *a, double *x, int *inc);

///	General banded matrix solver from CLAPACK
int dgbsv_(int *n, int *kl, int *ku, int *nrhs, double *ab, int *ldab, int *ipiv, double *b, int *ldb, int *info); 

/// LU decomposition from CLAPACK from CLAPACK
int dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);

///	General matrix inverse.
int dgetri_(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);

///	Triangular matrix inverse
int dtrtri_(char *uplo, char *diag, int *n, double *a, int *lda, int *info);

///	QR decomposition
int dgeqrf_(int *m, int *n, double *a, int *lda, double *tau, double *work, int *lwork, int *info);

///	Reconstruction of matrix Q in QR decomposition
int dorgqr_(int *m, int *n, int *k, double * a, int *lda, double *tau, double *work, int *lwork, int *info);

}

#endif

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		This class wraps the LAPACK linear algebra functions.
///	</summary>
class LAPACK {

///////////////////////////////////////////////////////////////////////////////
// LEVEL 1 BLAS ROUTINES
///////////////////////////////////////////////////////////////////////////////

public:
	///	<summary>
	///		Swap two vectors.
	///	</summary>
	inline static void DSWAP(
		DataVector<double> & dX,
		DataVector<double> & dY,
		int nIncX = 1,
		int nIncY = 1
	) {
		int nRowsX = dX.GetRows();
		int nRowsY = dY.GetRows();

		if (nRowsX != nRowsY) {
			_EXCEPTIONT("Incompatible vectors.");
		}

#if defined USEACML || defined USEESSL
		dswap(nRowsX, &(dX[0]), nIncX, &(dY[0]), nIncY);
#endif
#if defined USEVECLIB || defined USEMKL
		dswap_(&nRowsX, &(dX[0]), &nIncX, &(dY[0]), &nIncY);
#endif
	}

	///	<summary>
	///		Scale a vector by a constant.
	///	</summary>
	inline static void DSCAL(
		DataVector<double> & dX,
		double dAlpha,
		int nInc = 1
	) {
		int nRows = dX.GetRows();

#if defined USEACML || defined USEESSL
		dscal(nRows, dAlpha, &(dX[0]), nInc);
#endif
#if defined USEVECLIB || defined USEMKL
		dscal_(&nRows, &dAlpha, &(dX[0]), &nInc);
#endif
	}

	///	<summary>
	///		Copy a vector X on top of a vector Y.
	///	</summary>
	inline static void DCOPY(
		DataVector<double> & dX,
		DataVector<double> & dY,
		int nIncX = 1,
		int nIncY = 1
	) {
		int nRowsX = dX.GetRows();
		int nRowsY = dY.GetRows();

		if (nRowsX != nRowsY) {
			dY.Initialize(nRowsX);
		}

#if defined USEACML || defined USEESSL
		dcopy(nRowsX, &(dX[0]), nIncX, &(dY[0]), nIncY);
#endif
#if defined USEVECLIB || defined USEMKL
		dcopy_(&nRowsX, &(dX[0]), &nIncX, &(dY[0]), &nIncY);
#endif
	}

	///	<summary>
	///		Copy a vector into an array of doubles.
	///	</summary>
	inline static void DCOPY_A(
		DataVector<double> & dX,
		double * dY,
		int nIncX = 1,
		int nIncY = 1
	) {
		int nRowsX = dX.GetRows();

#if defined USEACML || defined USEESSL
		dcopy(nRowsX, &(dX[0]), nIncX, dY, nIncY);
#endif
#if defined USEVECLIB || defined USEMKL
		dcopy_(&nRowsX, &(dX[0]), &nIncX, dY, &nIncY);
#endif
	}

	///	<summary>
	///		Compute alpha * X + Y, where alpha is a constant and X and Y are
	///		vectors.  Store the result in Y.
	///	</summary>
	inline static void DAXPY(
		double dAlpha,
		DataVector<double> & dX,
		DataVector<double> & dY,
		int nIncX = 1,
		int nIncY = 1
	) {
		int nRowsX = dX.GetRows();

#if defined USEACML || defined USEESSL
		daxpy(nRowsX, dAlpha, &(dX[0]), nIncX, &(dY[0]), nIncY);
#endif
#if defined USEVECLIB || defined USEMKL
		daxpy_(&nRowsX, &dAlpha, &(dX[0]), &nIncX, &(dY[0]), &nIncY);
#endif
	}

	///	<summary>
	///		Compute alpha * X + Y, where alpha is a constant and X and Y are
	///		vectors.  Store the result in Y.
	///	</summary>
	inline static void DAXPY_A(
		double dAlpha,
		double * dX,
		DataVector<double> & dY,
		int nIncX = 1,
		int nIncY = 1
	) {
		int nRowsY = dY.GetRows();

#if defined USEACML || defined USEESSL
		daxpy(nRowsY, dAlpha, dX, nIncX, &(dY[0]), nIncY);
#endif
#if defined USEVECLIB || defined USEMKL
		daxpy_(&nRowsY, &dAlpha, dX, &nIncX, &(dY[0]), &nIncY);
#endif
	}

	///	<summary>
	///		Dot product of two vectors.
	///	</summary>
	inline static double DDOT(
		DataVector<double> & dX,
		DataVector<double> & dY,
		int nIncX = 1,
		int nIncY = 1
	) {
		int nRowsX = dX.GetRows();
		int nRowsY = dY.GetRows();

		if (nRowsX != nRowsY) {
			_EXCEPTIONT("Incompatible vectors.");
		}

#if defined USEACML || defined USEESSL
		return ddot(nRowsX, &(dX[0]), nIncX, &(dY[0]), nIncY);
#endif
#if defined USEVECLIB || defined USEMKL
		return ddot_(&nRowsX, &(dX[0]), &nIncX, &(dY[0]), &nIncY);
#endif
	}

	///	<summary>
	///		Dot product of two vectors.
	///	</summary>
	inline static double DDOT_A(
		DataVector<double> & dX,
		double * dY
	) {
		int nRowsX = dX.GetRows();
		int nIncX = 1;
		int nIncY = 1;

#if defined USEACML || defined USEESSL
		return ddot(nRowsX, &(dX[0]), nIncX, dY, nIncY);
#endif
#if defined USEVECLIB || defined USEMKL
		return ddot_(&nRowsX, &(dX[0]), &nIncX, dY, &nIncY);
#endif
	}

	///	<summary>
	///		Compute the norm of a vector.
	///	</summary>
	inline static double DNORM2(
		DataVector<double> & dX
	) {
		int nRows = static_cast<int>(dX.GetRows());
		int nIncX = 1;

#if defined USEACML || defined USEESSL
		return dnrm2(nRows, &(dX[0]), nIncX);
#endif
#if defined USEVECLIB || defined USEMKL
		return dnrm2_(&nRows, &(dX[0]), &nIncX);
#endif
	}

///////////////////////////////////////////////////////////////////////////////
// LEVEL 2 BLAS ROUTINES
///////////////////////////////////////////////////////////////////////////////

public:
	///	<summary>
	///		Solve a linear system of the form A * X = B.
	///	</summary>
	static int DGESV(
		DataMatrix<double> & dA,
		DataVector<double> & dBX,
		DataVector<int> & iPIV
	);

	///	<summary>
	///		Solve a linear system of the form A * X = B, with multiple RHS
	///		vectors.
	///	</summary>
	static int DGESV(
		DataMatrix<double> & dA,
		DataMatrix<double> & dBX,
		DataVector<int> & iPIV
	);

	///	<summary>
	///		Solve a linear system of the form A * X = B.
	///	</summary>
	static int DGBSV(
		DataMatrix<double> & dA,
		DataVector<double> & dBX,
		DataVector<int> & iPIV,
		int nKL,
		int nKU
	);

	///	<summary>
	///		Solve a linear system of the form A * X = B, where A is
	///		a triangular matrix.
	///	</summary>
	///	<parameters>
	///		chUpperLower - 'U' if A is upper triangular
	///		               'L' if A is lower triangular
	///		chTrans      - 'N' to solve for A * X = B
	///		               'T' to solve for (A^T) * X = B
	///		chDiagonal   - 'N' if A is non-unit triangular
	///		               'U' if A is unit triangular
	///		nN           - Order of the array dA
	///	</parameters>
	static int DTPSV(
		char chUpperLower,
		char chTrans,
		char chDiagonal,
		int nN,
		DataVector<double> & dA,
		DataVector<double> & dX
	);

	///	<summary>
	///		Calculate the LU decomposition of a given general matrix.
	///	</summary>
	static int DGETRF(
		DataMatrix<double> & dA,
		DataVector<int> & iPIV
	);

	///	<summary>
	///		Calculate the inverse of a given general matrix.
	///	</summary>
	static int DGETRI(
		DataMatrix<double> & dA,
		DataVector<int> & iPIV,
		DataVector<double> & dWork
	);

	///	<summary>
	///		Calculate the inverse of a given triangular matrix.
	///	</summary>
	///	<parameters>
	///		chUpperLower - 'U' if A is upper triangular
	///		               'L' if A is lower triangular
	///		chDiagonal   - 'N' if A is non-unit triangular
	///		               'U' if A is unit triangular
	///	</parameters>
	static int DTRTRI(
		char chUpperLower,
		char chDiagonal,
		DataMatrix<double> & dA
	);

	///	<summary>
	///		Perform a QR decomposition on the given matrix.
	///	</summary>
	static int DGEQRF(
		DataMatrix<double> & dA,
		DataVector<double> & dTau,
		DataVector<double> & dWork
	);

	///	<summary>
	///		Reconstruct the matrix Q from a QR decomposition.
	///	</summary>
	static int DORGQR(
		DataMatrix<double> & dA,
		DataVector<double> & dTau,
		DataVector<double> & dWork
	);

	///	<summary>
	///		Calculate the generalized inverse of a given matrix A.
	///	</summary>
	///	<parameters>
	///		dA   - Input matrix A in Fortran order.
	///		dOut - Generalized inverse of A in Fortran order.
	///	</parameters>
	static void GeneralizedInverse(
		DataMatrix<double> & dA,
		DataMatrix<double> & dOut
	);
};

///////////////////////////////////////////////////////////////////////////////

#endif

