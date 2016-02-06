///////////////////////////////////////////////////////////////////////////////
///
///	\file    LinearAlgebra.h
///	\author  Paul Ullrich
///	\version July 26, 2010
///
///	<summary>
///		This file provides wrappers for performing linear algebra operations
///		on DataArray2D and DataArray1D objects.
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

#include "DataArray1D.h"
#include "DataArray2D.h"
#include "MathHelper.h"

#include <cstdlib>
#include <iostream>

///////////////////////////////////////////////////////////////////////////////

#ifdef TEMPEST_LAPACK_ACML_INTERFACE
#include <acml.h>
#endif

#ifdef TEMPEST_LAPACK_ESSL_INTERFACE
#include <essl.h>
#endif

#ifdef TEMPEST_LAPACK_FORTRAN_INTERFACE

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

/// General matrix matrix multiply
int dgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);

///	General matrix solver from CLAPACK
int dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);

/// Triangular matrix solver from CLAPACK
int dtpsv_(char *uplo, char *trace, char *diag, int *n, double *a, double *x, int *inc);

///	General banded matrix solver from CLAPACK
int dgbsv_(int *n, int *kl, int *ku, int *nrhs, double *ab, int *ldab, int *ipiv, double *b, int *ldb, int *info); 

/// LU decomposition from CLAPACK
int dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);

/// General banded matrix LU decomposition from CLAPACK
int dgbtrf_(int *m, int *n, int *kl, int *ku, double *ab, int *ldab, int *ipiv, int *info);

///	Solve a matrix system using LU decomposition from dgetrf
int dgetrs_(char *trans, int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);

///	Solve a general banded matrix system using LU decomposition from dgbtrf
int dgbtrs_(char *trans, int *n, int *kl, int *ku, int *nrhs, double *ab, int *ldab, int *ipiv, double *b, int *ldb, int *info);

///	General matrix inverse.
int dgetri_(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);

///	Triangular matrix inverse
int dtrtri_(char *uplo, char *diag, int *n, double *a, int *lda, int *info);

///	QR decomposition
int dgeqrf_(int *m, int *n, double *a, int *lda, double *tau, double *work, int *lwork, int *info);

/// Singular value decomposition
int dgesvd_(char *jobu, char *jobvt, int *m, int *n, double *a, int *lda, double *s, double *u, int *ldu, double *vt, int *ldvt, double *work, int *lwork, int *info);

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
		DataArray1D<double> & dX,
		DataArray1D<double> & dY,
		int nIncX = 1,
		int nIncY = 1
	) {
		int nRowsX = dX.GetRows();
		int nRowsY = dY.GetRows();

		if (nRowsX != nRowsY) {
			_EXCEPTIONT("Incompatible vectors.");
		}

#if defined TEMPEST_LAPACK_ACML_INTERFACE \
 || defined TEMPEST_LAPACK_ESSL_INTERFACE
		dswap(nRowsX, &(dX[0]), nIncX, &(dY[0]), nIncY);
#endif
#ifdef TEMPEST_LAPACK_FORTRAN_INTERFACE
		dswap_(&nRowsX, &(dX[0]), &nIncX, &(dY[0]), &nIncY);
#endif
	}

	///	<summary>
	///		Scale a vector by a constant.
	///	</summary>
	inline static void DSCAL(
		DataArray1D<double> & dX,
		double dAlpha,
		int nInc = 1
	) {
		int nRows = dX.GetRows();

#if defined TEMPEST_LAPACK_ACML_INTERFACE \
 || defined TEMPEST_LAPACK_ESSL_INTERFACE
		dscal(nRows, dAlpha, &(dX[0]), nInc);
#endif
#ifdef TEMPEST_LAPACK_FORTRAN_INTERFACE
		dscal_(&nRows, &dAlpha, &(dX[0]), &nInc);
#endif
	}

	///	<summary>
	///		Copy a vector X on top of a vector Y.
	///	</summary>
	inline static void DCOPY(
		DataArray1D<double> & dX,
		DataArray1D<double> & dY,
		int nIncX = 1,
		int nIncY = 1
	) {
		int nRowsX = dX.GetRows();
		int nRowsY = dY.GetRows();

		if (nRowsX != nRowsY) {
			dY.Allocate(nRowsX);
		}

#if defined TEMPEST_LAPACK_ACML_INTERFACE \
 || defined TEMPEST_LAPACK_ESSL_INTERFACE
		dcopy(nRowsX, &(dX[0]), nIncX, &(dY[0]), nIncY);
#endif
#ifdef TEMPEST_LAPACK_FORTRAN_INTERFACE
		dcopy_(&nRowsX, &(dX[0]), &nIncX, &(dY[0]), &nIncY);
#endif
	}

	///	<summary>
	///		Copy a vector into an array of doubles.
	///	</summary>
	inline static void DCOPY_A(
		DataArray1D<double> & dX,
		double * dY,
		int nIncX = 1,
		int nIncY = 1
	) {
		int nRowsX = dX.GetRows();

#if defined TEMPEST_LAPACK_ACML_INTERFACE \
 || defined TEMPEST_LAPACK_ESSL_INTERFACE
		dcopy(nRowsX, &(dX[0]), nIncX, dY, nIncY);
#endif
#ifdef TEMPEST_LAPACK_FORTRAN_INTERFACE
		dcopy_(&nRowsX, &(dX[0]), &nIncX, dY, &nIncY);
#endif
	}

	///	<summary>
	///		Compute alpha * X + Y, where alpha is a constant and X and Y are
	///		vectors.  Store the result in Y.
	///	</summary>
	inline static void DAXPY(
		double dAlpha,
		DataArray1D<double> & dX,
		DataArray1D<double> & dY,
		int nIncX = 1,
		int nIncY = 1
	) {
		int nRowsX = dX.GetRows();

#if defined TEMPEST_LAPACK_ACML_INTERFACE \
 || defined TEMPEST_LAPACK_ESSL_INTERFACE
		daxpy(nRowsX, dAlpha, &(dX[0]), nIncX, &(dY[0]), nIncY);
#endif
#ifdef TEMPEST_LAPACK_FORTRAN_INTERFACE
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
		DataArray1D<double> & dY,
		int nIncX = 1,
		int nIncY = 1
	) {
		int nRowsY = dY.GetRows();

#if defined TEMPEST_LAPACK_ACML_INTERFACE \
 || defined TEMPEST_LAPACK_ESSL_INTERFACE
		daxpy(nRowsY, dAlpha, dX, nIncX, &(dY[0]), nIncY);
#endif
#ifdef TEMPEST_LAPACK_FORTRAN_INTERFACE
		daxpy_(&nRowsY, &dAlpha, dX, &nIncX, &(dY[0]), &nIncY);
#endif
	}

	///	<summary>
	///		Dot product of two vectors.
	///	</summary>
	inline static double DDOT(
		DataArray1D<double> & dX,
		DataArray1D<double> & dY,
		int nIncX = 1,
		int nIncY = 1
	) {
		int nRowsX = dX.GetRows();
		int nRowsY = dY.GetRows();

		if (nRowsX != nRowsY) {
			_EXCEPTIONT("Incompatible vectors.");
		}

#if defined TEMPEST_LAPACK_ACML_INTERFACE \
 || defined TEMPEST_LAPACK_ESSL_INTERFACE
		return ddot(nRowsX, &(dX[0]), nIncX, &(dY[0]), nIncY);
#endif
#ifdef TEMPEST_LAPACK_FORTRAN_INTERFACE
		return ddot_(&nRowsX, &(dX[0]), &nIncX, &(dY[0]), &nIncY);
#endif
	}

	///	<summary>
	///		Dot product of two vectors.
	///	</summary>
	inline static double DDOT_A(
		DataArray1D<double> & dX,
		double * dY
	) {
		int nRowsX = dX.GetRows();
		int nIncX = 1;
		int nIncY = 1;

#if defined TEMPEST_LAPACK_ACML_INTERFACE \
 || defined TEMPEST_LAPACK_ESSL_INTERFACE
		return ddot(nRowsX, &(dX[0]), nIncX, dY, nIncY);
#endif
#ifdef TEMPEST_LAPACK_FORTRAN_INTERFACE
		return ddot_(&nRowsX, &(dX[0]), &nIncX, dY, &nIncY);
#endif
	}

	///	<summary>
	///		Compute the norm of a vector.
	///	</summary>
	inline static double DNORM2(
		DataArray1D<double> & dX
	) {
		int nRows = static_cast<int>(dX.GetRows());
		int nIncX = 1;

#if defined TEMPEST_LAPACK_ACML_INTERFACE \
 || defined TEMPEST_LAPACK_ESSL_INTERFACE
		return dnrm2(nRows, &(dX[0]), nIncX);
#endif
#ifdef TEMPEST_LAPACK_FORTRAN_INTERFACE
		return dnrm2_(&nRows, &(dX[0]), &nIncX);
#endif
	}

	///	<summary>
	///		Perform a multiply-add on a matrix.
	///	</summary>
	static void DGEMM(
		DataArray2D<double> & dC,
		DataArray2D<double> & dA,
		DataArray2D<double> & dB,
		double dAlpha,
		double dBeta
	);

///////////////////////////////////////////////////////////////////////////////
// LEVEL 2 BLAS ROUTINES
///////////////////////////////////////////////////////////////////////////////

public:
	///	<summary>
	///		Solve a linear system of the form A * X = B.
	///	</summary>
	static int DGESV(
		DataArray2D<double> & dA,
		DataArray1D<double> & dBX,
		DataArray1D<int> & iPIV
	);

	///	<summary>
	///		Solve a linear system of the form A * X = B, with multiple RHS
	///		vectors.
	///	</summary>
	static int DGESV(
		DataArray2D<double> & dA,
		DataArray2D<double> & dBX,
		DataArray1D<int> & iPIV
	);

	///	<summary>
	///		Solve a linear system of the form A * X = B.
	///	</summary>
	static int DGBSV(
		DataArray2D<double> & dA,
		DataArray1D<double> & dBX,
		DataArray1D<int> & iPIV,
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
		DataArray1D<double> & dA,
		DataArray1D<double> & dX
	);

	///	<summary>
	///		Calculate the LU decomposition of a given general matrix.
	///	</summary>
	static int DGETRF(
		DataArray2D<double> & dA,
		DataArray1D<int> & iPIV
	);

	///	<summary>
	///		Calculate the LU decomposition of a given banded matrix.
	///	</summary>
	static int DGBTRF(
		DataArray2D<double> & dA,
		DataArray1D<int> & iPIV,
		int iKL,
		int iKU
	);

	///	<summary>
	///		Solve the matrix system using LU decomposition from DGETRF.
	///	</summary>
	static int DGETRS(
		char chTrans,
		DataArray2D<double> & dA,
		DataArray1D<double> & dB,
		DataArray1D<int> & iPIV
	);

	///	<summary>
	///		Solve the banded matrix system using LU decomposition from DGBTRF.
	///	</summary>
	static int DGBTRS(
		char chTrans,
		DataArray2D<double> & dA,
		DataArray1D<double> & dB,
		DataArray1D<int> & iPIV,
		int iKL,
		int iKU
	);

	///	<summary>
	///		Calculate the inverse of a given general matrix.
	///	</summary>
	static int DGETRI(
		DataArray2D<double> & dA,
		DataArray1D<int> & iPIV,
		DataArray1D<double> & dWork
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
		DataArray2D<double> & dA
	);

	///	<summary>
	///		Perform a QR decomposition on the given matrix.
	///	</summary>
	static int DGEQRF(
		DataArray2D<double> & dA,
		DataArray1D<double> & dTau,
		DataArray1D<double> & dWork
	);

	///	<summary>
	///		Reconstruct the matrix Q from a QR decomposition.
	///	</summary>
	static int DORGQR(
		DataArray2D<double> & dA,
		DataArray1D<double> & dTau,
		DataArray1D<double> & dWork
	);

	///	<summary>
	///		Calculate the generalized inverse of a given matrix A.
	///	</summary>
	///	<parameters>
	///		dA   - Input matrix A in Fortran order.
	///		dOut - Generalized inverse of A in Fortran order.
	///	</parameters>
	static void GeneralizedInverse(
		DataArray2D<double> & dA,
		DataArray2D<double> & dOut
	);

	///	<summary>
	///		Calculate the generalized inverse of a given matrix A
	///		using singular value decomposition (SVD).
	///	</summary>
	///	<parameters>
	///		dA   - Input matrix A in Fortran order.
	///		dOut - Generalized inverse of A in Fortran order.
	///	</parameters>
	static void GeneralizedInverseSVD(
		DataArray2D<double> & dA,
		DataArray2D<double> & dOut
	);

};

///////////////////////////////////////////////////////////////////////////////

#endif

