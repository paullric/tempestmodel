///////////////////////////////////////////////////////////////////////////////
///
///	\file    Polynomial.h
///	\author  Paul Ullrich
///	\version March 30, 2013
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

#ifndef _POLYNOMIAL_H_
#define _POLYNOMIAL_H_

#include "Exception.h"

///////////////////////////////////////////////////////////////////////////////

class LagrangeWeights {
public:
	///	<summary>
	///		Constructor.
	///	</summary>
	LagrangeWeights();

private:
	///	<summary>
	///		Number of weights.
	///	</summary>
	int m_nWeightCount;

	///	<summary>
	///		Array of weights.
	///	</summary>
	double * m_dWeights;
};

///////////////////////////////////////////////////////////////////////////////

class Polynomial {

public:
	///	<summary>
	///		Determine the weights of a Lagrangian polynomial through the
	///		specified points dX which is sampled at point dXsample.
	///	</summary>
	static void LagrangianPolynomialCoeffs(
		int nPoints,
		const double * dX,
		double * dWeights,
		double dXsample
	);

	///	<summary>
	///		Determine the weights of the first derivative of a Lagrangian
	///		polynomial through the specified points dX which is sampled at
	///		point dXsample.
	///	</summary>
	static void DiffLagrangianWeights(
		int nPoints,
		const double * dX,
		double * dWeights,
		double dXsample
	);

	///	<summary>
	///		Determine the weights of the second derivative of a Lagrangian
	///		polynomial through the specified points dX which is sampled at
	///		point dXsample.
	///	</summary>
	static void Diff2LagrangianCoeffs(
		int nPoints,
		const double * dX,
		double * dWeights,
		double dXsample
	);

	///	<summary>
	///		Determine the weights of the third derivative of a Lagrangian
	///		polynomial through the specified points dX which is sampled at
	///		point dXsample.
	///	</summary>
	static void Diff3LagrangianCoeffs(
		int nPoints,
		const double * dX,
		double * dWeights,
		double dXsample
	);

	///	<summary>
	///		Interpolate a polynomial through the given (X,Y) points and sample
	///		at point dXsample.  This method is faster and more computationally
	///		stable than InterpolateCoeffs.
	///	</summary>
	static double Interpolate(
		int nPoints,
		const double * dX,
		const double * dY,
		double dXsample
	);

	///	<summary>
	///		Obtain the coefficients a_i of a polynomial interpolated through the
	///		points (X,Y), where
	///		    p(x) = a_0 + a_1 (x - dXmid) + ... + a_(n-1) (x - dXmid)^(n-1)
	///	</summary>
	///	<parameter name = "dWorkingSpace">
	///		An external workspace of at least of size nPoints^2.
	///	</parameter>
	///	<parameter name = "iPivot">
	///		An external workspace of at least size nPoints.
	///	</parameter>
	static void InterpolateCoeffs(
		int nPoints,
		const double * dX,
		const double * dY,
		double * dA,
		double dXmid = 0.0,
		double * dWorkspace = NULL,
		int * iPivot = NULL
	);

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	Polynomial(int nCoeffCount);

	///	<summary>
	///		Copy constructor.
	///	</summary>
	Polynomial(const Polynomial & poly);

	///	<summary>
	///		Destructor.
	///	</summary>
	~Polynomial();

public:
	///	<summary>
	///		Assignment operator.
	///	</summary>
	Polynomial & operator= (const Polynomial & poly);

	///	<summary>
	///		Reserve space for the specified number of coefficients.
	///	</summary>
	void Reserve(int nCoeffCount);

public:
	///	<summary>
	///		Evaluate the polynomial at the specified point.
	///	</summary>
	double Evaluate(double dX) const;

public:
	///	<summary>
	///		Interpolate a polynomial through the given points.
	///	</summary>
	void Interpolate(
		const double * dX,
		const double * dY,
		double dXmid = 0.0
	);

	///	<summary>
	///		Differentiate the polynomial in place.
	///	</summary>
	void Differentiate();

	///	<summary>
	///		Integrate the polynomial in place.
	///	</summary>
	void Integrate();

private:
	///	<summary>
	///		Degree of the polynomial.
	///	</summary>
	int m_nCoeffCount;

	///	<summary>
	///		Array of polynomial coefficients.
	///	</summary>
	double * m_dCoeffs;

	///	<summary>
	///		Reference point of the polynomial.
	///	</summary>
	double m_dXMid;
};

///////////////////////////////////////////////////////////////////////////////

#endif

