///////////////////////////////////////////////////////////////////////////////
///
///	\file    JacobianFreeNewtonKrylov.h
///	\author  Paul Ullrich
///	\version May 20, 2013
///
///	<remarks>
///		Copyright 2000-2013 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _JACOBIANFREENEWTONKRYLOV_H_
#define _JACOBIANFREENEWTONKRYLOV_H_

///////////////////////////////////////////////////////////////////////////////

#include "DataVector.h"
#include "DataMatrix.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		An implementation of Jacobian-free Newton-Krylov type iteration for
///		the solution of nonlinear systems.
///	</summary>
class JacobianFreeNewtonKrylov {

public:
	///	<summary>
	///		Default constructor.
	///	</summary>
	JacobianFreeNewtonKrylov();

	///	<summary>
	///		Constructor.
	///	</summary>
	JacobianFreeNewtonKrylov(
		int nEquationCount,
		int nIterPerRestart,
		double dEpsilon = 1.0e-5
	);

public:
	///	<summary>
	///		Initializer.
	///	</summary>
	void InitializeJFNK(
		int nEquationCount,
		int nIterPerRestart,
		double dEpsilon = 1.0e-5
	);

private:
	///	<summary>
	///		Calculate the Givens rotation matrix with specified parameters.
	///	</summary>
	void GivensRotationMatrix(
		double dA,
		double dB,
		double & dC,
		double & dS
	) const;

public:
	///	<summary>
	///		Perform several Newton iterations until converged.
	///	</summary>
	bool PerformJFNK(
		double * dX,
		int nMaxIter,
		double dSolnAbsTolerance,
		int nGMRESMaxIter,
		double dGMRESTolerance
	);

	///	<summary>
	///		Perform one step of Newton iteration.
	///	</summary>
	///	<returns>
	///		The L2 norm of the difference between the initial state and
	///		updated state.
	///	</returns>
	double PerformJFNK_NewtonStep(
		double * dX,
		int nMaxIter,
		double dTolerance
	);

	///	<summary>
	///		Perform one step of Newton iteration (safer version).
	///	</summary>
	///	<returns>
	///		The L2 norm of the difference between the initial state and
	///		updated state.
	///	</returns>
	inline double PerformJFNK_NewtonStep_Safe(
		DataVector<double> & dX,
		int nMaxIter,
		double dTolerance
	) {
		// Check input size
		if (dX.GetRows() != m_dG.GetRows()) {
			_EXCEPTION2(
				"Input vector length mismatch: Expected length %i, received %i",
				m_dG.GetRows(), dX.GetRows());
		}
		
		return PerformJFNK_NewtonStep(dX, nMaxIter, dTolerance);
	}

	///	<summary>
	///		Perform one step of Newton iteration.
	///	</summary>
	///	<returns>
	///		The L2 norm of the difference between the initial state and
	///		updated state.
	///	</returns>
	double PerformBICGSTAB_NewtonStep(
		double * dX,
		int nMaxIter,
		double dTolerance
	);

	///	<summary>
	///		Perform one step of Newton iteration (safer version).
	///	</summary>
	///	<returns>
	///		The L2 norm of the difference between the initial state and
	///		updated state.
	///	</returns>
	inline double PerformBICGSTAB_NewtonStep_Safe(
		DataVector<double> & dX,
		int nMaxIter,
		double dTolerance
	) {
		// Check input size
		if (dX.GetRows() != m_dG.GetRows()) {
			_EXCEPTION2(
				"Input vector length mismatch: Expected length %i, received %i",
				m_dG.GetRows(), dX.GetRows());
		}
		
		return PerformBICGSTAB_NewtonStep(dX, nMaxIter, dTolerance);
	}

public:
	///	<summary>
	///		Function to evaluate.
	///	</summary>
	virtual void Evaluate(
		const double * dX,
		double * dF
	) = 0;

private:
	///	<summary>
	///		Flag indicating initialization of the solver.
	///	</summary>
	bool m_fInitialized;

	///	<summary>
	///		Number of equations to solve.
	///	</summary>
	int m_nEquationCount;

	///	<summary>
	///		Number of iterations to perform before a GMRES restart.
	///	</summary>
	int m_nIterPerRestart;

	///	<summary>
	///		Perturbation size.
	///	</summary>
	double m_dEpsilon;

	///	<summary>
	///		GMRES workspace variables.
	///	</summary>
	DataVector<double> m_dG;
	DataVector<double> m_dFX;
	DataMatrix<double> m_dV;
	DataVector<double> m_dH;
	DataVector<double> m_dCS;
	DataVector<double> m_dSN;
	DataVector<double> m_dS;
	DataVector<double> m_dT;
	DataVector<double> m_dY;
	DataVector<double> m_dW;
	DataVector<double> m_dPertX;
	DataVector<double> m_dR;
};

///////////////////////////////////////////////////////////////////////////////

#endif

