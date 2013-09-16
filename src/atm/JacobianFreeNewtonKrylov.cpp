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

#include "JacobianFreeNewtonKrylov.h"

#include "LinearAlgebra.h"

#include <cmath>

///////////////////////////////////////////////////////////////////////////////

JacobianFreeNewtonKrylov::JacobianFreeNewtonKrylov() :
	m_fInitialized(false)
{ }

///////////////////////////////////////////////////////////////////////////////

JacobianFreeNewtonKrylov::JacobianFreeNewtonKrylov(
	int nEquationCount,
	int nIterPerRestart,
	double dEpsilon
) {
	InitializeJFNK(nEquationCount, nIterPerRestart, dEpsilon);
}

///////////////////////////////////////////////////////////////////////////////

void JacobianFreeNewtonKrylov::InitializeJFNK(
	int nEquationCount,
	int nIterPerRestart,
	double dEpsilon
) {
	if (m_fInitialized) {
		_EXCEPTIONT("JacobianFreeNewtonKrylov already initialized");
	}

	if (nIterPerRestart < nEquationCount) {
		_EXCEPTION2("Number of iterations per restart (%i) must be at least "
			"equal to the number of equations (%i).",
				nIterPerRestart,
				nEquationCount);
	}

	// Initialize parameters
	m_nEquationCount = nEquationCount;
	m_nIterPerRestart = nIterPerRestart;
	m_dEpsilon = dEpsilon;

	// Initialize buffer arrays
	m_dG.Initialize(m_nEquationCount);

	m_dFX.Initialize(m_nEquationCount);

	m_dV.Initialize(m_nIterPerRestart+1, m_nEquationCount);

	m_dH.Initialize(m_nIterPerRestart * (m_nIterPerRestart + 1) / 2);

	m_dCS.Initialize(m_nIterPerRestart);

	m_dSN.Initialize(m_nIterPerRestart);

	m_dS.Initialize(m_nIterPerRestart);

	m_dT.Initialize(m_nIterPerRestart);

	m_dY.Initialize(m_nEquationCount);

	m_dW.Initialize(m_nEquationCount);

	m_dPertX.Initialize(m_nEquationCount);

	m_dR.Initialize(m_nEquationCount);

	m_fInitialized = true;
}

///////////////////////////////////////////////////////////////////////////////

void JacobianFreeNewtonKrylov::GivensRotationMatrix(
	double dA,
	double dB,
	double & dC,
	double & dS
) const {
	double dTemp;

	if (dB == 0.0) {
		dC = 1.0;
		dS = 0.0;

	} else if (fabs(dB) > fabs(dA)) {
		dTemp = dA / dB;
		dS = 1.0 / sqrt(1.0 + dTemp * dTemp);
		dC = dTemp * dS;

	} else {
		dTemp = dB / dA;
		dC = 1.0 / sqrt(1.0 + dTemp * dTemp);
		dS = dTemp * dC;
	}
}

///////////////////////////////////////////////////////////////////////////////

bool JacobianFreeNewtonKrylov::PerformJFNK(
	double * dX,
	int nMaxIter,
	double dSolnAbsTolerance,
	int nGMRESMaxIter,
	double dGMRESTolerance
) {
	// Check that JFNK is initialized
	if (!m_fInitialized) {
		_EXCEPTIONT(
			"Attempting to run uninitialized JacobianFreeNewtonKrylov");
	}

	// Loop over all iterations
	int k;
	for (k = 0; k < nMaxIter; k++) {

		// Perform one Newton step
		double dGNorm =
			PerformJFNK_NewtonStep(dX, nGMRESMaxIter, dGMRESTolerance);

		// Check the step difference against the tolerance
		if (dGNorm < dSolnAbsTolerance) {
			break;
		}
	}

	if (k == nMaxIter) {
		return false;
	}
	return true;
}

///////////////////////////////////////////////////////////////////////////////

double JacobianFreeNewtonKrylov::PerformJFNK_NewtonStep(
	double * dX,
	int nMaxIter,
	double dTolerance
) {
	// Check that JFNK is initialized
	if (!m_fInitialized) {
		_EXCEPTIONT(
			"Attempting to run uninitialized JacobianFreeNewtonKrylov");
	}

	// Info from LAPACK
	int iInfo;

	// Buffer variables
	double dTemp;
	double dError;

	// Set G to zero
	m_dG.Zero();

	// Total state vector size
	unsigned int nN = m_dG.GetRows();

	// Calculate the base G vector (RHS)
	m_dFX.Zero();
	Evaluate(dX, m_dFX);

	// Calculate the norm of the RHS
	m_dR = m_dFX;

	double dRNorm = LAPACK::DNORM2(m_dFX);
	if (fabs(dRNorm) < 1.0e-14) {
		return (0.0);
	}

	double dBNorm = dRNorm;

	// Iterate
	for (int iIter = 0; iIter < nMaxIter; iIter++) {

		// Reset H
		m_dH.Zero();

		// Store the residual in V(:,1)
		LAPACK::DSCAL(m_dR, 1.0 / dRNorm);
		LAPACK::DCOPY_A(m_dR, &(m_dV[0][0]));

		// Construct the orthonormal basis using Gram-Schmidt
		m_dS.Zero();
		m_dS[0] = dRNorm;

		dError = m_dS[0] / dBNorm;
		if (dError <= dTolerance) {
			break;
		}

		// Loop
		int u0i = 0;

		for (int i = 0; i < m_nIterPerRestart; i++) {

			// w = A*V(:,i);
			memcpy(m_dPertX, dX, m_dPertX.GetRows() * sizeof(double));
			LAPACK::DAXPY_A(m_dEpsilon, &(m_dV[i][0]), m_dPertX);

			m_dW.Zero();
			Evaluate(m_dPertX, m_dW);
			for (int j = 0; j < nN; j++) {
				m_dW[j] = (m_dW[j] - m_dFX[j]) / m_dEpsilon;
			}

			// Apply Gram-Schmidt
			for (int k = 0; k <= i; k++) {
				// H(k,i) = w'*V(:,k);
				m_dH[u0i + k] = LAPACK::DDOT_A(m_dW, &(m_dV[k][0]));

				// w = w - H(k,i)*V(:,k);
				LAPACK::DAXPY_A( -m_dH[u0i + k], &(m_dV[k][0]), m_dW);
			}

			// H(i+1,i) = norm( w );
			double dHx = LAPACK::DNORM2(m_dW);

			// V(:,i+1) = w / H(i+1,i);
			LAPACK::DSCAL(m_dW, 1.0 / dHx);
			LAPACK::DCOPY_A(m_dW, &(m_dV[i+1][0]));

			// Apply Givens rotation
			for (int k = 0; k <= i-1; k++) {
				int uki = u0i + k;

				dTemp       =   m_dCS[k] * m_dH[uki] + m_dSN[k] * m_dH[uki+1];
				m_dH[uki+1] = - m_dSN[k] * m_dH[uki] + m_dCS[k] * m_dH[uki+1];
				m_dH[uki]   =   dTemp;
			}

			// Form i-th rotation matrix, approximate residual norm

			// [cs(i),sn(i)] = rotmat( H(i,i), H(i+1,i) );
			GivensRotationMatrix(m_dH[u0i+i], dHx, m_dCS[i], m_dSN[i]);

			dTemp     =   m_dCS[i] * m_dS[i];
			m_dS[i+1] = - m_dSN[i] * m_dS[i];
			m_dS[i]   =   dTemp;

			// H(i,i) = cs(i)*H(i,i) + sn(i)*H(i+1,i);
			m_dH[u0i+i] = m_dCS[i] * m_dH[u0i+i] + m_dSN[i] * dHx;

			// error  = abs(s(i+1)) / bnrm2;
			dError = fabs(m_dS[i+1]) / dBNorm;

			// Update approximation and exit
			if (dError <= dTolerance) {
				m_dY = m_dS;

				iInfo = LAPACK::DTPSV('U', 'N', 'N', i+1, m_dH, m_dY);
				if (iInfo != 0) {
					_EXCEPTION1(
						"LAPACK error (%d)  No matrix solution found.", iInfo);
				}

				for (int j = 0; j <= i; j++) {
					LAPACK::DAXPY_A(m_dY[j], &(m_dV[j][0]), m_dG);
				}

				break;
			}

			// Increment upper triangular index counter
			u0i += (i + 1);
		}

		if (dError <= dTolerance) {
			break;
		}

		// y = H(1:m,1:m) \ s(1:m);
		m_dY = m_dS;

		iInfo = LAPACK::DTPSV('U', 'N', 'N', m_nIterPerRestart, m_dH, m_dY);
		if (iInfo != 0) {
			_EXCEPTION1(
				"LAPACK error (%d)  No matrix solution found.", iInfo);
		}

		for (int j = 0; j <= m_nIterPerRestart; j++) {
			LAPACK::DAXPY_A(m_dY[j], &(m_dV[j][0]), m_dG);
		}

		// Calculate the perturbed state
		memcpy(m_dPertX, dX, m_dPertX.GetRows() * sizeof(double));
		LAPACK::DAXPY_A(m_dEpsilon, m_dG, m_dPertX);

		// Calculate the residual (R = B-A*G)
		Evaluate(m_dPertX, m_dR);
		for (int j = 0; j < nN; j++) {
			m_dR[j] = m_dFX[j] - (m_dR[j] - m_dFX[j]) / m_dEpsilon;
		}

		dRNorm = LAPACK::DNORM2(m_dR);

		// error = s(i+1) / bnrm2;
		dError = dRNorm / dBNorm;

		// Check error
		if (dError <= dTolerance) {
			break;
		}
	}

	// Check for convergence
	if (dError > dTolerance) {
		_EXCEPTIONT("Convergence failure in GMRES.");
	}

	// Perturbation size
	//double dGNorm = LAPACK::DNORM2(m_dG);

	// Update the state
	for (int j = 0; j < m_dG.GetRows(); j++) {
		dX[j] -= m_dG[j];
	}

	return dError;
}

///////////////////////////////////////////////////////////////////////////////

double JacobianFreeNewtonKrylov::PerformBICGSTAB_NewtonStep(
	double * dX,
	int nMaxIter,
	double dTolerance
) {
	// Check that JFNK is initialized
	if (!m_fInitialized) {
		_EXCEPTIONT(
			"Attempting to run uninitialized JacobianFreeNewtonKrylov");
	}

	// Info from LAPACK
	int iInfo;

	// Buffer variables
	double dRhoOld = 1.0;
	double dRhoNew;
	double dAlpha = 1.0;
	double dBeta;
	double dOmega = 1.0;
	double dError;

	// Set G to zero
	m_dG.Zero();

	// Total state vector size
	unsigned int nN = m_dG.GetRows();

	// Calculate the base G vector (RHS)
	Evaluate(dX, m_dFX);

	// Calculate the norm of the RHS
	LAPACK::DCOPY(m_dFX, m_dR);

	// Iterate
	for (int i = 0; i < nMaxIter; i++) {

		// rho_i = (r[0], r[i-1])
		dRhoNew = LAPACK::DDOT(m_dFX, m_dR);

		// beta = (rho[i] / rho[i-1]) * (alpha / omega[i-1])
		dBeta = (dRhoNew / dRhoOld) * (dAlpha / dOmega);

		// p[i] = r[i-1] + beta * (p[i-1] - omega[i-1] * v[i-1])
		if (i != 0) {
			LAPACK::DAXPY(-dOmega, m_dY, m_dW);
			LAPACK::DSCAL(m_dW, dBeta);
			LAPACK::DAXPY(1.0, m_dR, m_dW);
		} else {
			LAPACK::DCOPY(m_dR, m_dW);
		}

		// v[i] = A p[i]
		memcpy(m_dPertX, dX, nN * sizeof(double));
		LAPACK::DAXPY(m_dEpsilon, m_dW, m_dPertX);
		Evaluate(m_dPertX, m_dY);
		for (int j = 0; j < nN; j++) {
			m_dY[j] = (m_dY[j] - m_dFX[j]) / m_dEpsilon;
		}

		// alpha = rho[i] / (r[0], v[i])
		dAlpha = dRhoNew / LAPACK::DDOT(m_dFX, m_dY);

		// s = r[i-1] - alpha * v[i]
		LAPACK::DCOPY(m_dR, m_dS);
		LAPACK::DAXPY(-dAlpha, m_dY, m_dS);

		// t = A s
		memcpy(m_dPertX, dX, nN * sizeof(double));
		LAPACK::DAXPY(m_dEpsilon, m_dS, m_dPertX);
		Evaluate(m_dPertX, m_dT);
		for (int j = 0; j < nN; j++) {
			m_dT[j] = (m_dT[j] - m_dFX[j]) / m_dEpsilon;
		}

		// omega[i] = (t,s)/(t,t)
		double dNormT = LAPACK::DDOT(m_dT, m_dT);
		if (dNormT == 0.0) {
			dOmega = 0.0;
		} else {
			dOmega = LAPACK::DDOT(m_dT, m_dS) / dNormT;
		}

		// x[i] = x[i-1] + alpha * p[i] + omega[i] * s
		LAPACK::DAXPY(dAlpha, m_dW, m_dG);
		LAPACK::DAXPY(dOmega, m_dS, m_dG);

		// r[i] = s - omega[i] * t
		LAPACK::DCOPY(m_dS, m_dR);
		LAPACK::DAXPY(-dOmega, m_dT, m_dR);

		dError = LAPACK::DDOT(m_dR, m_dR);

		// Check if the residual is less than the tolerance
		if (dError < dTolerance) {
			break;
		}

		// Update rho
		dRhoOld = dRhoNew;
	}

	// Update the state
	for (int j = 0; j < m_dG.GetRows(); j++) {
		dX[j] -= m_dG[j];
	}

	return dError;
}

///////////////////////////////////////////////////////////////////////////////

