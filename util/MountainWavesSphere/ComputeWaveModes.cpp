///////////////////////////////////////////////////////////////////////////////
///
///	\file    MountainWavesSphere.cpp
///	\author  Paul Ullrich
///	\version April 11, 2014
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

#include "DataVector.h"
#include "DataMatrix.h"
#include "CommandLine.h"
#include "Announce.h"
#include "Parameters.h"
#include "PolynomialInterp.h"

#include <vector>
#include <map>
#include <cmath>

#include <mpi.h>

#include <netcdfcpp.h>

#if defined USEVECLIB || defined USEMKL

extern "C" {
	void dggev_(char *jobvl, char *jobvr, int *n,  double  *a,  int
               *lda,  double  *b,  int *ldb, double *alphar, double
               *alphai, double *beta, double *vl, int *ldvl,
			   double *vr, int *ldvr, double *work, int *lwork, int *info);
}

#endif

///////////////////////////////////////////////////////////////////////////////

void GenerateEvolutionMatrix(
	int nK,
	const Parameters & param,
	DataMatrix<double> & matM,
	DataMatrix<double> & matB,
	double & dInvRo,
	double & dFr
) {
	int nMatrixSize = 5 * param.nPhiElements - 1;

	matM.Initialize(nMatrixSize, nMatrixSize);
	matB.Initialize(nMatrixSize, nMatrixSize);

	// Radius of the Earth
	const double ParamEarthRadius = 6.37122e6;

	// Ideal gas constant
	const double ParamRd = 287.0;

	// Inverse Rossby number
	dInvRo = 2.0 * ParamEarthRadius * param.dOmega * param.dXscale / param.dU0;

	// Scale height
	double dH = ParamRd * param.dT0 / param.dG;

	// Froude number
	dFr = param.dU0 / sqrt(param.dG * dH);

	// Squared Froude number
	double dFr2 = dFr * dFr;

	// Aspect ratio
	double dAs = dH / (ParamEarthRadius / param.dXscale);

	// Velocity ratio
	double dAv = dAs;

	// Wave number
	double dK2 = static_cast<double>(nK * nK);

	// Gamma
	double dInvGamma = 1.0 / param.dGamma;

	// Grid spacing
	double dDeltaPhi = param.vecNode[1] - param.vecNode[0];

	// Loop through all elements
	for (int j = 0; j < param.nPhiElements; j++) {

		// Index of the start of this block
		int ix = 4 * j;
		int ixU = ix;
		int ixP = ix + 1;
		int ixW = ix + 2;
		int ixR = ix + 3;

		int ixVL = 4 * param.nPhiElements + j - 1;
		int ixVR = 4 * param.nPhiElements + j;

		// Phi at this node
		double dPhi = param.vecNode[j];

		// Cosine Phi
		double dCosPhi = cos(dPhi);

		// Sin Phi
		double dSinPhi = sin(dPhi);

		// Tan Phi
		double dTanPhi = tan(dPhi);

		// U evolution equations
		matM[ixU][ixU] = dFr2 * dCosPhi * dCosPhi;
		matM[ixP][ixU] = 1.0;

		if (j != 0) {
			matM[ixVL][ixU] = - 0.5 * dFr2 * (2.0 + dInvRo) * dSinPhi * dCosPhi;
		}
		if (j != param.nPhiElements - 1) {
			matM[ixVR][ixU] = - 0.5 * dFr2 * (2.0 + dInvRo) * dSinPhi * dCosPhi;
		}

		// V evolution equations
		if (j != 0) {
			int ixV = ixVL;

			int ixUL = ix - 4;
			int ixPL = ix - 3;
			int ixRL = ix - 1;

			int ixUR = ix;
			int ixPR = ix + 1;
			int ixRR = ix + 3;

			double dPhiS = param.vecEdge[j];
			double dSinPhiS = sin(dPhiS);
			double dCosPhiS = cos(dPhiS);

			matM[ixUL][ixV] = 0.5 * dFr2 * (2.0 + dInvRo) * dSinPhiS * dCosPhiS;
			matM[ixUR][ixV] = 0.5 * dFr2 * (2.0 + dInvRo) * dSinPhiS * dCosPhiS;

			matM[ixV][ixV] = - dK2 * dFr2;

			matM[ixPL][ixV] =
				- 0.5 * dFr2 * (1.0 + dInvRo) * dSinPhiS * dCosPhiS
				- 1.0 / dDeltaPhi;
			matM[ixPR][ixV] =
				- 0.5 * dFr2 * (1.0 + dInvRo) * dSinPhiS * dCosPhiS
				+ 1.0 / dDeltaPhi;

			matM[ixRL][ixV] = 0.5 * dFr2 * (1.0 + dInvRo) * dSinPhiS * dCosPhiS;
			matM[ixRR][ixV] = 0.5 * dFr2 * (1.0 + dInvRo) * dSinPhiS * dCosPhiS;
		}

		// P evolution equations
		matM[ixU][ixP] = dCosPhi;
		matM[ixR][ixP] = dCosPhi;

		if (j != 0) {
			matM[ixVL][ixP] =
				- 0.5 * dFr2 * (1.0 + dInvRo) * dSinPhi * dCosPhi * dCosPhi
				- 0.5 * dSinPhi
				- dCosPhi / dDeltaPhi;
		}
		if (j != param.nPhiElements - 1) {
			matM[ixVR][ixP] =
				- 0.5 * dFr2 * (1.0 + dInvRo) * dSinPhi * dCosPhi * dCosPhi
				- 0.5 * dSinPhi
				+ dCosPhi / dDeltaPhi;
		}

		// W evolution equation
		matM[ixW][ixW] = - dK2 * dAs * dAv * dFr2;
		matM[ixR][ixW] = 1.0;

		// R evolution equation
		matM[ixP][ixR] = dInvGamma / (1.0 - dInvGamma);
		matM[ixW][ixR] = dAv / dAs;
		matM[ixR][ixR] = - 1.0 / (1.0 - dInvGamma);

		if (j != 0) {
			matM[ixVL][ixR] = 0.5 * dFr2 * (1.0 + dInvRo) * dSinPhi * dCosPhi;
		}
		if (j != param.nPhiElements - 1) {
			matM[ixVR][ixR] = 0.5 * dFr2 * (1.0 + dInvRo) * dSinPhi * dCosPhi;
		}

		// B matrix coefficients
		matB[ixP][ixW] = -1.0;
		matB[ixW][ixP] = -1.0;
	}
}

///////////////////////////////////////////////////////////////////////////////

void SolveEvolutionMatrix(
	DataMatrix<double> & matM,
	DataMatrix<double> & matB,
	DataVector<double> & vecAlphaR,
	DataVector<double> & vecAlphaI,
	DataVector<double> & vecBeta,
	DataMatrix<double> & matVR
) {
	int iInfo = 0;

	char jobvl = 'N';
	char jobvr = 'V';

	int n = matM.GetRows();
	int lda = matM.GetRows();
	int ldb = matM.GetRows();
	int ldvl = 1;
	int ldvr = matVR.GetRows();

	DataVector<double> vecWork;
	vecWork.Initialize(8 * n);
	int lwork = vecWork.GetRows();

	dggev_(
		&jobvl,
		&jobvr,
		&n,
		&(matM[0][0]),
		&lda,
		&(matB[0][0]),
		&ldb,
		&(vecAlphaR[0]),
		&(vecAlphaI[0]),
		&(vecBeta[0]),
		NULL,
		&ldvl,
		&(matVR[0][0]),
		&ldvr,
		&(vecWork[0]),
		&lwork,
		&iInfo);

	int nCount = 0;
	for (int i = 0; i < n; i++) {
		if (vecBeta[i] != 0.0) {
			//printf("%i %1.5e %1.5e\n", i, vecAlphaR[i] / vecBeta[i], vecAlphaI[i] / vecBeta[i]);
			nCount++;
		}
	}
/*
	for (int i = 0; i < 40; i++) {
		printf("%1.5e %1.5e\n", matVR[11][4*i+2], matVR[12][4*i+2]);
	}
*/
	Announce("%i total eigenvalues found", nCount);
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char ** argv) {

	MPI_Init(&argc, &argv);

try {

	// Parameters
	Parameters param;

	// Output filename
	std::string strOutputFile;

	// Horizontal minimum wave number
	int nKmin;

	// Horizontal maximum wave number
	int nKmax;

	// Parse the command line
	BeginCommandLine()
		CommandLineInt(param.nPhiElements, "n", 40);
		CommandLineInt(nKmin, "kmin", 1);
		CommandLineInt(nKmax, "kmax", 20);
		CommandLineDouble(param.dXscale, "X", 1.0);
		CommandLineDouble(param.dT0, "T0", 300.0);
		CommandLineDouble(param.dU0, "U0", 20.0);
		CommandLineDouble(param.dG, "G", 9.80616);
		CommandLineDouble(param.dOmega, "omega", 7.29212e-5);
		CommandLineDouble(param.dGamma, "gamma", 1.4);
		CommandLineString(strOutputFile, "out", "wave.nc");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Generate latitude values
	param.GenerateLatituteArray(param.nPhiElements);

	// Open NetCDF file
	NcFile ncdf_out(strOutputFile.c_str(), NcFile::Replace);

	NcDim *dimK = ncdf_out.add_dim("k", nKmax - nKmin + 1);
	NcDim *dimLat = ncdf_out.add_dim("lat", param.nPhiElements);
	NcDim *dimEig = ncdf_out.add_dim("eig", param.nPhiElements);

	// Write parameters and latitudes to file
	param.WriteToNcFile(ncdf_out, dimLat);

	// Wave numbers 
	NcVar *varK = ncdf_out.add_var("k", ncInt, dimK);

	DataVector<int> vecK;
	vecK.Initialize(nKmax - nKmin + 1);
	for (int nK = nKmin; nK <= nKmax; nK++) { 
		vecK[nK - nKmin] = nK;
	}

	varK->set_cur((long)0);
	varK->put(vecK, nKmax - nKmin + 1);

	// Eigenvalues
	NcVar *varMR = ncdf_out.add_var("mR", ncDouble, dimK, dimEig);
	NcVar *varMI = ncdf_out.add_var("mI", ncDouble, dimK, dimEig);

	NcVar *varUR = ncdf_out.add_var("uR", ncDouble, dimK, dimEig, dimLat);
	NcVar *varUI = ncdf_out.add_var("uI", ncDouble, dimK, dimEig, dimLat);

	NcVar *varVR = ncdf_out.add_var("vR", ncDouble, dimK, dimEig, dimLat);
	NcVar *varVI = ncdf_out.add_var("vI", ncDouble, dimK, dimEig, dimLat);

	NcVar *varPR = ncdf_out.add_var("pR", ncDouble, dimK, dimEig, dimLat);
	NcVar *varPI = ncdf_out.add_var("pI", ncDouble, dimK, dimEig, dimLat);

	NcVar *varWR = ncdf_out.add_var("wR", ncDouble, dimK, dimEig, dimLat);
	NcVar *varWI = ncdf_out.add_var("wI", ncDouble, dimK, dimEig, dimLat);

	NcVar *varRhoR = ncdf_out.add_var("rhoR", ncDouble, dimK, dimEig, dimLat);
	NcVar *varRhoI = ncdf_out.add_var("rhoI", ncDouble, dimK, dimEig, dimLat);

	// Allocate temporary arrays
	DataVector<double> dUR;
	dUR.Initialize(param.nPhiElements);

	DataVector<double> dUI;
	dUI.Initialize(param.nPhiElements);

	DataVector<double> dVR;
	dVR.Initialize(param.nPhiElements);

	DataVector<double> dVI;
	dVI.Initialize(param.nPhiElements);

	DataVector<double> dPR;
	dPR.Initialize(param.nPhiElements);

	DataVector<double> dPI;
	dPI.Initialize(param.nPhiElements);

	DataVector<double> dWR;
	dWR.Initialize(param.nPhiElements);

	DataVector<double> dWI;
	dWI.Initialize(param.nPhiElements);

	DataVector<double> dRhoR;
	dRhoR.Initialize(param.nPhiElements);

	DataVector<double> dRhoI;
	dRhoI.Initialize(param.nPhiElements);

	// Loop over all horizontal wave numbers
	for (int nK = nKmin; nK <= nKmax; nK++) {

		// Build matrices
		char szMessage[100];
		sprintf(szMessage, "Building evolution matrices (k = %i)", nK);
		AnnounceStartBlock(szMessage);

		DataMatrix<double> matM;
		DataMatrix<double> matB;

		double dInvRo;
		double dFr;
		GenerateEvolutionMatrix(nK, param, matM, matB, dInvRo, dFr);

		if (nK == nKmin) {
			ncdf_out.add_att("InvRo", dInvRo);
			ncdf_out.add_att("Fr", dFr);
		}

		AnnounceEndBlock("Done");

		// Solve the matrices
		AnnounceStartBlock("Solving evolution matrices");

		DataVector<double> vecAlphaR;
		DataVector<double> vecAlphaI;
		DataVector<double> vecBeta;
		DataMatrix<double> matVR;

		vecAlphaR.Initialize(matM.GetRows());
		vecAlphaI.Initialize(matM.GetRows());
		vecBeta  .Initialize(matM.GetRows());
		matVR    .Initialize(matM.GetRows(), matM.GetColumns());

		SolveEvolutionMatrix(
			matM,
			matB,
			vecAlphaR,
			vecAlphaI,
			vecBeta,
			matVR);

		// Sort eigenvalues
		std::multimap<double, int> mapSortedRows;
		for (int i = 0; i < vecBeta.GetRows(); i++) {
			if (vecBeta[i] != 0.0) {

				double dLambdaR = vecAlphaR[i] / vecBeta[i];
				double dLambdaI = vecAlphaI[i] / vecBeta[i];

				double dMR = dLambdaI;
				double dMI = - dLambdaR - 1.0;

				double dEvalMagnitude = fabs(dMR);

				mapSortedRows.insert(
					std::pair<double, int>(dEvalMagnitude, i));

				// Only store one entry for complex-conjugate pairs
				if (vecAlphaI[i] != 0.0) {
					i++;
				}
			}
		}
/*
		Announce("%i eigenmodes found to satisfy entropy condition",
			mapSortedRows.size());
*/
/*
		if (mapSortedRows.size() != param.nPhiElements) {
			_EXCEPTIONT("Mismatch between eigenmode count and latitude count");
		}
*/
		AnnounceEndBlock("Done");

		// Write the matrices to a file
		AnnounceStartBlock("Writing results");

		int iKix = nK - nKmin;
		int iWave = 0;
		std::map<double, int>::const_iterator it;
		for (it = mapSortedRows.begin(); it != mapSortedRows.end(); it++) {
			int i = it->second;

			double dLambdaR = vecAlphaR[i] / vecBeta[i];
			double dLambdaI = vecAlphaI[i] / vecBeta[i];

			double dMR = dLambdaI;
			double dMI = - dLambdaR - 1.0;

			//printf("%1.5e %1.5e", dMR, dMI);

			// Real eigenmode must decay with height
			if ((dLambdaI == 0.0) && (dMI < -1.0e-9)) {
				//printf("\n");
				continue;
			} else {
				//printf(" (OK)\n");
			}

			// Store real part of eigenfunction
			for (int j = 0; j < param.nPhiElements; j++) {
				dUR  [j] = matVR[i][4*j  ];
				dPR  [j] = matVR[i][4*j+1];
				dWR  [j] = matVR[i][4*j+2];
				dRhoR[j] = matVR[i][4*j+3];
			}

			dVR.Zero();
			for (int j = 0; j < param.nPhiElements-1; j++) {
				dVR[j  ] += 0.5 * matVR[i][4 * param.nPhiElements + j];
				dVR[j+1] += 0.5 * matVR[i][4 * param.nPhiElements + j];
			}

			// Complex eigenvalue / eigenfunction pair
			// Only store eigenfunction with positive root
			if (dLambdaI != 0.0) {

				// Dump eigenvalue to NetCDF file
				if (dMR > 0.0) {
					varMR->set_cur(iKix, iWave);
					varMR->put(&dMR, 1, 1);

					// Store imaginary component of vector
					for (int j = 0; j < param.nPhiElements; j++) {
						dUI  [j] = matVR[i+1][4*j  ];
						dPI  [j] = matVR[i+1][4*j+1];
						dWI  [j] = matVR[i+1][4*j+2];
						dRhoI[j] = matVR[i+1][4*j+3];
					}

					dVI.Zero();
					for (int j = 0; j < param.nPhiElements-1; j++) {
						dVI[j  ] += 0.5 * matVR[i+1][4 * param.nPhiElements + j];
						dVI[j+1] += 0.5 * matVR[i+1][4 * param.nPhiElements + j];
					}

				} else {
					dMR = - dMR;

					varMR->set_cur(iKix, iWave);
					varMR->put(&dMR, 1, 1);

					// Store imaginary component of vector
					for (int j = 0; j < param.nPhiElements; j++) {
						dUI  [j] = - matVR[i+1][4*j  ];
						dPI  [j] = - matVR[i+1][4*j+1];
						dWI  [j] = - matVR[i+1][4*j+2];
						dRhoI[j] = - matVR[i+1][4*j+3];
					}

					dVI.Zero();
					for (int j = 0; j < param.nPhiElements-1; j++) {
						dVI[j  ] -= 0.5 * matVR[i+1][4 * param.nPhiElements + j];
						dVI[j+1] -= 0.5 * matVR[i+1][4 * param.nPhiElements + j];
					}

				}

				varMI->set_cur(iKix, iWave);
				varMI->put(&dMI, 1, 1);

			// Real eigenvalue / eigenvector pair
			} else {
				dUI.Zero();
				dPI.Zero();
				dWI.Zero();
				dRhoI.Zero();
				dVI.Zero();

				varMR->set_cur(iKix, iWave);
				varMR->put(&dMR, 1, 1);

				varMI->set_cur(iKix, iWave);
				varMI->put(&dMI, 1, 1);
			}

			// Dump the eigenfunction to the file
			varUR->set_cur(iKix, iWave, 0);
			varUR->put(dUR, 1, 1, param.nPhiElements);

			varVR->set_cur(iKix, iWave, 0);
			varVR->put(dVR, 1, 1, param.nPhiElements);

			varPR->set_cur(iKix, iWave, 0);
			varPR->put(dPR, 1, 1, param.nPhiElements);

			varWR->set_cur(iKix, iWave, 0);
			varWR->put(dWR, 1, 1, param.nPhiElements);

			varRhoR->set_cur(iKix, iWave, 0);
			varRhoR->put(dRhoR, 1, 1, param.nPhiElements);

			varUI->set_cur(iKix, iWave, 0);
			varUI->put(dUI, 1, 1, param.nPhiElements);

			varVI->set_cur(iKix, iWave, 0);
			varVI->put(dVI, 1, 1, param.nPhiElements);

			varPI->set_cur(iKix, iWave, 0);
			varPI->put(dPI, 1, 1, param.nPhiElements);

			varWI->set_cur(iKix, iWave, 0);
			varWI->put(dWI, 1, 1, param.nPhiElements);

			varRhoI->set_cur(iKix, iWave, 0);
			varRhoI->put(dRhoI, 1, 1, param.nPhiElements);

			// Increment wave index
			iWave++;
		}

		Announce("%i eigenmodes found to satisfy entropy condition", iWave);

		AnnounceEndBlock("Done");

		AnnounceBanner();
	}

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}
	MPI_Finalize();
}


///////////////////////////////////////////////////////////////////////////////

