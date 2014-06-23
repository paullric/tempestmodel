///////////////////////////////////////////////////////////////////////////////
///
///	\file    GenerateMountainResponse.cpp
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
#include "DataMatrix3D.h"
#include "CommandLine.h"
#include "Announce.h"
#include "Parameters.h"

#include <map>

#include <netcdfcpp.h>

#include "LinearAlgebra.h"

/*
#if defined USEVECLIB || defined USEMKL

extern "C" {
	int dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);
}

#endif
*/
///////////////////////////////////////////////////////////////////////////////

int main(int argc, char ** argv) {

try {

/*
	// Horizontal minimum wave number
	int nKmin;

	// Horizontal maximum wave number
	int nKmax;
*/
	// Number of vertical levels
	int nLev;

	// Model cap
	double dZtop;

	// Input filename
	std::string strWaveFile;

	// Topography filename
	std::string strTopoFile;

	// Output filename
	std::string strOutputFile;

	// Parse the command line
	BeginCommandLine()
		//CommandLineInt(nKmin, "kmin", -1);
		//CommandLineInt(nKmax, "kmax", -1);
		CommandLineInt(nLev, "levels", 30);
		CommandLineDouble(dZtop, "ztop", 12000.0);
		CommandLineString(strWaveFile, "wave", "wave.nc");
		CommandLineString(strTopoFile, "topo", "topo.nc");
		CommandLineString(strOutputFile, "out", "out.nc");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Open NetCDF topography file
	AnnounceStartBlock("Loading topography file");
	NcFile ncdf_topo(strTopoFile.c_str(), NcFile::ReadOnly);

	// Load latitude array
	NcVar *varTopoLat = ncdf_topo.get_var("lat");
	int nLat = varTopoLat->get_dim(0)->size();

	DataVector<double> dLat;
	dLat.Initialize(nLat);
	varTopoLat->set_cur((long)0);
	varTopoLat->get(&(dLat[0]), nLat);

	// Load longitude array
	NcVar *varTopoLon = ncdf_topo.get_var("lon");
	int nLon = varTopoLon->get_dim(0)->size();

	DataVector<double> dLon;
	dLon.Initialize(nLon);
	varTopoLon->set_cur((long)0);
	varTopoLon->get(&(dLon[0]), nLon);

	// Load topography
	NcVar *varZs = ncdf_topo.get_var("Zs");
	DataMatrix<double> dZs;
	dZs.Initialize(nLat, nLon);
	varZs->set_cur((long)0, (long)0);
	varZs->get(&(dZs[0][0]), nLat, nLon);

	ncdf_topo.close();

	// Generate level array
	DataVector<double> dLev;
	dLev.Initialize(nLev);
	for (int l = 0; l < nLev; l++) {
		dLev[l] = dZtop / static_cast<double>(nLev-1)
			* static_cast<double>(l);
	}

	AnnounceEndBlock("Done");

	// Open NetCDF wave file
	AnnounceStartBlock("Processing wave file");
	NcFile ncdf_wave(strWaveFile.c_str(), NcFile::ReadOnly);

	// Get latitude count
	NcVar *varWaveLat = ncdf_wave.get_var("lat");
	int nWaveLat = varWaveLat->get_dim(0)->size();

	if (nWaveLat != nLat) {
		_EXCEPTION4("Mismatch in number of latitudes:\n"
			"Topography file \"%s\" (%i)\n"
			"Wave file \"%s\" (%i)",
			strTopoFile.c_str(), nLat,
			strWaveFile.c_str(), nWaveLat);
	}

	// Get parameters
	Parameters param;
	param.ReadFromNcFile(ncdf_wave);

	// Ideal gas constant
	const double ParamRd = 287.0;

	// Scale height
	double dH = ParamRd * param.dT0 / param.dG;

	// Non-dimensionalize the topography
	for (int j = 0; j < nLat; j++) {
	for (int i = 0; i < nLon; i++) {
		dZs[j][i] /= dH;
	}
	}

	// Get the wave number vector
	NcVar *varK = ncdf_wave.get_var("k");
	NcDim *dimK = ncdf_wave.get_dim("k");

	DataVector<int> vecK;
	vecK.Initialize(dimK->size());

	varK->get(&(vecK[0]), dimK->size());

	// Grid spacing in longitude
	double dLonIntCoeff = (dLon[1] - dLon[0]) / (2.0 * M_PI);

	// Wave mode arrays
	DataMatrix<double> dUWave;
	DataMatrix<double> dPWave;
	DataMatrix<double> dWWave;
	DataMatrix<double> dRhoWave;
	DataMatrix<double> dVWave;

	dUWave.Initialize(nLat, 2 * nLat);
	dVWave.Initialize(nLat, 2 * nLat);
	dPWave.Initialize(nLat, 2 * nLat);
	dWWave.Initialize(nLat, 2 * nLat);
	dRhoWave.Initialize(nLat, 2 * nLat);

	// Prognostic variable arrays
	DataMatrix3D<double> dU;
	DataMatrix3D<double> dV;
	DataMatrix3D<double> dW;
	DataMatrix3D<double> dP;
	DataMatrix3D<double> dRho;

	dU  .Initialize(nLev, nLat, nLon);
	dV  .Initialize(nLev, nLat, nLon);
	dW  .Initialize(nLev, nLat, nLon);
	dP  .Initialize(nLev, nLat, nLon);
	dRho.Initialize(nLev, nLat, nLon);

	// Loop through desired wave numbers
	for (int k = 0; k < vecK.GetRows(); k++) {
	//for (int k = 4; k <= 4; k++) {

		// Current value of K
		int nK = vecK[k];
		double dK = static_cast<double>(nK);

		Announce("Processing wave number %i", nK);

		// Perform Fourier transform of topography
		DataVector<double> dFZs;
		dFZs.Initialize(2 * nLat);

		for (int j = 0; j < nLat; j++) {
		for (int i = 0; i < nLon; i++) {
			dFZs[     j] += dLonIntCoeff * dZs[j][i] * cos(dK * dLon[i]);
			dFZs[nLat+j] -= dLonIntCoeff * dZs[j][i] * sin(dK * dLon[i]);
		}
		}

		// Number of eigenvalues
		int nEig = nLat;

		// Load in eigenvalue
		DataVector<double> dMR;
		dMR.Initialize(nEig);

		DataVector<double> dMI;
		dMI.Initialize(nEig);

		NcVar *varMR = ncdf_wave.get_var("mR");
		NcVar *varMI = ncdf_wave.get_var("mI");

		varMR->set_cur(k, 0);
		varMR->get(&(dMR[0]), 1, nEig);

		varMI->set_cur(k, 0);
		varMI->get(&(dMI[0]), 1, nEig);

		// Load in transformed vertical velocity for boundary condition
		NcVar *varUR = ncdf_wave.get_var("uR");
		NcVar *varUI = ncdf_wave.get_var("uI");

		NcVar *varVR = ncdf_wave.get_var("vR");
		NcVar *varVI = ncdf_wave.get_var("vI");

		NcVar *varPR = ncdf_wave.get_var("pR");
		NcVar *varPI = ncdf_wave.get_var("pI");

		NcVar *varWR = ncdf_wave.get_var("wR");
		NcVar *varWI = ncdf_wave.get_var("wI");

		NcVar *varRhoR = ncdf_wave.get_var("rhoR");
		NcVar *varRhoI = ncdf_wave.get_var("rhoI");

		for (int n = 0; n < nEig; n++) {

			varUR->set_cur(k, n, 0);
			varUR->get(&(dUWave[n][0]), 1, 1, nLat);

			varUI->set_cur(k, n, 0);
			varUI->get(&(dUWave[n][nLat]), 1, 1, nLat);

			varVR->set_cur(k, n, 0);
			varVR->get(&(dVWave[n][0]), 1, 1, nLat);

			varVI->set_cur(k, n, 0);
			varVI->get(&(dVWave[n][nLat]), 1, 1, nLat);

			varPR->set_cur(k, n, 0);
			varPR->get(&(dPWave[n][0]), 1, 1, nLat);

			varPI->set_cur(k, n, 0);
			varPI->get(&(dPWave[n][nLat]), 1, 1, nLat);

			varWR->set_cur(k, n, 0);
			varWR->get(&(dWWave[n][0]), 1, 1, nLat);

			varWI->set_cur(k, n, 0);
			varWI->get(&(dWWave[n][nLat]), 1, 1, nLat);

			varRhoR->set_cur(k, n, 0);
			varRhoR->get(&(dRhoWave[n][0]), 1, 1, nLat);

			varRhoI->set_cur(k, n, 0);
			varRhoI->get(&(dRhoWave[n][nLat]), 1, 1, nLat);
		}

		// Combined W matrix
		DataMatrix<double> dWWaveCombined;
		dWWaveCombined.Initialize(2 * nEig, 2 * nEig);
		for (int n = 0; n < nEig; n++) {
		for (int j = 0; j < nEig; j++) {
			dWWaveCombined[     n][     j] = + dWWave[n][     j];
			dWWaveCombined[     n][nEig+j] = + dWWave[n][nEig+j];
			dWWaveCombined[nEig+n][     j] = - dWWave[n][nEig+j];
			dWWaveCombined[nEig+n][nEig+j] = + dWWave[n][     j];
		}
		}
/*
		FILE * fp = fopen("w.dat", "w");
		for (int i = 0; i < 2*nEig; i++) {
			for (int j = 0; j < 2*nEig; j++) {
				fprintf(fp, "%1.15e ", dWWaveCombined[i][j]);
			}
			fprintf(fp, "\n");
		}
		fclose(fp);

		FILE * fpwr = fopen("wr.dat", "w");
		FILE * fpwi = fopen("wi.dat", "w");
		for (int i = 0; i < nEig; i++) {
			for (int j = 0; j < nEig; j++) {
				fprintf(fpwr, "%1.15e ", dWWave[i][j]);
				fprintf(fpwi, "%1.15e ", dWWave[i][nEig+j]);
			}
			fprintf(fpwr, "\n");
			fprintf(fpwi, "\n");
		}
		fclose(fpwr);
		fclose(fpwi);

		FILE * fph = fopen("h.dat", "w");
		for (int i = 0; i < 2 * nEig; i++) {
			fprintf(fph, "%1.15e\n", dFZs[i]);
		}
		fclose(fph);
*/
/*
		printf("%i\n", nInfo);
		for (int j = 0; j < nEig; j++) {
			printf("B: %1.10e\n", dFZs[j]);
		}
		_EXCEPTION();
*/

		// Solve the matrix system
		int nDim = 2 * nEig;
		int nRHS = 1;
		int nLDA = 2 * nEig;
		int nLDB = 2 * nEig;
		int nInfo;

		DataVector<double> dC;
		dC = dFZs;

		DataVector<int> iPIV;
		iPIV.Initialize(2 * nEig, 2 * nEig);

		dgesv_(
			&nDim, &nRHS,
			&(dWWaveCombined[0][0]), &nLDA,
			&(iPIV[0]),
			&(dC[0]), &nLDB,
			&nInfo);

		if (nInfo != 0) {
			_EXCEPTION1("DGESV exited with message: %i\n", nInfo);
		}

/*
		// Solve the matrix system using the generalized inverse
		DataMatrix<double> dWWaveInv;
		LAPACK::GeneralizedInverse(dWWaveCombined, dWWaveInv);

		FILE * fpwinv = fopen("winv.dat", "w");
		for (int i = 0; i < 2*nEig; i++) {
			for (int j = 0; j < 2*nEig; j++) {
				fprintf(fpwinv, "%1.15e ", dWWaveInv[i][j]);
			}
			fprintf(fpwinv, "\n");
		}
		fclose(fpwinv);

		DataVector<double> dC;
		dC.Initialize(2 * nEig);

		for (int i = 0; i < 2 * nEig; i++) {
		for (int j = 0; j < 2 * nEig; j++) {
			dC[j] += dWWaveInv[i][j] * dFZs[i];
		}
		}
*/
/*
		for (int i = 0; i < nEig; i++) {
			printf("%1.5e %1.5e\n", dC[i], dC[i+nEig]);
		}
*/
		// Add contributions to the array
		for (int n = 0; n < nEig; n++) {

			// Skip rapidly growing wave modes
			if (dMI[n] < -1.0) {
				continue;
			}

			for (int l = 0; l < nLev; l++) {
			for (int j = 0; j < nLat; j++) {
			for (int i = 0; i < nLon; i++) {
				double dLambda = dLon[i];
				double dPhi = dLat[j];
				double dZ = dLev[l] / dH;

				double dPsi = dK * dLambda + dMR[n] * dZ;

				double dAmp = exp(- dMI[n] * dZ);

				// Zonal velocity
				dU[l][j][i] += dC[n] * cos(dPhi) * dAmp * (
					+ dUWave[n][     j] * cos(dPsi)
					- dUWave[n][nLat+j] * sin(dPsi));

				dU[l][j][i] -= dC[nEig+n] * cos(dPhi) * dAmp * (
					+ dUWave[n][     j] * sin(dPsi)
					+ dUWave[n][nLat+j] * cos(dPsi));

				// Meridional velocity
				dV[l][j][i] -= dC[n] * dK * dAmp * (
					+ dVWave[n][     j] * sin(dPsi)
					+ dVWave[n][nLat+j] * cos(dPsi));

				dV[l][j][i] += dC[nEig+n] * dK * dAmp * (
					- dVWave[n][     j] * cos(dPsi)
					+ dVWave[n][nLat+j] * sin(dPsi));

				// Pressure
				dP[l][j][i] += dC[n] * dAmp * (
					+ dPWave[n][     j] * cos(dPsi)
					- dPWave[n][nLat+j] * sin(dPsi));

				dP[l][j][i] -= dC[nEig+n] * dAmp * (
					+ dPWave[n][     j] * sin(dPsi)
					+ dPWave[n][nLat+j] * cos(dPsi));

				// Vertical velocity
				dW[l][j][i] -= dC[n] * dK * dAmp * (
					+ dWWave[n][     j] * sin(dPsi)
					+ dWWave[n][nLat+j] * cos(dPsi));

				dW[l][j][i] += dC[nEig+n] * dK * dAmp * (
					- dWWave[n][     j] * cos(dPsi)
					+ dWWave[n][nLat+j] * sin(dPsi));

				// Density
				dRho[l][j][i] += dC[n] * dAmp * (
					+ dRhoWave[n][     j] * cos(dPsi)
					- dRhoWave[n][nLat+j] * sin(dPsi));

				dRho[l][j][i] -= dC[nEig+n] * dAmp * (
					+ dRhoWave[n][     j] * sin(dPsi)
					+ dRhoWave[n][nLat+j] * cos(dPsi));

			}
			}
			}
		}
	}

	AnnounceEndBlock("Done!");

	// Output variables
	AnnounceStartBlock("Writing results");

	NcFile ncdf_out(strOutputFile.c_str(), NcFile::Replace);

	// Output dimensions and corresponding variables
	NcDim * dimOutLat = ncdf_out.add_dim("lat", nLat);
	NcDim * dimOutLon = ncdf_out.add_dim("lon", nLon);
	NcDim * dimOutLev = ncdf_out.add_dim("lev", nLev);

	NcVar * varOutLat = ncdf_out.add_var("lat", ncDouble, dimOutLat);
	NcVar * varOutLon = ncdf_out.add_var("lon", ncDouble, dimOutLon);
	NcVar * varOutLev = ncdf_out.add_var("lev", ncDouble, dimOutLev);

	varOutLat->set_cur((long)0);
	varOutLat->put(&(dLat[0]), nLat);

	varOutLon->set_cur((long)0);
	varOutLon->put(&(dLon[0]), nLon);

	varOutLev->set_cur((long)0);
	varOutLev->put(&(dLev[0]), nLev);

	// Output topography
	NcVar * varZS = ncdf_out.add_var("zs", ncDouble, dimOutLat, dimOutLon);
	varZS->set_cur(0, 0);
	varZS->put(&(dZs[0][0]), nLat, nLon);

	// Output variables
	NcVar * varOutU =
		ncdf_out.add_var("u", ncDouble, dimOutLev, dimOutLat, dimOutLon);
	NcVar * varOutV =
		ncdf_out.add_var("v", ncDouble, dimOutLev, dimOutLat, dimOutLon);
	NcVar * varOutP =
		ncdf_out.add_var("p", ncDouble, dimOutLev, dimOutLat, dimOutLon);
	NcVar * varOutW =
		ncdf_out.add_var("w", ncDouble, dimOutLev, dimOutLat, dimOutLon);
	NcVar * varOutRho =
		ncdf_out.add_var("rho", ncDouble, dimOutLev, dimOutLat, dimOutLon);

	varOutU->set_cur(0, 0, 0);
	varOutU->put(&(dU[0][0][0]), nLev, nLat, nLon);

	varOutV->set_cur(0, 0, 0);
	varOutV->put(&(dV[0][0][0]), nLev, nLat, nLon);

	varOutP->set_cur(0, 0, 0);
	varOutP->put(&(dP[0][0][0]), nLev, nLat, nLon);

	varOutW->set_cur(0, 0, 0);
	varOutW->put(&(dW[0][0][0]), nLev, nLat, nLon);

	varOutRho->set_cur(0, 0, 0);
	varOutRho->put(&(dRho[0][0][0]), nLev, nLat, nLon);

	ncdf_out.close();

	AnnounceEndBlock("Done");

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}
}

///////////////////////////////////////////////////////////////////////////////

