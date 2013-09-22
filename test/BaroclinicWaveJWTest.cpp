///////////////////////////////////////////////////////////////////////////////
///
///	\file    BaroclinicWaveTest.cpp
///	\author  Paul Ullrich
///	\version May 23, 2013
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

#include "CommandLine.h"
#include "Announce.h"
#include "STLStringHelper.h"

#include "Model.h"
#include "TimestepSchemeARK4.h"
#include "HorizontalDynamicsFEM.h"
#include "VerticalDynamicsFEM.h"

#include "PhysicalConstants.h"
#include "TestCase.h"
#include "OutputManagerComposite.h"
#include "OutputManagerReference.h"
#include "OutputManagerChecksum.h"
#include "GridData4D.h"
#include "EquationSet.h"

#include "GridCSGLL.h"

#include "mpi.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Jablonowski and Williamson (2006) Baroclinic wave test
///	</summary>
class BaroclinicWaveJWTest : public TestCase {

public:
	///	<summary>
	///		Auxiliary eta.
	///	</summary>
	static const double ParamEta0 = 0.252;

	///	<summary>
	///		Tropopause level (in eta coordinates)
	///	</summary>
	static const double ParamTropopauseEta = 0.2;

	///	<summary>
	///		Horizontal-mean temperature (K)
	///	</summary>
	static const double ParamT0 = 288.0;

	///	<summary>
	///		Empirical temperature difference (K)
	///	</summary>
	static const double ParamDeltaT = 4.8e5;

	///	<summary>
	///		Temperature lapse rate (K/m)
	///	</summary>
	static const double ParamLapseRate = 0.005;

	///	<summary>
	///		Maximum zonal wind (m/s)
	///	</summary>
	static const double ParamU0 = 35.0;

protected:
	///	<summary>
	///		Alpha parameter.
	///	</summary>
	double m_dAlpha;

	///	<summary>
	///		Whether to use an auxilliary tracer field.
	///	</summary>
	bool m_fTracerOn;

	///	<summary>
	///		Model height cap.
	///	</summary>
	double m_dZtop;

	///	<summary>
	///		Flag indicating that the reference profile should be used.
	///	</summary>
	bool m_fNoReferenceState;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	BaroclinicWaveJWTest(
		double dAlpha,
		bool fTracerOn,
		double dZtop,
		bool fNoReferenceState
	) :
		m_dAlpha(dAlpha),
		m_fTracerOn(fTracerOn),
		m_dZtop(dZtop),
		m_fNoReferenceState(fNoReferenceState)
	{ }

public:
	///	<summary>
	///		Number of tracers used in this test.
	///	</summary>
	virtual int GetTracerCount() const {
		return (m_fTracerOn)?(1):(0);
	}

	///	<summary>
	///		Get the altitude of the model cap.
	///	</summary>
	virtual double GetZtop() const {
		return m_dZtop;
	}

	///	<summary>
	///		Flag indicating that a reference state is available.
	///	</summary>
	virtual bool HasReferenceState() const {
		return !m_fNoReferenceState;
	}

	///	<summary>
	///		Obtain test case specific physical constants.
	///	</summary>
	virtual void EvaluatePhysicalConstants(
		PhysicalConstants & phys
	) const {
		// Set the alpha parameter
		phys.SetAlpha(m_dAlpha);
	}

	///	<summary>
	///		Evaluate the topography at the given point.
	///	</summary>
	virtual double EvaluateTopography(
		const PhysicalConstants & phys,
		double dLon,
		double dLat
	) const {

		// Calculate auxiliary eta
		double dAuxEta = 0.5 * M_PI * (1.0 - ParamEta0);

		// Various powers of trigonometric functions of latitude
		double dSinLat = sin(dLat);
		double dSinLat2 = dSinLat * dSinLat;
		double dSinLat3 = dSinLat * dSinLat2;
		double dSinLat4 = dSinLat * dSinLat3;
		double dSinLat5 = dSinLat * dSinLat4;
		double dSinLat6 = dSinLat * dSinLat5;

		double dCosLat = cos(dLat);
		double dCosLat2 = dCosLat * dCosLat;
		double dCosLat3 = dCosLat * dCosLat2;

		// Reference profiles
		double dRefProfile1 =
			ParamU0 * pow(cos(dAuxEta), 1.5)
				* (-2.0 * dSinLat6 * (dCosLat2 + 1.0 / 3.0) + 10.0 / 63.0);

		double dRefProfile2 =
			phys.GetEarthRadius() * phys.GetOmega()
				* (8.0 / 5.0 * dCosLat3 * (dSinLat2 + 2.0 / 3.0) - 0.25 * M_PI);

		double dSurfGeopotential =
			ParamU0 * pow(cos(dAuxEta), 1.5) *
				(dRefProfile1 + dRefProfile2);

		return (dSurfGeopotential / phys.GetG());
	}

	///	<summary>
	///		Calculate the geopotential and temperature at the given point.
	///	</summary>
	void CalculateGeopotentialTemperature(
		const PhysicalConstants & phys,
		double dEta,
		double dLon,
		double dLat,
		double & dGeopotential,
		double & dTemperature
	) const {

		// Calculate auxiliary eta
		double dAuxEta = 0.5 * M_PI * (dEta - ParamEta0);

		// Horizontally averaged temperature
		double dAvgTemperature =
			ParamT0 * pow(dEta, phys.GetR() * ParamLapseRate / phys.GetG());

		if (dEta < ParamTropopauseEta) {
			dAvgTemperature +=
				ParamDeltaT * pow(ParamTropopauseEta - dEta, 5.0);
		}

		// Various powers of trigonometric functions of latitude
		double dSinLat = sin(dLat);
		double dSinLat2 = dSinLat * dSinLat;
		double dSinLat3 = dSinLat * dSinLat2;
		double dSinLat4 = dSinLat * dSinLat3;
		double dSinLat5 = dSinLat * dSinLat4;
		double dSinLat6 = dSinLat * dSinLat5;

		double dCosLat = cos(dLat);
		double dCosLat2 = dCosLat * dCosLat;
		double dCosLat3 = dCosLat * dCosLat2;

		// Reference profiles
		double dRefProfile1 =
			ParamU0 * pow(cos(dAuxEta), 1.5)
				* (-2.0 * dSinLat6 * (dCosLat2 + 1.0 / 3.0) + 10.0 / 63.0);

		double dRefProfile2 =
			phys.GetEarthRadius() * phys.GetOmega()
				* (8.0 / 5.0 * dCosLat3 * (dSinLat2 + 2.0 / 3.0) - 0.25 * M_PI);

		// Total temperature distribution
		dTemperature = 2.0 * dRefProfile1 + dRefProfile2;

		dTemperature =
			dAvgTemperature
			+ 0.75 * dEta * M_PI * ParamU0 / phys.GetR()
				* sin(dAuxEta) * sqrt(cos(dAuxEta)) * dTemperature;

		// Geopotential distribution
		double dAvgGeopotential =
			ParamT0 * phys.GetG() / ParamLapseRate
				* (1.0 - pow(dEta, phys.GetR() * ParamLapseRate / phys.GetG()));

		if (dEta < ParamTropopauseEta) {
			double dEta2 = dEta * dEta;
			double dEta3 = dEta * dEta2;
			double dEta4 = dEta * dEta3;
			double dEta5 = dEta * dEta4;

			double dTropoEta  = ParamTropopauseEta;
			double dTropoEta2 = dTropoEta * dTropoEta;
			double dTropoEta3 = dTropoEta * dTropoEta2;
			double dTropoEta4 = dTropoEta * dTropoEta3;
			double dTropoEta5 = dTropoEta * dTropoEta4;

			dAvgGeopotential -= phys.GetR() * ParamDeltaT * (
				(log(dEta / ParamTropopauseEta) + 137.0 / 60.0) * dTropoEta5
				- 5.0 * dTropoEta4 * dEta
				+ 5.0 * dTropoEta3 * dEta2
				- (10.0 / 3.0) * dTropoEta2 * dEta3
				+ 5.0 / 4.0 * dTropoEta * dEta4
				- 1.0 / 5.0 * dEta5);
		}

		dGeopotential = dAvgGeopotential
			+ ParamU0 * pow(cos(dAuxEta), 1.5)
				* (dRefProfile1 + dRefProfile2);

	}

	///	<summary>
	///		Calculate eta at the given point via Newton iteration.  The
	///		geopotential and temperature at this point are also returned via
	///		command-line parameters.
	///	</summary>
	double EtaFromRLL(
		const PhysicalConstants &phys,
		double dZ,
		double dLon,
		double dLat,
		double & dGeopotential,
		double & dTemperature
	) const {
		const int MaxIterations  = 25;
		const double InitialEta  = 1.0e-7;
		const double Convergence = 1.0e-13;

		// Buffer variables
		double dEta = InitialEta;
		double dNewEta;

		double dF;
		double dDiffF;

		// Iterate until convergence is achieved
		int i = 0;
		for (; i < MaxIterations; i++) {

			CalculateGeopotentialTemperature(
				phys, dEta, dLon, dLat, dGeopotential, dTemperature);

			dF     = - phys.GetG() * dZ + dGeopotential;
			dDiffF = - phys.GetR() / dEta * dTemperature;

			dNewEta = dEta - dF / dDiffF;

			if (fabs(dEta - dNewEta) < Convergence) {
				return dNewEta;
			}

			dEta = dNewEta;
		}

		// Check for convergence failure
		if (i == MaxIterations) {
			_EXCEPTIONT("Maximum number of iterations exceeded.");
		}

		if ((dEta > 1.0) || (dEta < 0.0)) {
			_EXCEPTIONT("Invalid Eta value");
		}

		return dEta;
	}


	///	<summary>
	///		Evaluate the reference state at the given point.
	///	</summary>
	virtual void EvaluateReferenceState(
		const PhysicalConstants & phys,
		double dZ,
		double dLon,
		double dLat,
		double * dState
	) const {

		// Pressure coordinate
		double dGeopotential;
		double dTemperature;

		double dEta = EtaFromRLL(
			phys, dZ, dLon, dLat, dGeopotential, dTemperature);


		// Calculate zonal velocity
		double dUlon =
			ParamU0 * pow(cos(0.5 * M_PI * (dEta - ParamEta0)), 1.5)
				* sin(2.0 * dLat) * sin(2.0 * dLat);

		dState[0] = dUlon;

		// Calculate rho and theta
		double dPressure = phys.GetP0() * dEta;

		double dRho = dPressure / (phys.GetR() * dTemperature);

		double dRhoTheta = phys.RhoThetaFromPressure(dPressure);

		dState[2] = dRhoTheta / dRho;
		dState[4] = dRho;
	}

	///	<summary>
	///		Evaluate the state vector at the given point.
	///	</summary>
	virtual void EvaluatePointwiseState(
		const PhysicalConstants & phys,
		double dTime,
		double dZ,
		double dLon,
		double dLat,
		double * dState,
		double * dTracer
	) const {

		// Evaluate the reference state at this point
		EvaluateReferenceState(phys, dZ, dLon, dLat, dState);
	}

};

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	// Initialize MPI
	MPI_Init(&argc, &argv);

try {
	// Output directory
	std::string strOutputDir;

	// Output file prefix
	std::string strOutputPrefix;

	// Resolution
	int nResolution;

	// Number of vertical levels
	int nLevels;

	// Order
	int nHorizontalOrder;

	// Vertical order
	int nVerticalOrder;

	// Model height cap
	double dZtop;

	// Grid rotation angle
	double dAlpha;

	// Use reference state flag
	bool fNoReferenceState;

	// Include tracer field
	bool fTracersOn;

	// Output time
	double dOutputDeltaT;

	// Numerical method
	std::string strHorizontalDynamics;

	// Use hyperdiffusion
	bool fNoHyperviscosity;

	// Model parameters
	ModelParameters params;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strOutputDir, "output_dir", "out");
		CommandLineString(strOutputPrefix, "output_prefix", "out");
		CommandLineInt(nResolution, "resolution", 20);
		CommandLineInt(nLevels, "levels", 10);
		CommandLineInt(nHorizontalOrder, "order", 4);
		CommandLineInt(nVerticalOrder, "vertorder", 1);
		CommandLineDouble(dZtop, "ztop", 10000.0);
		CommandLineDouble(dAlpha, "alpha", 0.0);
		CommandLineBool(fNoReferenceState, "norefstate");
		CommandLineBool(fTracersOn, "with_tracer");
		CommandLineDouble(params.m_dDeltaT, "dt", 200.0);
		CommandLineDouble(params.m_dEndTime, "endtime", 200.0);
		CommandLineDouble(dOutputDeltaT, "outputtime", 21600.0);
		CommandLineStringD(strHorizontalDynamics, "method", "SE", "(SE | DG)");
		CommandLineBool(fNoHyperviscosity, "nohypervis");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner("INITIALIZATION");

	// Construct a model
	AnnounceStartBlock("Creating model");
	Model model(EquationSet::PrimitiveNonhydrostaticEquations);
	AnnounceEndBlock("Done");

	// Set the parameters for the model
	AnnounceStartBlock("Initializing parameters");
	model.SetParameters(&params);
	AnnounceEndBlock("Done");

	// Set the timestep scheme
	TimestepSchemeARK4 timestep(model);
	AnnounceStartBlock("Initializing timestep scheme");
	model.SetTimestepScheme(&timestep);
	AnnounceEndBlock("Done");

	// Set the horizontal dynamics
	HorizontalDynamicsFEM::Type eHorizontalDynamicsType;
	STLStringHelper::ToLower(strHorizontalDynamics);
	if (strHorizontalDynamics == "se") {
		eHorizontalDynamicsType = HorizontalDynamicsFEM::SpectralElement;
	} else if (strHorizontalDynamics == "dg") {
		eHorizontalDynamicsType = HorizontalDynamicsFEM::DiscontinuousGalerkin;
	} else {
		_EXCEPTIONT("Invalid method: Expected \"SE\" or \"DG\"");
	}

	HorizontalDynamicsFEM hdyn(
		model, nHorizontalOrder, eHorizontalDynamicsType, fNoHyperviscosity);
	AnnounceStartBlock("Initializing horizontal dynamics");
	model.SetHorizontalDynamics(&hdyn);
	AnnounceEndBlock("Done");

	// Set the vertical dynamics
	VerticalDynamicsFEM vdyn(model, nHorizontalOrder, nVerticalOrder);
	AnnounceStartBlock("Initializing vertical dynamics");
	model.SetVerticalDynamics(&vdyn);
	AnnounceEndBlock("Done");

	// Generate a new cubed-sphere GLL grid
	AnnounceStartBlock("Initializing grid");
	GridCSGLL grid(
		model,
		nResolution,
		4,
		nHorizontalOrder,
		nVerticalOrder,
		nLevels);

	model.SetGrid(&grid);
	AnnounceEndBlock("Done");

	// Set the test case for the model
	BaroclinicWaveJWTest test(
		dAlpha,
		fTracersOn,
		dZtop,
		fNoReferenceState);

	AnnounceStartBlock("Initializing data");
	model.SetTestCase(&test);
	AnnounceEndBlock("Done");

	// Set the reference output manager for the model
	AnnounceStartBlock("Creating reference output manager");
	OutputManagerReference outmanRef(
		grid, dOutputDeltaT, strOutputDir, strOutputPrefix,
		360, 180, false, false, true, true);
	model.AttachOutputManager(&outmanRef);

	outmanRef.InitializeNcOutput("ref.nc");
	AnnounceEndBlock("Done");
/*
	// Set the composite output manager for the model
	AnnounceStartBlock("Creating composite output manager");
	OutputManagerComposite outmanComp(
		grid, dOutputDeltaT, strOutputDir, strOutputPrefix);
	model.AttachOutputManager(&outmanComp);

	outmanComp.InitializeNcOutput("comp.nc");
	AnnounceEndBlock("Done");
*/
	// Set the checksum output manager for the model
	AnnounceStartBlock("Creating checksum output manager");
	OutputManagerChecksum outmanChecksum(grid, dOutputDeltaT);
	model.AttachOutputManager(&outmanChecksum);
	AnnounceEndBlock("Done");

	// Begin execution
	AnnounceBanner("SIMULATION");
	model.Go();

	// Execution complete
	//AnnounceBanner("Execution completed successfully");

	// Compute error norms
	AnnounceBanner("RESULTS");
	model.ComputeErrorNorms();
	AnnounceBanner();

} catch(Exception & e) {
	std::cout << e.ToString() << std::endl;
}

	// Finalize MPI
	MPI_Finalize();
}

///////////////////////////////////////////////////////////////////////////////

