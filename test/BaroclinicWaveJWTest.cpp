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
#include "VerticalDynamicsStub.h"

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

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	BaroclinicWaveJWTest(
		double dAlpha,
		bool fTracerOn
	) :
		m_dAlpha(dAlpha),
		m_fTracerOn(fTracerOn)
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
		return 30000.0;
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
	) {

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
/*
		if (fabs(dLat - 0.5 * M_PI) < 1.0e-12) {
			dLat -= 1.0e-12;
		}
		if (fabs(dLat + 0.5 * M_PI) < 1.0e-12) {
			dLat += 1.0e-12;
		}

		// Height field
		double dHTrig =
			- cos(dLon) * cos(dLat) * sin(m_dAlpha)
			+ sin(dLat) * cos(m_dAlpha);

		dState[0] = m_dH0
			- (phys.GetEarthRadius() * phys.GetOmega() + 0.5 * m_dU0)
				* m_dU0 * dHTrig * dHTrig / phys.GetG();

		// Pointwise zonal component of velocity
		dState[1] = m_dU0 * cos(dLat)
			* (cos(m_dAlpha) + cos(dLon) * tan(dLat) * sin(m_dAlpha));

		// Pointwise meridional component of velocity
		dState[2] = - m_dU0 * sin(dLon) * sin(m_dAlpha);

		// Cosine bell tracer field
		if (m_fTracerOn) {
			const double ParamLonC = 1.5 * M_PI;
			const double ParamLatC = 0.0;
			const double ParamH0 = 1000.0;
			const double ParamR = phys.GetEarthRadius() / 3.0;

			double dR = phys.GetEarthRadius()
				* acos(sin(ParamLatC) * sin(dLat)
					+ cos(ParamLatC) * cos(dLat) * cos(dLon - ParamLonC));

			if (dR < ParamR) {
				dTracer[0] = (ParamH0 / 2.0) * (1.0 + cos(M_PI * dR / ParamR));
			} else {
				dTracer[0] = 0.0;
			}
		}
*/
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

	// Order
	int nOrder;

	// Model height cap
	double dZtop;

	// Grid rotation angle
	double dAlpha;

	// Include tracer field
	bool fTracersOn;

	// Output time
	double dOutputDeltaT;

	// Numerical method
	std::string strHorizontalDynamics;

	// Use hyperdiffusion
	bool fUseHyperviscosity;

	// Model parameters
	ModelParameters params;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strOutputDir, "output_dir", "out");
		CommandLineString(strOutputPrefix, "output_prefix", "out");
		CommandLineInt(nResolution, "resolution", 40);
		CommandLineInt(nOrder, "order", 4);
		CommandLineDouble(dZtop, "ztop", 30000.0);
		CommandLineDouble(dAlpha, "alpha", 0.0);
		CommandLineBool(fTracersOn, "with_tracer");
		CommandLineDouble(params.m_dDeltaT, "dt", 200.0);
		CommandLineDouble(params.m_dEndTime, "endtime", 200.0);//86400.0 * 5.0);
		CommandLineDouble(dOutputDeltaT, "outputtime", 86400.0);
		CommandLineStringD(strHorizontalDynamics, "method", "SE", "(SE | DG)");
		CommandLineBool(fUseHyperviscosity, "hypervis");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner("INITIALIZATION");

	// Construct a model
	AnnounceStartBlock("Creating model");
	Model model(EquationSet::PrimitiveNonhydrostaticEquations);
	AnnounceEndBlock("Done");

	// Generate a new cubed-sphere GLL grid
	AnnounceStartBlock("Creating grid");
	GridCSGLL grid(
		model,
		nResolution,
		4,
		nOrder,
		nOrder,
		1);

	AnnounceEndBlock("Done");

	// Set the parameters for the model
	AnnounceStartBlock("Initializing parameters");
	model.SetParameters(&params);
	AnnounceEndBlock("Done");

	// Set the grid for the model
	AnnounceStartBlock("Initializing grid");
	model.SetGrid(&grid);
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
		model, nOrder, eHorizontalDynamicsType, fUseHyperviscosity);
	AnnounceStartBlock("Initializing horizontal dynamics");
	model.SetHorizontalDynamics(&hdyn);
	AnnounceEndBlock("Done");

	// Set the vertical dynamics
	VerticalDynamicsStub vdyn(model);
	AnnounceStartBlock("Initializing vertical dynamics");
	model.SetVerticalDynamics(&vdyn);
	AnnounceEndBlock("Done");

	// Set the test case for the model
	BaroclinicWaveJWTest test(dAlpha, fTracersOn);
	AnnounceStartBlock("Initializing data");
	model.SetTestCase(&test);
	AnnounceEndBlock("Done");

	// Set the reference output manager for the model
	AnnounceStartBlock("Creating reference output manager");
	OutputManagerReference outmanRef(
		grid, dOutputDeltaT, strOutputDir, strOutputPrefix,
		720, 360, false, false);
	model.AttachOutputManager(&outmanRef);

	outmanRef.InitializeNcOutput("ref.nc");
	AnnounceEndBlock("Done");

	// Set the composite output manager for the model
	AnnounceStartBlock("Creating composite output manager");
	OutputManagerComposite outmanComp(
		grid, dOutputDeltaT, strOutputDir, strOutputPrefix);
	model.AttachOutputManager(&outmanComp);

	outmanComp.InitializeNcOutput("comp.nc");
	AnnounceEndBlock("Done");

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

