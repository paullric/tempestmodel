///////////////////////////////////////////////////////////////////////////////
///
///	\file    BaldaufGravityWaveTest.cpp
///	\author  Paul Ullrich
///	\version June 23, 2013
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

#include "Tempest.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Baldauf and Brdar (2013) inertia-gravity wave test (modified).
///	</summary>
class BaldaufGravityWaveTest : public TestCase {

protected:
	///	<summary>
	///		Background temperature.
	///	</summary>
	double ParamT0;

	///	<summary>
	///		Model cap.
	///	</summary>
	double ParamZtop;

	///	<summary>
	///		Potential temperature perturbation.
	///	</summary>
	double ParamPertMagnitude;

protected:
	///	<summary>
	///		Earth radius scaling parameter.
	///	</summary>
	double m_dEarthRadiusScaling;

	///	<summary>
	///		Model height cap.
	///	</summary>
	double m_dZtop;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	BaldaufGravityWaveTest(
		double dEarthRadiusScaling
	) :
		ParamT0(300.0),
		ParamZtop(10000.0),
		ParamPertMagnitude(1.0),

		m_dEarthRadiusScaling(dEarthRadiusScaling)
	{ }

public:
	///	<summary>
	///		Number of tracers used in this test.
	///	</summary>
	virtual int GetTracerCount() const {
		return 0;
	}

	///	<summary>
	///		Get the altitude of the model cap.
	///	</summary>
	virtual double GetZtop() const {
		return ParamZtop;
	}

	///	<summary>
	///		Flag indicating that a reference state is available.
	///	</summary>
	virtual bool HasReferenceState() const {
		return true;
	}

	///	<summary>
	///		Obtain test case specific physical constants.
	///	</summary>
	virtual void EvaluatePhysicalConstants(
		PhysicalConstants & phys
	) const {
		phys.SetOmega(0.0);
		phys.SetEarthRadius(phys.GetEarthRadius() / m_dEarthRadiusScaling);
	}

	///	<summary>
	///		Evaluate the topography at the given point.
	///	</summary>
	virtual double EvaluateTopography(
		const PhysicalConstants & phys,
		double dLon,
		double dLat
	) const {
		return (0.0);
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

		// Calculate the isothermal pressure
		double dPressure =
			phys.GetP0() * exp(- phys.GetG() * dZ / phys.GetR() / ParamT0);

		// Calculate exact density
		double dRho = dPressure / (phys.GetR() * ParamT0);

		// Store the state
		dState[0] = 0.0;
		dState[1] = 0.0;
		dState[2] = phys.RhoThetaFromPressure(dPressure) / dRho;
		dState[3] = 0.0;
		dState[4] = dRho;
	}

	///	<summary>
	///		Evaluate the state vector at the given point.
	///	</summary>
	virtual void EvaluatePointwiseState(
		const PhysicalConstants & phys,
		const Time & time,
		double dZ,
		double dLon,
		double dLat,
		double * dState,
		double * dTracer
	) const {

		// Calculate the isothermal pressure
		double dPressure =
			phys.GetP0() * exp(- phys.GetG() * dZ / (phys.GetR() * ParamT0));

		// Calculate exact density
		double dRho = dPressure / (phys.GetR() * ParamT0);

		// Calculate exact temperature
		double dTemperature = dPressure / (dRho * phys.GetR());

		dTemperature += ParamPertMagnitude
			* exp(-100.0 * dLat * dLat)
			* sin(M_PI * dZ / ParamZtop);

		// Recalculate Rho
		dRho = dPressure / (phys.GetR() * dTemperature);

/*
		// Add the perturbation
		double dTb = ParamPertMagnitude
			* exp(100.0 * (sin(dLat) - 1.0))
			* sin(M_PI * dZ / ParamZtop);

		double dRhob =
			phys.GetP0() / (phys.GetR() * ParamT0) * (-dTb / ParamT0);

		dRho += exp(- 0.5 * phys.GetG() * dZ / phys.GetR() / ParamT0) * dRhob;
*/
		// Calculate the potential temperature
		double dTheta = phys.RhoThetaFromPressure(dPressure) / dRho;
/*
		// Add the perturbation
		dTheta += ParamPertMagnitude
			* exp(-100.0 * dLat * dLat)
			* sin(M_PI * dZ / ParamZtop);
*/
		// Store the state
		dState[0] = 0.0;
		dState[1] = 0.0;
		dState[2] = dTheta;
		dState[3] = 0.0;
		dState[4] = dRho;
	}
};

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	// Initialize Tempest
	TempestInitialize(&argc, &argv);

try {
	// Earth radius scaling
	double dEarthRadiusScaling;

	// Parse the command line
	BeginTempestCommandLine("BaldaufGravityWaveTest");
		SetDefaultResolution(20);
		SetDefaultLevels(10);
		SetDefaultOutputDeltaT("200s");
		SetDefaultDeltaT("200s");
		SetDefaultEndTime("200s");
		SetDefaultHorizontalOrder(4);
		SetDefaultVerticalOrder(1);

		CommandLineDouble(dEarthRadiusScaling, "radiusscale", 1.0);

		ParseCommandLine(argc, argv);
	EndTempestCommandLine(argv)

	// Setup the Model
	AnnounceBanner("MODEL SETUP");

	Model model(EquationSet::PrimitiveNonhydrostaticEquations);

	TempestSetupCubedSphereModel(model);

	// Set the test case for the model
	AnnounceStartBlock("Initializing test case");
	model.SetTestCase(
		new BaldaufGravityWaveTest(dEarthRadiusScaling));
	AnnounceEndBlock("Done");

	// Begin execution
	AnnounceBanner("SIMULATION");
	model.Go();

	// Compute error norms
	AnnounceBanner("RESULTS");
	model.ComputeErrorNorms();
	AnnounceBanner();

} catch(Exception & e) {
	std::cout << e.ToString() << std::endl;
}

	// Deinitialize Tempest
	TempestDeinitialize();
}

///////////////////////////////////////////////////////////////////////////////

