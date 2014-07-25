///////////////////////////////////////////////////////////////////////////////
///
///	\file    StationaryMountainFlowTest.cpp
///	\author  Paul Ullrich
///	\version April 25, 2014
///
///	<remarks>
///		Copyright 2000-2014 Paul Ullrich
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
///		Pressure gradient errors induced by a Schar-type mountain
///		(DCMIP 2012 test 2-0-0)
///	</summary>
class StationaryMountainFlowTest : public TestCase {

protected:
	///	<summary>
	///		Model height cap.
	///	</summary>
	double m_dZtop;

	///	<summary>
	///		Model scaling parameter
	///	</summary>
	const double m_dEarthScaling;

	///	<summary>
	///		Rotation rate of the Earth with X = 1.
	///	</summary>
	double m_dOmega;

	///	<summary>
	///		Reference temperature.
	///	</summary>
	double m_dT0;

	///	<summary>
	///		Temperature lapse rate.
	///	</summary>
	double m_dGamma;

	///	<summary>
	///		Longitude of Schar-type mountain centerpoint.
	///	</summary>
	double m_dLonM;

	///	<summary>
	///		Latitude of Schar-type mountain centerpoint.
	///	</summary>
	double m_dLatM;

	///	<summary>
	///		Maximum Schar-type mountain height.
	///	</summary>
	double m_dH0;

	///	<summary>
	///		Schar-type mountain radius (radians).
	///	</summary>
	double m_dRM;

	///	<summary>
	///		Schar-type mountain oscillation half-width (radians).
	///	</summary>
	double m_dZetaM;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	StationaryMountainFlowTest(
		double dZtop,
		double dEarthScaling,
		double dOmega,
		double dT0,
		double dGamma,
		double dLonM,
		double dLatM,
		double dH0,
		double dRM,
		double dZetaM
	) :
		m_dZtop(dZtop),
		m_dEarthScaling(dEarthScaling),
		m_dOmega(dOmega),
		m_dT0(dT0),
		m_dGamma(dGamma),
		m_dLonM(dLonM * M_PI / 180.0),
		m_dLatM(dLatM * M_PI / 180.0),
		m_dH0(dH0),
		m_dRM(dRM * M_PI / 180.0),
		m_dZetaM(dZetaM * M_PI / 180.0)
	{
	}

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
		return m_dZtop;
	}

	///	<summary>
	///		Obtain test case specific physical constants.
	///	</summary>
	virtual void EvaluatePhysicalConstants(
		PhysicalConstants & phys
	) const {
		phys.SetOmega(m_dOmega * m_dEarthScaling);
		phys.SetEarthRadius(phys.GetEarthRadius() / m_dEarthScaling);
	}

	///	<summary>
	///		Evaluate the topography at the given point.
	///	</summary>
	virtual double EvaluateTopography(
		const PhysicalConstants & phys,
		double dLon,
		double dLat
	) const {

		// Great circle distance from mountain centerpoint (radians)
		double dR = acos(sin(m_dLatM) * sin(dLat)
				+ cos(m_dLatM) * cos(dLat) * cos(dLon - m_dLonM));

		// Topography height
		double dCosTerm = cos(M_PI * dR / m_dZetaM);

		double dBellTerm = 0.5 * (1.0 + cos(M_PI * dR / m_dRM));

		if (dR >= m_dRM) {
			dBellTerm = 0.0;
		}

		return (m_dH0 * dBellTerm * dCosTerm * dCosTerm);
	}

	///	<summary>
	///		Flag indicating whether or not Rayleigh friction strength is given.
	///	</summary>
	virtual bool HasRayleighFriction() const {
		return false;
	}

	///	<summary>
	///		Evaluate the Rayleigh friction strength at the given point.
	///	</summary>
	virtual double EvaluateRayleighStrength(
		double dZ,
		double dXp,
		double dYp
	) const {
		return (0.0);
	}

	///	<summary>
	///		Flag indicating that a reference state is available.
	///	</summary>
	virtual bool HasReferenceState() const {
		return true;
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

		// 3D temperature
		double dT = m_dT0 - m_dGamma * dZ;

		// 3D pressure
		double dPressure = phys.GetP0()
			* pow(1.0 - m_dGamma / m_dT0 * dZ,
				phys.GetG() / (phys.GetR() * m_dGamma));

		// 3D density
		double dRho = dPressure / (phys.GetR() * dT);

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
		return EvaluateReferenceState(phys, dZ, dLon, dLat, dState);
	}
};

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	// Initialize MPI
	TempestInitialize(&argc, &argv);

try {
	// Model height cap.
	double dZtop;

	// Model scaling parameter
	double dEarthScaling;

	// Rotation rate of the Earth with X = 1.
	double dOmega;

	// Reference temperature.
	double dT0;

	// Temperature lapse rate.
	double dGamma;

	// Longitude of Schar-type mountain centerpoint.
	double dLonM;

	// Latitude of Schar-type mountain centerpoint.
	double dLatM;

	// Maximum Schar-type mountain height.
	double dH0;

	// Schar-type mountain radius (radians).
	double dRM;

	// Schar-type mountain oscillation half-width (radians).
	double dZetaM;

	// Parse the command line
	BeginTempestCommandLine("StationaryMountainFlowTest");
		SetDefaultResolution(30);
		SetDefaultLevels(30);
		SetDefaultOutputDeltaT("1d");
		SetDefaultDeltaT("300s");
		SetDefaultEndTime("6d");
		SetDefaultHorizontalOrder(4);
		SetDefaultVerticalOrder(1);

		CommandLineDouble(dZtop, "ztop", 30000.0);
		CommandLineDouble(dEarthScaling, "X", 1.0);
		CommandLineDouble(dOmega, "omega", 0.0);
		CommandLineDouble(dT0, "T0", 300.0);
		CommandLineDouble(dGamma, "Gamma", 0.0065);
		CommandLineDouble(dLonM, "lonm", 270.0);
		CommandLineDouble(dLatM, "latm", 0.0);
		CommandLineDouble(dH0, "h0", 2000.0);
		CommandLineDouble(dRM, "rm", 135.0);
		CommandLineDouble(dZetaM, "zetam", 11.25);

		ParseCommandLine(argc, argv);
	EndTempestCommandLine(argv)

	// Setup the Model
	AnnounceBanner("MODEL SETUP");

	Model model(EquationSet::PrimitiveNonhydrostaticEquations);

	TempestSetupCubedSphereModel(model);

	// Set the test case for the model
	AnnounceStartBlock("Initializing test case");

	model.SetTestCase(
		new StationaryMountainFlowTest(
			dZtop,
			dEarthScaling,
			dOmega,
			dT0,
			dGamma,
			dLonM,
			dLatM,
			dH0,
			dRM,
			dZetaM));

	AnnounceEndBlock("Done");

	// Set the reference length
	model.GetGrid()->SetReferenceLength(0.5 * M_PI / 30.0 * dEarthScaling);

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

