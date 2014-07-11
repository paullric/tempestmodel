///////////////////////////////////////////////////////////////////////////////
///
///	\file    MountainRossby3DTest.cpp
///	\author  Paul Ullrich
///	\version June 24, 2014
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
///		Mountain-induced Rossby wave train (DCMIP 2008 Test 5)
///	</summary>
class MountainRossby3DTest : public TestCase {

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
	///		Longitude of mountain centerpoint.
	///	</summary>
	double m_dLonC;

	///	<summary>
	///		Latitude of mountain centerpoint.
	///	</summary>
	double m_dLatC;

	///	<summary>
	///		Maximum mountain height.
	///	</summary>
	double m_dH0;

	///	<summary>
	///		Mountain half-width.
	///	</summary>
	double m_dD;

	///	<summary>
	///		Pole pressure.
	///	</summary>
	double m_dPp;

	///	<summary>
	///		Isothermal atmosphere temperature.
	///	</summary>
	double m_dT0;

	///	<summary>
	///		Reference zonal wind velocity.
	///	</summary>
	double m_dU0;

	///	<summay>
	///		True if Rayleigh damping is disabled.
	///	</summary>
	bool m_fNoRayleigh;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	MountainRossby3DTest(
		double dZtop,
		double dEarthScaling,
		double dLonC,
		double dLatC,
		double dH0,
		double dD,
		double dPp,
		double dT0,
		double dU0,
		bool fNoRayleigh
	) :
		m_dZtop(dZtop),
		m_dEarthScaling(dEarthScaling),
		m_dLonC(dLonC * M_PI / 180.0),
		m_dLatC(dLatC * M_PI / 180.0),
		m_dH0(dH0),
		m_dD(dD),
		m_dPp(dPp),
		m_dT0(dT0),
		m_dU0(dU0),
		m_fNoRayleigh(fNoRayleigh)
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
		phys.SetOmega(phys.GetOmega() * m_dEarthScaling);
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

		// Great circle distance from mountain centerpoint
		double dR = phys.GetEarthRadius() * acos(
			sin(m_dLatC) * sin(dLat)
				+ cos(m_dLatC) * cos(dLat) * cos(dLon - m_dLonC));

		// Topography height
		double dExpTerm = exp(- dR * dR / (m_dD * m_dD));

		return (m_dH0 * dExpTerm);
	}

	///	<summary>
	///		Flag indicating whether or not Rayleigh friction strength is given.
	///	</summary>
	virtual bool HasRayleighFriction() const {
		return !m_fNoRayleigh;
	}

	///	<summary>
	///		Evaluate the Rayleigh friction strength at the given point.
	///	</summary>
	virtual double EvaluateRayleighStrength(
		double dZ,
		double dLon,
		double dLat
	) const {
        const double dRayleighStrength = 4.0e-3;
        const double dRayleighDepth = 10000.0;

        double dNuDepth = 0.0;

        if (dZ > m_dZtop - dRayleighDepth) {
            double dNormZ = (m_dZtop - dZ) / dRayleighDepth;
            dNuDepth = 0.5 * dRayleighStrength * (1.0 + cos(M_PI * dNormZ));
        }

        return dNuDepth;
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

		double dSin2Lat = sin(dLat) * sin(dLat);

		// 3D temperature
		double dT = m_dT0;

		// 3D pressure
		double dPressure =
			m_dPp * exp(
				- m_dU0 / (2.0 * phys.GetR() * m_dT0) * (dSin2Lat - 1.0)
					* (m_dU0 + 2.0 * phys.GetOmega() * phys.GetEarthRadius())
				- phys.GetG() * dZ / (phys.GetR() * m_dT0));

		// 3D density
		double dRho = dPressure / (phys.GetR() * m_dT0);

		// Store the state
		dState[0] = m_dU0 * cos(dLat);
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

	// Longitude of Schar-type mountain centerpoint.
	double dLonC;

	// Latitude of Schar-type mountain centerpoint.
	double dLatC;

	// Maximum Schar-type mountain height.
	double dH0;

	// Schar-type mountain half width.
	double dD;

	// Polar pressure
	double dPp;

	// Isothermal atmosphere temperature
	double dT0;

	// Reference zonal wind velocity.
	double dU0;

	// No Rayleigh damping
	bool fNoRayleigh;

	// Parse the command line
	BeginTempestCommandLine("MountainRossby3DTest");
		SetDefaultResolution(30);
		SetDefaultLevels(30);
		SetDefaultOutputDeltaT("86400s");
		SetDefaultDeltaT("200s");
		SetDefaultEndTime("30d");
		SetDefaultHorizontalOrder(4);
		SetDefaultVerticalOrder(1);

		CommandLineDouble(dZtop, "ztop", 30000.0);
		CommandLineDouble(dEarthScaling, "X", 1.0);
		CommandLineDouble(dLonC, "lonc", 90.0);
		CommandLineDouble(dLatC, "latc", 30.0);
		CommandLineDouble(dH0, "h0", 2000.0);
		CommandLineDouble(dD, "d", 1.5e6);
		CommandLineDouble(dPp, "pp", 93000.0);
		CommandLineDouble(dT0, "t0", 288.0);
		CommandLineDouble(dU0, "u0", 20.0);
		CommandLineBool(fNoRayleigh, "norayleigh");

		ParseCommandLine(argc, argv);
	EndTempestCommandLine(argv)

	// Setup the Model
	AnnounceBanner("MODEL SETUP");

	Model model(EquationSet::PrimitiveNonhydrostaticEquations);

	TempestSetupCubedSphereModel(model);

	// Set the test case for the model
	AnnounceStartBlock("Initializing test case");

	model.SetTestCase(
		new MountainRossby3DTest(
			dZtop,
			dEarthScaling,
			dLonC,
			dLatC,
			dH0,
			dD,
			dPp,
			dT0,
			dU0,
			fNoRayleigh));

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

