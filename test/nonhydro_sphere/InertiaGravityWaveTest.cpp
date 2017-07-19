///////////////////////////////////////////////////////////////////////////////
///
///	\file    InertiaGravityWaveTest.cpp
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
///		DCMIP 2012 Test 3-0-0.  Non-hydrostatic inertia gravity waves.
///	</summary>
class InertiaGravityWaveTest : public TestCase {

protected:
	///	<summary>
	///		Model cap.
	///	</summary>
	double m_dZtop;

	///	<summary>
	///		Earth radius scaling parameter.
	///	</summary>
	double m_dEarthScaling;

	///	<summary>
	///		Earth rotation rate parameter.
	///	</summary>
	double m_dOmega;

	///	<summary>
	///		Background wind speed.
	///	</summary>
	double m_dU0;

	///	<summary>
	///		Background Brunt-Vaisala frequency.
	///	</summary>
	double m_dN;

	///	<summary>
	///		Surface temperature at the equator.
	///	</summary>
	double m_dTeq;

	///	<summary>
	///		Potential temperature perturbation width parameter (m).
	///	</summary>
	double m_dPertWidth;

	///	<summary>
	///		Longitudinal centerpoint of the potential temperature pert.
	///	</summary>
	double m_dPertLonC;

	///	<summary>
	///		Latitudinal centerpoint of the potential temperature pert.
	///	</summary>
	double m_dPertLatC;

	///	<summary>
	///		Magnitude of the potential temperature perturbation.
	///	</summary>
	double m_dPertMagnitude;

	///	<summary>
	///		Vertical wavelength of the potential temperature perturbation.
	///	</summary>
	double m_dPertVerticalWavelength;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	InertiaGravityWaveTest(
		double dZtop,
		double dEarthScaling,
		double dOmega,
		double dU0,
		double dN,
		double dTeq,
		double dPertWidth,
		double dPertLonC,
		double dPertLatC,
		double dPertMagnitude,
		double dPertVerticalWavelength
	) :
		m_dZtop(dZtop),
		m_dEarthScaling(dEarthScaling),
		m_dOmega(dOmega),
		m_dU0(dU0),
		m_dN(dN),
		m_dTeq(dTeq),
		m_dPertWidth(dPertWidth),
		m_dPertLonC(dPertLonC * M_PI / 180.0),
		m_dPertLatC(dPertLatC * M_PI / 180.0),
		m_dPertMagnitude(dPertMagnitude),
		m_dPertVerticalWavelength(dPertVerticalWavelength)
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
		return m_dZtop;
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
		return (0.0);
	}

	///	<summary>
	///		Evaluate the reference state at the given point.
	///	</summary>
	virtual void EvaluateReferenceState(
		const PhysicalConstants & phys,
		double dXi,
		double dZ,
		double dLon,
		double dLat,
		double * dState
	) const {

		// Reference temperature
		double dG = phys.GetG() * phys.GetG() / (m_dN * m_dN * phys.GetCp());

		// Surface temperature
		double dTsExpTerm =
			- m_dU0 * m_dN * m_dN / 4.0 / (phys.GetG() * phys.GetG())
				* (m_dU0 + 2.0 * phys.GetOmega() * phys.GetEarthRadius())
				* (cos(2.0 * dLat) - 1.0);

		double dTs = dG + (m_dTeq - dG) * exp(dTsExpTerm);

		// 3D temperature
		double dT = dG + (dTs - dG) * exp(m_dN * m_dN * dZ / phys.GetG());

		// Surface pressure
		double dPsTempScaling = pow(dTs / m_dTeq, 1.0 / phys.GetKappa());

		double dPsExpTerm =
			m_dU0 / (4.0 * dG * phys.GetR())
				* (m_dU0 + 2.0 * phys.GetOmega() * phys.GetEarthRadius())
				* (cos(2.0 * dLat) - 1.0);

		double dPs = phys.GetP0() * exp(dPsExpTerm) * dPsTempScaling;

		// 3D Pressure
		double dPVertTerm =
			dG / dTs * exp(- m_dN * m_dN * dZ / phys.GetG()) + 1.0 - dG / dTs;

		double dPressure = dPs * pow(dPVertTerm, 1.0 / phys.GetKappa());

		// Check for negative pressure
		if (dPressure < 0.0) {
			_EXCEPTIONT("Negative pressure detected");
		}

		// Calculate exact density
		double dRho = dPressure / (phys.GetR() * dT);

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
		double dXi,
		double dZ,
		double dLon,
		double dLat,
		double * dState,
		double * dTracer
	) const {

		// Calculate the reference state
		EvaluateReferenceState(phys, dXi, dZ, dLon, dLat, dState);

		// Add in the potential temperature perturbation
		double dR = phys.GetEarthRadius() * acos(
			sin(m_dPertLatC) * sin(dLat)
			+ cos(m_dPertLatC) * cos(dLat) * cos(dLon - m_dPertLonC));

		double dS = m_dPertWidth * m_dPertWidth /
			(m_dPertWidth * m_dPertWidth + dR * dR);

		dState[2] += m_dPertMagnitude
			* dS * sin(2.0 * M_PI * dZ / m_dPertVerticalWavelength);
	}
};

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	// Initialize Tempest
	TempestInitialize(&argc, &argv);

try {
	// Model cap.
	double dZtop;

	// Earth radius scaling parameter.
	double dEarthScaling;

	// Earth rotation rate parameter.
	double dOmega;

	// Background wind speed.
	double dU0;

	// Background Brunt-Vaisala frequency.
	double dN;

	// Surface temperature at the equator.
	double dTeq;

	// Potential temperature perturbation width parameter (m).
	double dPertWidth;

	// Longitudinal centerpoint of the potential temperature pert.
	double dPertLonC;

	// Latitudinal centerpoint of the potential temperature pert.
	double dPertLatC;

	// Magnitude of the potential temperature perturbation.
	double dPertMagnitude;

	// Vertical wavelength of the potential temperature perturbation.
	double dPertVerticalWavelength;

	// Parse the command line
	BeginTempestCommandLine("InertiaGravityWaveTest");
		SetDefaultResolution(20);
		SetDefaultLevels(10);
		SetDefaultOutputDeltaT("1500000u");
		SetDefaultDeltaT("1500000u");
		SetDefaultEndTime("1500000u");
		SetDefaultHorizontalOrder(4);
		SetDefaultVerticalOrder(1);

		CommandLineDouble(dZtop, "ztop", 10000.0);
		CommandLineDouble(dEarthScaling, "X", 125.0);
		CommandLineDouble(dOmega, "omega", 0.0);
		CommandLineDouble(dU0, "u0", 20.0);
		CommandLineDouble(dN, "N", 0.01);
		CommandLineDouble(dTeq, "Teq", 300.0);
		CommandLineDouble(dPertWidth, "d", 5000.0);
		CommandLineDouble(dPertLonC, "lon_c", 120.0);
		CommandLineDouble(dPertLatC, "lat_c", 0.0);
		CommandLineDouble(dPertMagnitude, "dtheta", 1.0);
		CommandLineDouble(dPertVerticalWavelength, "Lz", 20000.0);

		ParseCommandLine(argc, argv);
	EndTempestCommandLine(argv)

	// Setup the Model
	AnnounceBanner("MODEL SETUP");

	Model model(EquationSet::PrimitiveNonhydrostaticEquations);

	TempestSetupCubedSphereModel(model);

	// Set the test case for the model
	AnnounceStartBlock("Initializing test case");
	model.SetTestCase(
		new InertiaGravityWaveTest(
			dZtop,
			dEarthScaling,
			dOmega,
			dU0,
			dN,
			dTeq,
			dPertWidth,
			dPertLonC,
			dPertLatC,
			dPertMagnitude,
			dPertVerticalWavelength));
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
