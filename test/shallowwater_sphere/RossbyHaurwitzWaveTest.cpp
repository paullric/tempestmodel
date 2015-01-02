///////////////////////////////////////////////////////////////////////////////
///
///	\file    RossbyHaurwitzWave.cpp
///	\author  Paul Ullrich
///	\version October 7, 2014
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
///		Williamson et al. (1994) Test case 6
///
///		Rossby-Haurwitz wave.
///	</summary>
class RossbyHaurwitzWaveTest : public TestCase {

private:
	///	<summary>
	///		Background height field.
	///	</summary>
	double m_dH0;

	///	<summary>
	///		Wavenumber of Rossby Haurwitz wave.
	///	</summary>
	double m_dR;

	///	<summary>
	///		W parameter.
	///	</summary>
	double m_dW;

	///	<summary>
	///		K parameter.
	///	</summary>
	double m_dK;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	RossbyHaurwitzWaveTest(
		double dH0,
		double dR,
		double dW,
		double dK
	) :
		m_dH0(dH0),
		m_dR(dR),
		m_dW(dW),
		m_dK(dK)
	{ }

public:
	///	<summary>
	///		Number of tracers used in this test.
	///	</summary>
	int GetTracerCount() const {
		return 0;
	}

	///	<summary>
	///		Obtain test case specific physical constants.
	///	</summary>
	virtual void EvaluatePhysicalConstants(
		PhysicalConstants & phys
	) const {
	}

	///	<summary>
	///		Evaluate the topography at the given point.
	///	</summary>
	virtual double EvaluateTopography(
		const PhysicalConstants & phys,
		double dLon,
		double dLat
	) const {
		return 0.0;
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
		if (fabs(dLat - 0.5 * M_PI) < 1.0e-12) {
			dLat -= 1.0e-12;
		}
		if (fabs(dLat + 0.5 * M_PI) < 1.0e-12) {
			dLat += 1.0e-12;
		}

		// Rotation rate
		double dOmega = phys.GetOmega();

		// Gravity
		double dG = phys.GetG();

		// Earth radius
		double dEarthRadius = phys.GetEarthRadius();

		// Calculate coefficients
		double dCoeffA =
			0.5 * m_dW * (2.0 * dOmega + m_dW) * cos(dLat) * cos(dLat)
			+ 0.25 * m_dK * m_dK * pow(cos(dLat), 2.0 * m_dR) *
				((m_dR + 1.0) * cos(dLat) * cos(dLat)
				+ (2.0 * m_dR * m_dR - m_dR - 2.0)
				- 2.0 * m_dR * m_dR / (cos(dLat) * cos(dLat)));

		double dCoeffB =
			2.0 * (dOmega + m_dW) * m_dK / ((m_dR + 1.0) * (m_dR + 2.0))
			* pow(cos(dLat), m_dR) * (
				(m_dR * m_dR + 2.0 * m_dR + 2.0)
				- (m_dR + 1.0) * (m_dR + 1.0) * cos(dLat) * cos(dLat));

		double dCoeffC =
			0.25 * m_dK * m_dK * pow(cos(dLat), 2.0 * m_dR)
			* ((m_dR + 1.0) * cos(dLat) * cos(dLat) - (m_dR + 2.0));

		double dGH = dG * m_dH0
			+ dEarthRadius * dEarthRadius * dCoeffA
			+ dEarthRadius * dEarthRadius * dCoeffB * cos(m_dR * dLon)
			+ dEarthRadius * dEarthRadius * dCoeffC * cos(2.0 * m_dR * dLon);

		dState[2] = (dGH / dG);

		// Velocity field
		dState[0] =
			m_dW * cos(dLat)
			+ m_dK * pow(cos(dLat), m_dR - 1.0)
			* cos(m_dR * dLon) * (
				m_dR * sin(dLat) * sin(dLat) - cos(dLat) * cos(dLat));

		dState[1] =
			- m_dK * m_dR * pow(cos(dLat), m_dR - 1.0)
			* sin(dLat) * sin(m_dR * dLon);

		dState[0] *= dEarthRadius;
		dState[1] *= dEarthRadius;
	}

};

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	// Initialize Tempest
	TempestInitialize(&argc, &argv);

try {
	// Background height field
	double dH0;

	// Wavenumber
	double dR;

	// W parameter
	double dW;

	// K parameter
	double dK;

	// Parse the command line
	BeginTempestCommandLine("RossbyHaurwitzWaveTest")
		SetDefaultResolution(16);
		SetDefaultLevels(1);
		SetDefaultOutputDeltaT("1d");
		SetDefaultDeltaT("480s");
		SetDefaultEndTime("15d");
		SetDefaultHorizontalOrder(4);
		SetDefaultVerticalOrder(1);

		CommandLineDouble(dH0, "h0", 8000.0);
		CommandLineDouble(dR, "r", 4.0);
		CommandLineDouble(dW, "w", 7.848e-6);
		CommandLineDouble(dK, "k", 7.848e-6);

		ParseCommandLine(argc, argv);
	EndTempestCommandLine(argv)

	// Setup the Model
	AnnounceBanner("MODEL SETUP");

	Model model(EquationSet::ShallowWaterEquations);

	TempestSetupCubedSphereModel(model);

	// Set the test case for the model
	AnnounceStartBlock("Initializing test case");
	model.SetTestCase(
		new RossbyHaurwitzWaveTest(dH0, dR, dW, dK));
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

	// Deinitialize
	TempestDeinitialize();
}

///////////////////////////////////////////////////////////////////////////////

