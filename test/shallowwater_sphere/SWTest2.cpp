///////////////////////////////////////////////////////////////////////////////
///
///	\file    SWTest2.cpp
///	\author  Paul Ullrich
///	\version February 24, 2013
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
///		Williamson et al. (1994) Test case 2
///
///		Geostrophically balanced flow with a shallow-water model.
///	</summary>
class ShallowWaterTestCase2 : public TestCase {

private:
	///	<summary>
	///		Background height field.
	///	</summary>
	double m_dH0;

	///	<summary>
	///		Maximum velocity.
	///	</summary>
	double m_dU0;

	///	<summary>
	///		Grid inclination.
	///	</summary>
	double m_dAlpha;

	///	<summary>
	///		True if a cosine bell tracer should be included.
	///	</summary>
	bool m_fTracerOn;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	ShallowWaterTestCase2(
		double dH0,
		double dU0,
		double dAlpha
	) :
		m_dH0(dH0),
		m_dU0(dU0),
		m_dAlpha(dAlpha * M_PI / 180.0),
		m_fTracerOn(false)
	{ }

public:
	///	<summary>
	///		Number of tracers used in this test.
	///	</summary>
	int GetTracerCount() const {
		return (m_fTracerOn)?(1):(0);
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

		// Pointwise zonal component of velocity
		dState[0] = m_dU0 * cos(dLat)
			* (cos(m_dAlpha) + cos(dLon) * tan(dLat) * sin(m_dAlpha));

		// Pointwise meridional component of velocity
		dState[1] = - m_dU0 * sin(dLon) * sin(m_dAlpha);
/*
		dState[0] += cos(2.0 * dLon) * cos(dLat) * cos(dLat);
		dState[1] += sin(2.0 * dLat);
*/
		// Height field
		double dHTrig =
			- cos(dLon) * cos(dLat) * sin(m_dAlpha)
			+ sin(dLat) * cos(m_dAlpha);

		dState[2] = m_dH0
			- (phys.GetEarthRadius() * phys.GetOmega() + 0.5 * m_dU0)
				* m_dU0 * dHTrig * dHTrig / phys.GetG();
/*
		dState[0] = 100.0;

		double dXX = cos(dLat) * cos(dLon) - 1.0;
		double dYY = cos(dLat) * sin(dLon);
		double dZZ = sin(dLat);

		double dR = (dXX * dXX + dYY * dYY + dZZ * dZZ);

		dState[0] += 100.0 * exp(-2.0 * dR);
*/

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
	}

};

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

try {
	// Initialize Tempest
	TempestInitialize(&argc, &argv);

	// Background velocity
	double dU0;

	// Background height field
	double dH0;

	// Grid rotation angle
	double dAlpha;

	// Parse the command line
	BeginTempestCommandLine("SWTest2")
		SetDefaultResolution(40);
		SetDefaultLevels(1);
		SetDefaultOutputDeltaT("200s");
		SetDefaultDeltaT("200s");
		SetDefaultEndTime("200s");
		SetDefaultHorizontalOrder(4);
		SetDefaultVerticalOrder(1);

		CommandLineDouble(dU0, "u0", 38.61068277);
		CommandLineDouble(dH0, "h0", 2998.104995);
		CommandLineDouble(dAlpha, "alpha", 0.0);

		ParseCommandLine(argc, argv);
	EndTempestCommandLine(argv)

	// Setup the Model
	AnnounceBanner("MODEL SETUP");

	Model model(EquationSet::ShallowWaterEquations);

	TempestSetupCubedSphereModel(model);

	// Set the test case for the model
	AnnounceStartBlock("Initializing test case");
	model.SetTestCase(
		new ShallowWaterTestCase2(dH0, dU0, dAlpha));
	AnnounceEndBlock("Done");

	// Begin execution
	AnnounceBanner("SIMULATION");
	model.Go();

	// Compute error norms
	AnnounceBanner("RESULTS");
	model.ComputeErrorNorms();
	AnnounceBanner();

	// Deinitialize
	TempestDeinitialize();

} catch(Exception & e) {
	std::cout << e.ToString() << std::endl;
}
}

///////////////////////////////////////////////////////////////////////////////

