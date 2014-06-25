///////////////////////////////////////////////////////////////////////////////
///
///	\file    ShallowWaterEddyTest.cpp
///	\author  Paul Ullrich
///	\version October 4, 2013
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
///		Shallow water eddy test.
///	</summary>
class ShallowWaterEddyTest : public TestCase {

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
	///		Southern latitude of barotropic jet.
	///	</summary>
	double m_dTheta0;

	///	<summary>
	///		Northern latitude of barotropic jet.
	///	</summary>
	double m_dTheta1;

	///	<summary>
	///		Stretching parameter.
	///	</summary>
	double m_dXE;

	///	<summary>
	///		Parameter used for computing height field (??)
	///	</summary>
	double m_dHHat;

	///	<summary>
	///		Parameter used for computing height field (??)
	///	</summary>
	double m_dHPhi2;

	///	<summary>
	///		Parameter used for computing height field (??)
	///	</summary>
	double m_dHAlpha;

	///	<summary>
	///		Parameter used for computing height field (??)
	///	</summary>
	double m_dHBeta;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	ShallowWaterEddyTest(
		double dAlpha
	) :
		m_dH0(10158.18617045463179),
		m_dU0(80.0),
		m_dAlpha(dAlpha),
		m_dTheta0(M_PI / 7.0),
		m_dTheta1(M_PI / 2.0 - m_dTheta0),
		m_dXE(0.3),
		m_dHHat(120.0),
		m_dHPhi2(M_PI / 4.0),
		m_dHAlpha(1.0 / 3.0),
		m_dHBeta(1.0 / 15.0)
	{ }

public:
	///	<summary>
	///		Number of tracers used in this test.
	///	</summary>
	virtual int GetTracerCount() const {
		return 0;
	}

	///	<summary>
	///		Obtain test case specific physical constants.
	///	</summary>
	virtual void EvaluatePhysicalConstants(
		PhysicalConstants & phys
	) const {
		// Set the alpha parameter
		phys.SetAlpha(m_dAlpha);

		// Turn off Coriolis forces
		phys.SetOmega(0.0);
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
	///		Determined the corresponding point in rotated spherical coordinates.
	///	</summary>
	void CalculateRLLPrime(
		double dLon,
		double dLat,
		double & dLonP,
		double & dLatP
	) const {
		if (m_dAlpha == 0.0) {
			dLonP = dLon;
			dLatP = dLat;
			return;
		}

		// Calculate latitude and longitude in rotated coordinates
		dLatP = asin(sin(dLat) * cos(m_dAlpha)
			- cos(dLat) * cos(dLon) * sin(m_dAlpha));

		dLonP = asin(sin(dLon) * cos(dLat) / cos(dLatP));

		double dTemp =
			cos(m_dAlpha) * cos(dLon) * cos(dLat)
			+ sin(m_dAlpha) * sin(dLat);

		if (dTemp < 0.0) {
			dLonP = M_PI - dLonP;
		}
		if (dLonP < 0.0) {
			dLonP += 2.0 * M_PI;
		}
	}

	///	<summary>
	///		Evaluate the perturbed velocity field.
	///	</summary>
	double EvaluateUPrime(
		double dLonP,
		double dLatP
	) const {
		dLatP = fabs(dLatP);

		if (dLatP < m_dTheta0) {
			return 0.0;

		} else if (dLatP > m_dTheta1) {
			return 0.0;
		}

		double dNormalizer =
			exp(- 4.0 / (m_dTheta1 - m_dTheta0) / (m_dTheta1 - m_dTheta0));

		double dUp =
			exp(1.0 / (dLatP - m_dTheta0) / (dLatP - m_dTheta1));

		return (m_dU0 / dNormalizer * dUp);
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

		// Calculate rotated RLL coordinates of this point
		double dLonP;
		double dLatP;

		CalculateRLLPrime(dLon, dLat, dLonP, dLatP);

		// Calculate height field via numerical integration
		int nIntervals =
			static_cast<unsigned int>((dLatP + 0.5 * M_PI) / (1.0e-2));

		if (nIntervals < 1) {
			nIntervals = 1;
		}

		// Numerically integrate H using 4th order Gaussian quadrature
		DataVector<double> dLatX;
		dLatX.Initialize(nIntervals+1);

		for (int i = 0; i <= nIntervals; i++) {
			dLatX[i] = - 0.5 * M_PI +
				((dLatP + 0.5 * M_PI) / static_cast<double>(nIntervals))
					* static_cast<double>(i);
		}

		double dH = 0.0;

		double dXeval;
		double dU;

		for (int i = 0; i < nIntervals; i++) {
		for (int m = -1; m <= 1; m += 2) {
			dXeval =
				0.5 * (dLatX[i+1] + dLatX[i])
				+ static_cast<double>(m)
					* sqrt(1.0 / 3.0) * 0.5 * (dLatX[i+1] - dLatX[i]);

			dU = EvaluateUPrime(dLonP, dXeval);

			dH += (2.0 * phys.GetEarthRadius() * phys.GetOmega() * sin(dXeval)
				+ dU * tan(dXeval)) * dU;
		}
		}

		dH *= 0.5 * (dLatX[1] - dLatX[0]);

		dState[2] = m_dH0 - dH / phys.GetG();

		// Add perturbation
		if (dLon > M_PI) {
			dLon = dLon - 2.0 * M_PI;
		}
		if ((dLon < -M_PI) || (dLon > M_PI)) {
			_EXCEPTIONT("Invalid value of longitude.");
		}

		dState[2] += m_dHHat * cos(dLat)
			* exp(-((dLon * dLon) / (m_dHAlpha * m_dHAlpha)))
			* exp(-((m_dHPhi2 - dLat) * (m_dHPhi2 - dLat)
				/ (m_dHBeta * m_dHBeta)));

		// Evaluate the velocity field
		double dUP = EvaluateUPrime(dLonP, dLatP);

		double dUlat =
			(- dUP * sin(m_dAlpha) * sin(dLonP)) / cos(dLat);

		double dUlon;
		if (fabs(cos(dLon)) < 1.0e-13) {
			if (fabs(m_dAlpha) > 1.0e-13) {
				if (cos(dLon) > 0.0) {
					dUlon = - dUlat * cos(dLat) / tan(m_dAlpha);
				} else {
					dUlon = dUlat * cos(dLat) / tan(m_dAlpha);
				}
			} else {
				dUlon = dUP;
			}

		} else {
			dUlon = (dUlat * sin(dLat) * sin(dLon) + dUP * cos(dLonP))
				/ cos(dLon);
		}

		dState[0] = dUlon;
		dState[1] = dUlat;
	}

};

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	// Initialize Tempest
	TempestInitialize(&argc, &argv);

try {

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
		new ShallowWaterEddyTest(dAlpha));
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

