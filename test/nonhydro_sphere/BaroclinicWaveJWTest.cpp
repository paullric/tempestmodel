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

#include "Tempest.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Jablonowski and Williamson (2006) Baroclinic wave test
///	</summary>
class BaroclinicWaveJWTest : public TestCase {

public:
	///	<summary>
	///		Perturbation type.
	///	</summary>
	enum PerturbationType {
		PerturbationType_Default = 0,
		PerturbationType_None = PerturbationType_Default,
		PerturbationType_Exp = 1,
		PerturbationType_StreamFn = 2,
	};

public:
	///	<summary>
	///		Auxiliary eta.
	///	</summary>
	const double ParamEta0;

	///	<summary>
	///		Tropopause level (in eta coordinates)
	///	</summary>
	const double ParamTropopauseEta;

	///	<summary>
	///		Horizontal-mean temperature (K)
	///	</summary>
	const double ParamT0;

	///	<summary>
	///		Empirical temperature difference (K)
	///	</summary>
	const double ParamDeltaT;

	///	<summary>
	///		Temperature lapse rate (K/m)
	///	</summary>
	const double ParamLapseRate;

	///	<summary>
	///		Maximum zonal wind (m/s)
	///	</summary>
	const double ParamU0;

	///	<summary>
	///		Zonal wind perturbation (m / s)
	///	</summary>
	const double ParamUp;

	///	<summary>
	///		Perturbation longitude center (radians)
	///	</summary>
	const double ParamPertLon;

	///	<summary>
	///		Perturbation latitude center (radians)
	///	</summary>
	const double ParamPertLat;

	///	<summary>
	///		Perturbation radius (Earth radii)
	///	</summary>
	const double ParamPertR;

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
	///		Type of perturbation.
	///	</summary>
	PerturbationType m_ePerturbationType;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	BaroclinicWaveJWTest(
		double dAlpha,
		double dZtop,
		PerturbationType ePerturbationType = PerturbationType_None
	) :
		ParamEta0(0.252),
		ParamTropopauseEta(0.2),
		ParamT0(288.0),
		ParamDeltaT(4.8e5),
		ParamLapseRate(0.005),
		ParamU0(35.0),
		ParamUp(1.0),
		ParamPertLon(M_PI / 9.0),
		ParamPertLat(2.0 * M_PI / 9.0),
		ParamPertR(0.1),

		m_dAlpha(dAlpha),
		m_fTracerOn(false),
		m_dZtop(dZtop),
		m_ePerturbationType(ePerturbationType)
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
		return true;
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
		const Time & time,
		double dZ,
		double dLon,
		double dLat,
		double * dState,
		double * dTracer
	) const {

		// Evaluate the reference state at this point
		EvaluateReferenceState(phys, dZ, dLon, dLat, dState);

		// Add perturbation in zonal velocity
		if (m_ePerturbationType == PerturbationType_Exp) {
			double dGreatCircleR =
				acos(sin(ParamPertLat) * sin(dLat)
					+ cos(ParamPertLat) * cos(dLat) * cos(dLon - ParamPertLon));

			dGreatCircleR /= ParamPertR;

			if (dGreatCircleR < 1.0) {
				dState[0] += ParamUp * exp( - dGreatCircleR * dGreatCircleR);
			}
		}
	}

};

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	// Initialize MPI
	TempestInitialize(&argc, &argv);

try {
	// Model height cap
	double dZtop;

	// Grid rotation angle
	double dAlpha;

	// Tracers turned on
	bool fTracersOn;

	// Perturbation type
	std::string strPerturbationType;

	// Parse the command line
	BeginTempestCommandLine("BaroclinicWaveJW");
		SetDefaultResolution(20);
		SetDefaultLevels(10);
		SetDefaultOutputDeltaT("200s");
		SetDefaultDeltaT("200s");
		SetDefaultEndTime("200s");
		SetDefaultHorizontalOrder(4);
		SetDefaultVerticalOrder(1);

		CommandLineDouble(dZtop, "ztop", 10000.0);
		CommandLineDouble(dAlpha, "alpha", 0.0);
		CommandLineStringD(strPerturbationType, "pert",
			"None", "(None | Exp)");

		ParseCommandLine(argc, argv);
	EndTempestCommandLine(argv)

	// Setup the Model
	AnnounceBanner("MODEL SETUP");

	Model model(EquationSet::PrimitiveNonhydrostaticEquations);

	TempestSetupCubedSphereModel(model);

	// Set the test case for the model
	AnnounceStartBlock("Initializing test case");

	BaroclinicWaveJWTest::PerturbationType ePerturbationType;
	STLStringHelper::ToLower(strPerturbationType);
	if (strPerturbationType == "none") {
		ePerturbationType = BaroclinicWaveJWTest::PerturbationType_None;
	} else if (strPerturbationType == "exp") {
		ePerturbationType = BaroclinicWaveJWTest::PerturbationType_Exp;
	} else {
		_EXCEPTIONT("Invalid perturbation type:"
			" Expected \"None\" or \"Exp\"");
	}

	model.SetTestCase(
		new BaroclinicWaveJWTest(
			dAlpha,
			dZtop,
			ePerturbationType));

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

