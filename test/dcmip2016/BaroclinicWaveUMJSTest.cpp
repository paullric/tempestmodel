///////////////////////////////////////////////////////////////////////////////
///
///	\file    BaroclinicWaveUMJSTest.cpp
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

extern "C" {
	void tc_baroclinic_(
                      double * dLon,
                      double * dLat,
                      double * dP,
                      double * dZ,
                      double * dU,
                      double * dV,
                      double * dW,
                      double * dT,
                      double * dPhis,
                      double * dPs,
                      double * dRho,
                      double * dQ
                      );
}



///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Ullrich, Melvin, Jablonowski and Staniforth (2013) Baroclinic wave test
///	</summary>
class BaroclinicWaveUMJSTest : public TestCase {

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
	///		Model scaling parameter
	///	</summary>
	const double ParamEarthRadiusScaling;

	///	<summary>
	///		Model height limit
	///	</summary>
	const double ParamHeightLimit;

	///	<summary>
	///		Equatorial surface temperature (K)
	///	</summary>
	const double ParamT0E;

	///	<summary>
	///		Polar surface temperature (K)
	///	</summary>
	const double ParamT0P;

	///	<summary>
	///		Half-width parameter
	///	</summary>
	const double ParamB;

	///	<summary>
	///		Jet parameter
	///	</summary>
	const double ParamK;

	///	<summary>
	///		Temperature lapse rate
	///	</summary>
	const double ParamLapseRate;

	///	<summary>
	///		Stream function perturbation wind parameter
	///	</summary>
	const double ParamU0;

	///	<summary>
	///		Stream function perturbation radius (Earth radii)
	///	</summary>
	const double ParamPertR;

	///	<summary>
	///		Expontential zonal wind perturbation
	///	</summary>
	const double ParamUp;

	///	<summary>
	///		Exponential perturbation radius (Earth radii)
	///	</summary>
	const double ParamPertExpR;

	///	<summary>
	///		Perturbation longitude center (radians)
	///	</summary>
	const double ParamPertLon;

	///	<summary>
	///		Perturbation latitude center (radians)
	///	</summary>
	const double ParamPertLat;

	///	<summary>
	///		Height cap of the perturbation (polynomial tapering)
	///	</summary>
	const double ParamPertZ;

protected:
	///	<summary>
	///		Alpha parameter.
	///	</summary>
	double m_dAlpha;

	///	<summary>
	///		Deep atmosphere flag.
	///	</summary>
	bool m_fDeepAtmosphere;

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
	BaroclinicWaveUMJSTest(
		double dAlpha,
		bool fDeepAtmosphere,
		double dZtop,
		PerturbationType ePerturbationType = PerturbationType_None
	) :
		ParamEarthRadiusScaling(1.0),
		ParamHeightLimit(30000.0),
		ParamT0E(310.0),
		ParamT0P(240.0),
		ParamB(2.0),
		ParamK(3.0),
		ParamLapseRate(0.005),
		ParamU0(-0.5),
		ParamPertR(1.0 / 6.0),
		ParamUp(1.0),
		ParamPertExpR(0.1),
		ParamPertLon(M_PI / 9.0),
		ParamPertLat(2.0 * M_PI / 9.0),
		ParamPertZ(15000.0),

		m_dAlpha(dAlpha),
		m_fDeepAtmosphere(fDeepAtmosphere),
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
		return (0.0);
	}

	///	<summary>
	///		Evaluate the stream function perturbation at the given point (in
	///		latitude-longitude coordinates).
	///	</summary>
	double EvaluateStreamFunction(
		double dZ,
		double dLon,
		double dLat
	) const {

		// Verify that the stream function perturbation is being used
		if (m_ePerturbationType != PerturbationType_StreamFn) {
			_EXCEPTION();
		}

		// Calculate great circle distance
		double dGreatCircleR =
			acos(sin(ParamPertLat) * sin(dLat)
				+ cos(ParamPertLat) * cos(dLat) * cos(dLon - ParamPertLon));

		dGreatCircleR /= ParamPertR;

		// Tapered perturbation with height
		double dPertTaper = 0.0;
		if (dZ < ParamPertZ) {
			dPertTaper = 1.0
				- 3.0 * dZ * dZ / (ParamPertZ * ParamPertZ)
				+ 2.0 * dZ * dZ * dZ / (ParamPertZ * ParamPertZ * ParamPertZ);
		} else {
			dPertTaper = 0.0;
		}

		// Calculate stream function
		double dCosPert;
		if (dGreatCircleR < 1.0) {
			dCosPert = cos(0.5 * M_PI * dGreatCircleR);
		} else {
			dCosPert = 0.0;
		}

		double dStreamFn =
			ParamU0 * ParamPertR * dPertTaper
				* dCosPert * dCosPert * dCosPert * dCosPert;

		return dStreamFn;
	}

	///	<summary>
	///		Evaluate the pointwise perturbation at the given point.
	///	</summary>
	void EvaluatePointwisePerturbation(
		const PhysicalConstants & phys,
		double dZ,
		double dLon,
		double dLat,
		double & dUlon,
		double & dUlat
	) const {

		// A small value for numerical derivatives
		const double Epsilon = 1.0e-5;

		// No perturbation
		if (m_ePerturbationType == PerturbationType_None) {
			dUlon = 0.0;
			dUlat = 0.0;
		}

		// Exponential perturbation
		if (m_ePerturbationType == PerturbationType_Exp) {
			// Calculate great circle distance
			double dGreatCircleR =
				acos(sin(ParamPertLat) * sin(dLat)
					+ cos(ParamPertLat) * cos(dLat) * cos(dLon - ParamPertLon));

			dGreatCircleR /= ParamPertExpR;

	  	 	// Tapered perturbation with height
   			double dPertTaper = 0.0;
			if (dZ < ParamPertZ) {
				dPertTaper = 1.0
 					- 3.0 * dZ * dZ / (ParamPertZ * ParamPertZ)
					+ 2.0 * dZ * dZ * dZ / (ParamPertZ * ParamPertZ * ParamPertZ);
			} else {
				dPertTaper = 0.0;
			}

			// Apply perturbation in zonal velocity
			if (dGreatCircleR < 1.0) {
				dUlon = ParamUp * dPertTaper
					* exp(- dGreatCircleR * dGreatCircleR);
			} else {
				dUlon = 0.0;
			}

			dUlat = 0.0;
		}

		// Stream function perturbation
		if (m_ePerturbationType == PerturbationType_StreamFn) {
			// Evaluate the perturbation in zonal velocity
			dUlon = - 1.0 / (2.0 * Epsilon) * (
				  EvaluateStreamFunction(dZ, dLon, dLat + Epsilon)
				- EvaluateStreamFunction(dZ, dLon, dLat - Epsilon));

			dUlat = 1.0 / (2.0 * Epsilon * cos(dLat)) * (
				  EvaluateStreamFunction(dZ, dLon + Epsilon, dLat)
				- EvaluateStreamFunction(dZ, dLon - Epsilon, dLat));
		}
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

		// Radius
		double dR = dZ + phys.GetEarthRadius();

		// Calculate parameters
		double dT0 = 0.5 * (ParamT0E + ParamT0P);

		double dConstA = 1.0 / ParamLapseRate;

		double dConstB = (dT0 - ParamT0P) / (dT0 * ParamT0P);

		double dConstC = 0.5 * (ParamK + 2.0)
			* (ParamT0E - ParamT0P) / (ParamT0E * ParamT0P);

		double dConstH = phys.GetR() * dT0 / phys.GetG();

		// Computed quantities
		double dScaledZ = dZ / (ParamB * dConstH);

		// Calculate tau values
		double dTau1 =
			dConstA * ParamLapseRate / dT0
				* exp(ParamLapseRate / dT0 * dZ)
			+ dConstB
				* (1.0 - 2.0 * dScaledZ * dScaledZ)
				* exp(- dScaledZ * dScaledZ);

		double dTau2 =
			dConstC * (1.0 - 2.0 * dScaledZ * dScaledZ)
				* exp(- dScaledZ * dScaledZ);

		double dIntTau1 =
			dConstA * (exp(ParamLapseRate / dT0 * dZ) - 1.0)
			+ dConstB * dZ * exp(- dScaledZ * dScaledZ);

		double dIntTau2 =
			dConstC * dZ * exp(- dScaledZ * dScaledZ);

		// Calculate utility terms
		double dRRatio;
		if (m_fDeepAtmosphere) {
			dRRatio = dR / phys.GetEarthRadius();
		} else {
			dRRatio = 1.0;
		}

		double dInteriorTerm = pow(dRRatio * cos(dLat), ParamK)
			- ParamK / (ParamK + 2.0) * pow(dRRatio * cos(dLat), ParamK + 2.0);

		// Calculate temperature
		double dTemperature = 1.0
			/ (dRRatio * dRRatio)
			/ (dTau1 - dTau2 * dInteriorTerm);

		// Calculate hydrostatic pressure
		double dPressure = phys.GetP0() * exp(
			- phys.GetG() / phys.GetR()
				* (dIntTau1 - dIntTau2 * dInteriorTerm));
/*
		// Calculate pressure derivative
		double dDrPressure;
		if (m_fDeepAtmosphere) {
			dDrPressure =
				dPressure * phys.GetG() / phys.GetR() * (
					- dTau1 + dTau2 * dInteriorTerm
					+ dIntTau2 * ParamK / dR
						* (pow(dRRatio * cos(dLat), ParamK)
							- pow(dRRatio * cos(dLat), ParamK + 2.0)));

		} else {
			dDrPressure =
				dPressure * phys.GetG() / phys.GetR() * (
					- dTau1 + dTau2 * dInteriorTerm);
		}

		// Calculate hydrostatic density
		double dHydroRho = - dRRatio * dRRatio / phys.GetG() * dDrPressure;
*/
		// Calculate exact density
		double dRho = dPressure / (phys.GetR() * dTemperature);

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
        double dRho;
        double dU;
        double dV;
        double dW;
        double dP;
        double dT;
        double dPhis;
        double dPs;
        double dQ;
        
        tc_baroclinic_(&dLon, &dLat,&dP,&dZ,&dU,&dV,&dW,&dT,&dPhis,&dPs,&dRho,&dQ);

		// Radius
		double dR = dZ + phys.GetEarthRadius();

		// Calculate parameters
		double dT0 = 0.5 * (ParamT0E + ParamT0P);

		double dConstA = 1.0 / ParamLapseRate;

		double dConstB = (dT0 - ParamT0P) / (dT0 * ParamT0P);

		double dConstC = 0.5 * (ParamK + 2.0)
			* (ParamT0E - ParamT0P) / (ParamT0E * ParamT0P);

		double dConstH = phys.GetR() * dT0 / phys.GetG();

		// Computed quantities
		double dScaledZ = dZ / (ParamB * dConstH);

		// Calculate tau values
		double dTau1 =
			dConstA * ParamLapseRate / dT0
				* exp(ParamLapseRate / dT0 * dZ)
			+ dConstB
				* (1.0 - 2.0 * dScaledZ * dScaledZ)
				* exp(- dScaledZ * dScaledZ);

		double dTau2 =
			dConstC * (1.0 - 2.0 * dScaledZ * dScaledZ)
				* exp(- dScaledZ * dScaledZ);

		double dIntTau1 =
			dConstA * (exp(ParamLapseRate / dT0 * dZ) - 1.0)
			+ dConstB * dZ * exp(- dScaledZ * dScaledZ);

		double dIntTau2 =
			dConstC * dZ * exp(- dScaledZ * dScaledZ);

		// Calculate utility terms
		double dRRatio;
		if (m_fDeepAtmosphere) {
			dRRatio = dR / phys.GetEarthRadius();
		} else {
			dRRatio = 1.0;
		}

		double dInteriorTerm = pow(dRRatio * cos(dLat), ParamK)
			- ParamK / (ParamK + 2.0) * pow(dRRatio * cos(dLat), ParamK + 2.0);

		// Calculate temperature
		double dTemperature = 1.0
			/ (dRRatio * dRRatio)
			/ (dTau1 - dTau2 * dInteriorTerm);

		// Calculate hydrostatic pressure
		double dPressure = phys.GetP0() * exp(
			- phys.GetG() / phys.GetR()
				* (dIntTau1 - dIntTau2 * dInteriorTerm));

		// Calculate hydrostatic density
		/*double*/ dRho = dPressure / (phys.GetR() * dTemperature);

		// Velocity field
		double dInteriorTermU =
			  pow(dRRatio * cos(dLat), ParamK - 1.0)
			- pow(dRRatio * cos(dLat), ParamK + 1.0);

		double dBigU = phys.GetG() / phys.GetEarthRadius() * ParamK
			* dIntTau2 * dInteriorTermU * dTemperature;

		double dRCosLat;
		if (m_fDeepAtmosphere) {
			dRCosLat = dR * cos(dLat);
		} else {
			dRCosLat = phys.GetEarthRadius() * cos(dLat);
		}

		double dOmegaRCosLat = phys.GetOmega() * dRCosLat;

		if (dOmegaRCosLat * dOmegaRCosLat + dRCosLat * dBigU < 0.0) {
			_EXCEPTIONT("Negative discriminant detected.");
		}

		double dUlon = - dOmegaRCosLat +
			sqrt(dOmegaRCosLat * dOmegaRCosLat + dRCosLat * dBigU);

		double dUlat = 0.0;

		// Calculate velocity perturbation
		double dUlonPert;
		double dUlatPert;

		EvaluatePointwisePerturbation(
			phys, dZ, dLon, dLat,
			dUlonPert, dUlatPert);

		dUlon += dUlonPert;
		dUlat += dUlatPert;

		// Store the state
		dState[0] = dU;
		dState[1] = dV;
		dState[2] = dT*pow((100000.0/dP),(phys.GetR()/phys.GetCp()));
		dState[3] = 0.0;
		dState[4] = dRho;

	}
};

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	// Initialize Tempest
	TempestInitialize(&argc, &argv);

try {
	// Model height cap
	double dZtop;

	// Grid rotation angle
	double dAlpha;

	// Include tracer field
	bool fTracersOn;

	// Deep atmosphere flag
	bool fDeepAtmosphere;

	// Perturbation type
	std::string strPerturbationType;

	// Parse the command line
	BeginTempestCommandLine("BaroclinicWaveUMJS");
		SetDefaultResolution(20);
		SetDefaultLevels(10);
		SetDefaultOutputDeltaT("200s");
		SetDefaultDeltaT("200s");
		SetDefaultEndTime("200s");
		SetDefaultHorizontalOrder(4);
		SetDefaultVerticalOrder(1);

		CommandLineDouble(dZtop, "ztop", 10000.0);
		CommandLineDouble(dAlpha, "alpha", 0.0);
		CommandLineBool(fDeepAtmosphere, "deep_atmosphere");
		CommandLineStringD(strPerturbationType, "pert",
			"None", "(None | Exp | Sfn)");

		ParseCommandLine(argc, argv);
	EndTempestCommandLine(argv)

	// Setup the Model
	AnnounceBanner("MODEL SETUP");

	Model model(EquationSet::PrimitiveNonhydrostaticEquations);

	TempestSetupCubedSphereModel(model);

	// Set the test case for the model
	AnnounceStartBlock("Initializing test case");

	BaroclinicWaveUMJSTest::PerturbationType ePerturbationType;
	STLStringHelper::ToLower(strPerturbationType);
	if (strPerturbationType == "none") {
		ePerturbationType = BaroclinicWaveUMJSTest::PerturbationType_None;
	} else if (strPerturbationType == "exp") {
		ePerturbationType = BaroclinicWaveUMJSTest::PerturbationType_Exp;
	} else if (strPerturbationType == "sfn") {
		ePerturbationType = BaroclinicWaveUMJSTest::PerturbationType_StreamFn;
	} else {
		_EXCEPTIONT("Invalid perturbation type:"
			" Expected \"None\", \"Exp\" or \"SFn\"");
	}

	model.SetTestCase(
		new BaroclinicWaveUMJSTest(
			dAlpha,
			fDeepAtmosphere,
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

