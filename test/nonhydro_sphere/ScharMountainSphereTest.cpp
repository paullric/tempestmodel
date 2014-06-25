///////////////////////////////////////////////////////////////////////////////
///
///	\file    ScharMountainSphereTest.cpp
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
///		Mountain waves induced by a Schar-type mountain (DCMIP 2012 test 2-x)
///	</summary>
class ScharMountainSphereTest : public TestCase {

public:
	///	<summary>
	///		Perturbation type.
	///	</summary>
	enum MountainType {
		MountainType_Default = 0,
		MountainType_None = MountainType_Default,
		MountainType_Wave6 = 1,
	};

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
	///		Longitude of Schar-type mountain centerpoint.
	///	</summary>
	double m_dLonC;

	///	<summary>
	///		Latitude of Schar-type mountain centerpoint.
	///	</summary>
	double m_dLatC;

	///	<summary>
	///		Maximum Schar-type mountain height.
	///	</summary>
	double m_dH0;

	///	<summary>
	///		Schar-type mountain half width.
	///	</summary>
	double m_dD;

	///	<summary>
	///		Schar-type mountain wavelength.
	///	</summary>
	double m_dXi;

	///	<summary>
	///		Reference surface temperature at the equator.
	///	</summary>
	double m_dTeq;

	///	<summary>
	///		Reference zonal wind velocity.
	///	</summary>
	double m_dUeq;

	///	<summary>
	///		Equatorial surface wind shear.
	///	</summary>
	double m_dCs;

	///	<summary>
	///		Height of the Rayleigh damped layer.
	///	</summary>
	double m_dZh;

	///	<summary>
	///		Rayleigh friction time scale.
	///	</summary>
	double m_dTau0;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	ScharMountainSphereTest(
		double dZtop,
		double dEarthScaling,
		double dOmega,
		double dLonC,
		double dLatC,
		double dH0,
		double dD,
		double dXi,
		double dTeq,
		double dUeq,
		double dCs,
		double dZh,
		double dTau0
	) :
		m_dZtop(dZtop),
		m_dEarthScaling(dEarthScaling),
		m_dOmega(dOmega),
		m_dLonC(dLonC * M_PI / 180.0),
		m_dLatC(dLatC * M_PI / 180.0),
		m_dH0(dH0),
		m_dD(dD),
		m_dXi(dXi),
		m_dTeq(dTeq),
		m_dUeq(dUeq),
		m_dCs(dCs),
		m_dZh(dZh),
		m_dTau0(dTau0)
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

		// Great circle distance from mountain centerpoint
		double dR = phys.GetEarthRadius() * acos(
			sin(m_dLatC) * sin(dLat)
				+ cos(m_dLatC) * cos(dLat) * cos(dLon - m_dLonC));

		// Topography height
		double dCosTerm = cos(M_PI * dR / m_dXi);

		double dExpTerm = exp(- dR * dR / (m_dD * m_dD));

		return (m_dH0 * dExpTerm * dCosTerm * dCosTerm);
	}

	///	<summary>
	///		Flag indicating whether or not Rayleigh friction strength is given.
	///	</summary>
	virtual bool HasRayleighFriction() const {
		return true;
	}

	///	<summary>
	///		Evaluate the Rayleigh friction strength at the given point.
	///	</summary>
	virtual double EvaluateRayleighStrength(
		double dZ,
		double dXp,
		double dYp
	) const {
		double dNuDepth = 0.0;

		if (dZ > m_dZh) {
			double dNormZ = (dZ - m_dZh) / (m_dZtop - m_dZh);

			dNuDepth = sin(M_PI / 2.0 * dNormZ);

			dNuDepth = dNuDepth * dNuDepth;
		}

		return (dNuDepth / m_dTau0);
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
		double dT = m_dTeq
			* (1.0 - m_dCs * m_dUeq * m_dUeq / phys.GetG() * dSin2Lat);

		// 3D pressure
		double dPressure = phys.GetP0() * exp(
			- m_dUeq * m_dUeq / (2.0 * phys.GetR() * m_dTeq) * dSin2Lat
			- phys.GetG() * dZ / (phys.GetR() * dT));

		// 3D density
		double dRho = dPressure / (phys.GetR() * dT);

		// Zonal velocity
		double dU = m_dUeq * cos(dLat)
			* sqrt(2.0 * m_dTeq / dT * m_dCs * dZ + dT / m_dTeq);

		// Store the state
		dState[0] = dU;
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

	// Longitude of Schar-type mountain centerpoint.
	double dLonC;

	// Latitude of Schar-type mountain centerpoint.
	double dLatC;

	// Maximum Schar-type mountain height.
	double dH0;

	// Schar-type mountain half width.
	double dD;

	// Schar-type mountain wavelength.
	double dXi;

	// Reference surface temperature at the equator.
	double dTeq;

	// Reference zonal wind velocity.
	double dUeq;

	// Equatorial surface wind shear.
	double dCs;

	// Height of the Rayleigh damped layer.
	double dZh;

	// Rayleigh friction time scale.
	double dTau0;

	// Parse the command line
	BeginTempestCommandLine("ScharMountainSphereTest");
		SetDefaultResolution(20);
		SetDefaultLevels(60);
		SetDefaultOutputDeltaT("600s");
		SetDefaultDeltaT("500000u");
		SetDefaultEndTime("7200s");
		SetDefaultHorizontalOrder(4);
		SetDefaultVerticalOrder(1);

		CommandLineDouble(dZtop, "ztop", 30000.0);
		CommandLineDouble(dEarthScaling, "X", 500.0);
		CommandLineDouble(dOmega, "omega", 0.0);
		CommandLineDouble(dLonC, "lonc", 45.0);
		CommandLineDouble(dLatC, "latc", 0.0);
		CommandLineDouble(dH0, "h0", 250.0);
		CommandLineDouble(dD, "d", 5000.0);
		CommandLineDouble(dXi, "xi", 4000.0);
		CommandLineDouble(dTeq, "teq", 300.0);
		CommandLineDouble(dUeq, "ueq", 20.0);
		CommandLineDoubleD(dCs, "cs", 0.0, "(for sheared flow 2.5e-4)");
		CommandLineDouble(dZh, "zh", 20000.0);
		CommandLineDouble(dTau0, "tau0", 25.0);

		ParseCommandLine(argc, argv);
	EndTempestCommandLine(argv)

	// Setup the Model
	AnnounceBanner("MODEL SETUP");

	Model model(EquationSet::PrimitiveNonhydrostaticEquations);

	TempestSetupCubedSphereModel(model);

	// Set the test case for the model
	AnnounceStartBlock("Initializing test case");

	model.SetTestCase(
		new ScharMountainSphereTest(
			dZtop,
			dEarthScaling,
			dOmega,
			dLonC,
			dLatC,
			dH0,
			dD,
			dXi,
			dTeq,
			dUeq,
			dCs,
			dZh,
			dTau0));

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

