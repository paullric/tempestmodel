///////////////////////////////////////////////////////////////////////////////
///
///	\file    MountainWaveSphereTest.cpp
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
///		Mountain waves on the sphere test.
///	</summary>
class MountainWaveSphereTest : public TestCase {

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
	///		Background temperature (isothermal).
	///	</summary>
	double m_dT0;

	///	<summary>
	///		Background flow velocity (isothermal).
	///	</summary>
	double m_dU0;

	///	<summary>
	///		Flag indicating that rotation is disabled.
	///	</summary>
	bool m_fNoRotation;

	///	<summary>
	///		Type of mountain.
	///	</summary>
	MountainType m_eMountainType;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	MountainWaveSphereTest(
		double dZtop,
		double dEarthScaling,
		double dT0,
		double dU0,
		bool fNoRotation,
		MountainType eMountainType = MountainType_Default
	) :
		m_dZtop(dZtop),
		m_dEarthScaling(dEarthScaling),
		m_dT0(dT0),
		m_dU0(dU0),
		m_fNoRotation(fNoRotation),
		m_eMountainType(eMountainType)
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
		if (m_fNoRotation) {
			phys.SetOmega(0.0);
		}
	}

	///	<summary>
	///		Evaluate the topography at the given point.
	///	</summary>
	virtual double EvaluateTopography(
		const PhysicalConstants & phys,
		double dLon,
		double dLat
	) const {
		if (m_eMountainType == MountainType_None) {
			return 0.0;

		} else if (m_eMountainType == MountainType_Wave6) {
			return 10.0 * sin(6.0 * dLon) * cos(dLat) * cos(dLat);
		}

		_EXCEPTIONT("Invalid MountainType");
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
		const double dRayleighStrength = 8.0e-3;
		const double dRayleighDepth = 6000.0;

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
		double dXi,
		double dZ,
		double dLon,
		double dLat,
		double * dState
	) const {

		// Scale height
		double dH = phys.GetR() * m_dT0 / phys.GetG();

		// Froude number squared
		double dFr2 = m_dU0 * m_dU0 / (phys.GetG() * dH);

		// Inverse Rossby number
		double dInvRo = 2.0 * phys.GetEarthRadius() * phys.GetOmega() / m_dU0;

		// Pressure (isothermal)
		double dPressure =
			phys.GetP0()
				* exp(- dZ / dH)
				* exp(- 0.5 * dFr2 * (1.0 + dInvRo) * sin(dLat) * sin(dLat));

		// Density (isothermal)
		double dRho = dPressure / (phys.GetG() * dH);

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
		return EvaluateReferenceState(phys, dXi, dZ, dLon, dLat, dState);
	}
};

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	// Initialize MPI
	TempestInitialize(&argc, &argv);

try {
	// Model cap
	double dZtop;

	// Earth scaling parameter
	double dEarthScaling;

	// Background temperature
	double dT0;

	// Background wind speed
	double dU0;

	// No planetary rotation
	bool fNoRotation;

	// Mountain type
	std::string strMountainType;

	// Parse the command line
	BeginTempestCommandLine("MountainWaveSphereTest");
		SetDefaultResolution(20);
		SetDefaultLevels(10);
		SetDefaultOutputDeltaT("200s");
		SetDefaultDeltaT("200s");
		SetDefaultEndTime("200s");
		SetDefaultHorizontalOrder(4);
		SetDefaultVerticalOrder(1);

		CommandLineDouble(dZtop, "ztop", 10000.0);
		CommandLineDouble(dEarthScaling, "X", 1.0);
		CommandLineDouble(dT0, "T0", 300.0);
		CommandLineDouble(dU0, "U0", 20.0);
		CommandLineBool(fNoRotation, "noomega");
		CommandLineStringD(strMountainType, "mountaintype",
			"None", "(None | Wave6)");

		ParseCommandLine(argc, argv);
	EndTempestCommandLine(argv)

	// Setup the Model
	AnnounceBanner("MODEL SETUP");

	Model model(EquationSet::PrimitiveNonhydrostaticEquations);

	TempestSetupCubedSphereModel(model);

	// Set the test case for the model
	AnnounceStartBlock("Initializing test case");

	MountainWaveSphereTest::MountainType eMountainType;
	STLStringHelper::ToLower(strMountainType);
	if (strMountainType == "none") {
		eMountainType = MountainWaveSphereTest::MountainType_None;
	} else if (strMountainType == "wave6") {
		eMountainType = MountainWaveSphereTest::MountainType_Wave6;
	} else {
		_EXCEPTIONT("Invalid mountain type:"
			" Expected \"None\" or \"Wave6\"");
	}

	model.SetTestCase(
		new MountainWaveSphereTest(
			dZtop,
			dEarthScaling,
			dT0,
			dU0,
			fNoRotation,
			eMountainType));

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
