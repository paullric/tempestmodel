///////////////////////////////////////////////////////////////////////////////
///
///	\file    HeldSuarezTest.cpp
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

#include "HeldSuarezPhysics.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Ullrich, Melvin, Jablonowski and Staniforth (2013) Baroclinic wave test
///	</summary>
class HeldSuarezTest : public TestCase {

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
	///		Initial isothermal temperature (K)
	///	</summary>
	const double ParamT0;

protected:
	///	<summary>
	///		Model cap.
	///	</summary>
	double m_dZtop;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	HeldSuarezTest(
		double dZtop
	) :
		ParamEarthRadiusScaling(1.0),
		ParamT0(300.0),

		m_dZtop(dZtop)
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
		double dH = phys.GetR() * ParamT0 / phys.GetG();

		double dP = phys.GetP0() * exp(- dZ / dH);

		double dRho = dP / phys.GetG() / dH;

		// Store the state
		dState[0] = 0.0;
		dState[1] = 0.0;
		dState[2] = phys.RhoThetaFromPressure(dP) / dRho;
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
		double dH = phys.GetR() * ParamT0 / phys.GetG();

		double dP = phys.GetP0() * exp(- dZ / dH);

		double dRho = dP / phys.GetG() / dH;

		// Store the state
		dState[0] = 0.0;
		dState[1] = 0.0;
		dState[2] = phys.RhoThetaFromPressure(dP) / dRho;
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

	// Parse the command line
	BeginTempestCommandLine("HeldSuarezTest");
		SetDefaultResolution(20);
		SetDefaultLevels(10);
		SetDefaultOutputDeltaT("200s");
		SetDefaultDeltaT("200s");
		SetDefaultEndTime("200s");
		SetDefaultHorizontalOrder(4);
		SetDefaultVerticalOrder(1);

		CommandLineDouble(dZtop, "ztop", 30000.0);

		ParseCommandLine(argc, argv);
	EndTempestCommandLine(argv)

	// Setup the Model
	AnnounceBanner("MODEL SETUP");

	Model model(EquationSet::PrimitiveNonhydrostaticEquations);

	TempestSetupCubedSphereModel(model);

	// Set the test case for the model
	AnnounceStartBlock("Initializing test case");

	model.SetTestCase(new HeldSuarezTest(dZtop));

	AnnounceEndBlock("Done");

	// Add Held-Suarez physics
	model.AttachWorkflowProcess(
		new HeldSuarezPhysics(
			model,
			model.GetDeltaT()));

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
