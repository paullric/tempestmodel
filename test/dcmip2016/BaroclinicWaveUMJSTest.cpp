///////////////////////////////////////////////////////////////////////////////
///
///	\file	BaroclinicWaveUMJSTest.cpp
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

#include "DCMIPPhysics.h"
#include "TerminatorPhysics.h"

extern "C" {
	void baroclinic_wave_test(
		int * iDeep,
		int * iMoist,
		int * iPertT,
		double * dX,
		double * dLon,
		double * dLat,
		double * dP,
		double * dZ,
		int * iZCoords,
		double * dU,
		double * dV,
		double * dT,
		double * dThetaV,
		double * dPhis,
		double * dPs,
		double * dRho,
		double * dQ);

	void initial_value_Terminator(
		double * lat,
		double * lon,
		double * cl,
		double * cl2);
}



///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Ullrich, Melvin, Jablonowski and Staniforth (2013) Baroclinic wave test
///	</summary>
class BaroclinicWaveUMJSTest : public TestCase {

protected:
	///	<summary>
	///		Model cap.
	///	</summary>
	double m_dZtop;

	///	<summary>
	///		Earth radius scaling parameter.
	///	</summary>
	double m_dEarthScaling;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	BaroclinicWaveUMJSTest(
		double dZtop,
		double dEarthScaling
	) :
		m_dZtop(dZtop),
		m_dEarthScaling(dEarthScaling)
	{ }

public:
	///	<summary>
	///		Number of tracers used in this test.
	///	</summary>
	virtual int GetTracerCount() const {
		return 3;
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
		double dZ,
		double dLon,
		double dLat,
		double * dState
	) const {
		dState[0] = 0.0;
		dState[1] = 0.0;
		dState[2] = 0.0;
		dState[3] = 0.0;
		dState[4] = 0.0;
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
		int iDeep = 0;
		int iMoist = 1;
		int iPertT = 0;
		int iZCoords = 1;
		double dX = 1.0;

		double dRho;
		double dU;
		double dV;
		double dP;
		double dT;
		double dThetaV;
		double dPhis;
		double dPs;
		double dQ;
		
		baroclinic_wave_test(
			&iDeep,
			&iMoist,
			&iPertT,
			&dX,
			&dLon,
			&dLat,
			&dP,
			&dZ,
			&iZCoords,
			&dU,
			&dV,
			&dT,
			&dThetaV,
			&dPhis,
			&dPs,
			&dRho,
			&dQ);

		// Terminator test
		double dQCl = 0.0;
		double dQCl2 = 0.0;

		double dLatDeg = dLat * 180.0 / M_PI;
		double dLonDeg = dLon * 180.0 / M_PI;

		initial_value_Terminator(
			&dLatDeg,
			&dLonDeg,
			&dQCl,
			&dQCl2);

		// Store the state
		dState[0] = dU;
		dState[1] = dV;
		dState[2] = dThetaV;
		dState[3] = 0.0;
		dState[4] = dRho;
		
		dTracer[0] = dQ * dRho;
		dTracer[1] = 0;
		dTracer[2] = 0;

		dTracer[3] = dQCl * dRho;
		dTracer[4] = dQCl2 * dRho;
	}
};

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	// Initialize Tempest
	TempestInitialize(&argc, &argv);

try {
	// Model cap
	double dZtop;

	// Earth radius scaling parameter
	double dEarthScaling;

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
		CommandLineDouble(dEarthScaling, "X", 1.0);

		ParseCommandLine(argc, argv);
	EndTempestCommandLine(argv)

	// Setup the Model
	AnnounceBanner("MODEL SETUP");

	EquationSet eqn(EquationSet::PrimitiveNonhydrostaticEquations);
 
	eqn.InsertTracer("RhoQv", "RhoQv");
	eqn.InsertTracer("RhoQc", "RhoQc");
	eqn.InsertTracer("RhoQr", "RhoQr");

	eqn.InsertTracer("RhoQCl", "RhoQCl");
	eqn.InsertTracer("RhoQCl2", "RhoQCl2");

	UserDataMeta metaUserData;

	metaUserData.InsertDataItem2D("PRECT");

	Model model(eqn, metaUserData);

	TempestSetupCubedSphereModel(model);

	// Set the test case for the model
	AnnounceStartBlock("Initializing test case");

	model.SetTestCase(
		new BaroclinicWaveUMJSTest(
			dZtop,
			dEarthScaling));
	AnnounceEndBlock("Done");

	// Add DCMIP physics
	model.AttachWorkflowProcess(
		new DCMIPPhysics(
			model,
			model.GetDeltaT(),
			1,
			0,
			0));

	// Add Terminator physics
	model.AttachWorkflowProcess(
		new TerminatorPhysics(
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

