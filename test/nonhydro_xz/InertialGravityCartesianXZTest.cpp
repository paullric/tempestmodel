///////////////////////////////////////////////////////////////////////////////
///
///	\file    InertialGravityCartesianXZTest.cpp
///	\author  Paul Ullrich, Jorge Guerra
///	\version October 2, 2013
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
///		Giraldo et al. (2007)
///
///		Intertia-gravity waves test case.
///	</summary>
class InertialGravityCartesianXZTest : public TestCase {

public:
	/// <summary>
	///		Grid dimension array (FOR CARTESIAN GRIDS).
	///	</summary>
	double m_dGDim[6];

private:
	///	<summary>
	///		Background height field.
	///	</summary>
	double m_dH0;

	///	<summary>
	///		Uniform +X flow field.
	///	</summary>
	double m_dU0;

	///	<summary>
	///		Brunt-Vaisala frequency
	///	</summary>
	double m_dNbar;

	///	<summary>
	///		Reference pontential temperature
	///	</summary>
	double m_dTheta0;

	///	<summary>
	///		Parameter factor for temperature disturbance
	///	</summary>
	double m_dThetaC;

	///	<summary>
	///		Parameter reference height for temperature disturbance
	///	</summary>
	double m_dhC;

	///	<summary>
	///		Parameter reference length a for temperature disturbance
	///	</summary>
	double m_daC;

	///	<summary>
	///		Parameter reference length x for temperature disturbance
	///	</summary>
	double m_dxC;

	///	<summary>
	///		Parameter Archimede's Constant (essentially Pi but to some digits)
	///	</summary>
	double m_dpiC;

	///	<summary>
	///		Flag indicating that the reference profile should be used.
	///	</summary>
	bool m_fNoReferenceState;

public:
	///	<summary>
	///		Constructor. (with physical constants defined privately here)
	///	</summary>
	InertialGravityCartesianXZTest(
		bool fNoReferenceState
	) :
		m_dH0(10000.),
		m_dU0(20.),
		m_dNbar(0.01),
		m_dTheta0(300.0),
		m_dThetaC(1.0),
		m_dhC(10000.),
		m_daC(5000.),
		m_dxC(1.0E+5),
		m_dpiC(3.14159265),
		m_fNoReferenceState(fNoReferenceState)
	{
		// Set the dimensions of the box
		m_dGDim[0] = 0.0;
		m_dGDim[1] = 300000.0;
		m_dGDim[2] = -100000.0;
		m_dGDim[3] = 100000.0;
		m_dGDim[4] = 0.0;
		m_dGDim[5] = 10000.0;
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
		return m_dGDim[5];
	}

	///	<summary>
	///		Flag indicating that a reference state is available.
	///	</summary>
	virtual bool HasReferenceState() const {
		return !m_fNoReferenceState;
	}

	///	<summary>
	///		Obtain test case specific physical constants.
	///	</summary>
	virtual void EvaluatePhysicalConstants(
		PhysicalConstants & phys
	) const {
		// Do nothing to the PhysicalConstants for global simulations
	}

	///	<summary>
	///		Evaluate the topography at the given point. (cartesian version)
	///	</summary>
	virtual double EvaluateTopography(
	   double dxp,
	   double dyp
	) const {
		// This test case has no topography associated with it
		return 0.0;
	}

	///	<summary>
	///		Evaluate the potential temperature field perturbation.
	///	</summary>
	double EvaluateTPrime(
		const PhysicalConstants & phys,
		double dXp,
		double dZp
	) const {
		double dG = phys.GetG();

		// Potential temperature perturbation
		double dThetaHat1 = m_dThetaC * sin(m_dpiC * dZp / m_dhC);
		double argX = (dXp - m_dxC)/m_daC;
		double dThetaHat2 = (1.0 + argX * argX);
		double dThetaHat = dThetaHat1 / dThetaHat2;

		return dThetaHat;
	}

	///	<summary>
	///		Evaluate the reference state at the given point.
	///	</summary>
	virtual void EvaluateReferenceState(
		const PhysicalConstants & phys,
		double dZp,
		double dXp,
		double dYp,
		double * dState
	) const {
		// Base potential temperature field
		const double dG = phys.GetG();
		const double dCv = phys.GetCv();
		const double dCp = phys.GetCp();
		const double dRd = phys.GetR();
		const double dP0 = phys.GetP0();
		double dThetaBar = m_dTheta0 * exp(m_dNbar * m_dNbar / dG * dZp);

		// Set the uniform U, V, W field for all time
		dState[0] = m_dU0;
		dState[1] = 0.0;
		dState[3] = 0.0;

		// Set the initial potential temperature field
		dState[2] = dThetaBar;

		// Set the initial density based on the Exner pressure
		double dExnerP = (dG * dG) / (dCp * m_dTheta0 * (m_dNbar * m_dNbar));
		dExnerP *= (exp(-pow(m_dNbar,2.0)/dG * dZp) - 1.0);
		dExnerP += 1.0;
		double dRho = dP0 / (dRd * dThetaBar) * pow(dExnerP,(dCv / dRd));
		dState[4] = dRho;
	}

	///	<summary>
	///		Evaluate the state vector at the given point.
	///	</summary>
	virtual void EvaluatePointwiseState(
		const PhysicalConstants & phys,
		const Time & time,
		double dZp,
		double dXp,
		double dYp,
		double * dState,
		double * dTracer
	) const {
		// Base potential temperature field
		const double dG = phys.GetG();
		const double dCv = phys.GetCv();
		const double dCp = phys.GetCp();
		const double dRd = phys.GetR();
		const double dP0 = phys.GetP0();
		double dThetaBar = m_dTheta0 * exp(m_dNbar * m_dNbar / dG * dZp);

		// Set the uniform U, V, W field for all time
		dState[0] = m_dU0;
		dState[1] = 0.0;
		dState[3] = 0.0;

		// Set the initial potential temperature field
		dState[2] = dThetaBar + EvaluateTPrime(phys, dXp, dZp);

		// Set the initial density based on the Exner pressure
		double dExnerP =
			(dG * dG) / (dCp * m_dTheta0 * (m_dNbar * m_dNbar));
		dExnerP *= (exp(-(m_dNbar * m_dNbar)/dG * dZp) - 1.0);
		dExnerP += 1.0;

		double dRho =
			dP0 / (dRd * dThetaBar) * pow(dExnerP, (dCv / dRd));
		dState[4] = dRho;
	}
};

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	// Initialize Tempest
	TempestInitialize(&argc, &argv);

try {
	// No Rayleigh friction
	bool fNoRayleighFriction;

	// Parse the command line
	BeginTempestCommandLine("InertialGravityCartesianXZTest");
		SetDefaultResolutionX(40);
		SetDefaultResolutionY(1);
		SetDefaultLevels(48);
		SetDefaultOutputDeltaT("250s");
		SetDefaultDeltaT("500000u");
		SetDefaultEndTime("3000s");
		SetDefaultHorizontalOrder(4);
		SetDefaultVerticalOrder(3);

		CommandLineBool(fNoRayleighFriction, "norayleigh");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	// Create a new instance of the test
	InertialGravityCartesianXZTest * test =
		new InertialGravityCartesianXZTest(fNoRayleighFriction);

	// Setup the Model
	AnnounceBanner("MODEL SETUP");

	Model model(EquationSet::PrimitiveNonhydrostaticEquations);

	TempestSetupCartesianModel(model, test->m_dGDim, 0.0);

	// Set the test case for the model
	AnnounceStartBlock("Initializing test case");
	model.SetTestCase(test);
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

