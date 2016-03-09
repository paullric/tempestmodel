///////////////////////////////////////////////////////////////////////////////
///
///	\file    ThermalBubbleCartesianTest.cpp
///	\author  Paul Ullrich, Jorge Guerra
///	\version December 18, 2013
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
///		Thermal rising bubble test case.
///	</summary>
class ThermalBubbleCartesianTest : public TestCase {

public:
        /// <summary>
	///		Lateral BC array (FOR CARTESIAN GRIDS).
	///	</summary>
	int m_iLatBC[4];

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
	///		Reference constant background pontential temperature
	///	</summary>
	double m_dThetaBar;

	///	<summary>
	///		Parameter factor for temperature disturbance
	///	</summary>
	double m_dThetaC;

	///	<summary>
	///		Parameter reference bubble radius
	///	</summary>
	double m_drC;

	///	<summary>
	///		Parameter reference length x for temperature disturbance
	///	</summary>
	double m_dxC;

	///	<summary>
	///		Parameter reference length z for temperature disturbance
	///	</summary>
	double m_dzC;

	///	<summary>
	///		Parameter Archimede's Constant (essentially Pi but to some digits)
	///	</summary>
	double m_dpiC;

public:
	///	<summary>
	///		Constructor. (with physical constants defined privately here)
	///	</summary>
	ThermalBubbleCartesianTest(
		double dThetaBar,
		double dThetaC,
		double drC,
		double dxC,
		double dzC,
		double dpiC
	) :
		m_dThetaBar(dThetaBar),
		m_dThetaC(dThetaC),
		m_drC(drC),
		m_dxC(dxC),
		m_dzC(dzC),
		m_dpiC(dpiC)
	{
		// Set the dimensions of the box
		m_dGDim[0] = 0.0;
		m_dGDim[1] = 1000.0;
		m_dGDim[2] = -500.0;
		m_dGDim[3] = 500.0;
		m_dGDim[4] = 0.0;
		m_dGDim[5] = 1000.0;

		// Set the boundary conditions for this test (no-flux in Y)
		m_iLatBC[0] = Grid::BoundaryCondition_Periodic;
		m_iLatBC[1] = Grid::BoundaryCondition_Periodic;
		m_iLatBC[2] = Grid::BoundaryCondition_Periodic;
		m_iLatBC[3] = Grid::BoundaryCondition_Periodic;
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
		return true;
	}

	///	<summary>
	///		Obtain test case specific physical constants.
	///	</summary>
	virtual void EvaluatePhysicalConstants(
		PhysicalConstants & phys
	) const {
		// No Coriolis
		phys.SetOmega(0.0);
	}

	///	<summary>
	///		Evaluate the topography at the given point. (cartesian version)
	///	</summary>
	virtual double EvaluateTopography(
	   double dXp,
	   double dYp
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

		// Potential temperature perturbation bubble using radius
		double xL2 = (dXp - m_dxC) * (dXp - m_dxC);
		double zL2 = (dZp - m_dzC) * (dZp - m_dzC);
		double dRp = sqrt(xL2 + zL2);

		double dThetaHat = 1.0;
		if (dRp <= m_drC) {
			dThetaHat = 0.5 * m_dThetaC * (1.0 + cos(m_dpiC * dRp / m_drC));
		} else if (dRp > m_drC) {
			dThetaHat = 0.0;
		}

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
	    const double dG = phys.GetG();
		const double dCv = phys.GetCv();
		const double dCp = phys.GetCp();
		const double dRd = phys.GetR();
		const double dP0 = phys.GetP0();
		// Set the uniform U, V, W field for all time
		dState[0] = 0.0;
		dState[1] = 0.0;
		dState[3] = 0.0;

		// Set the initial potential temperature field
		dState[2] = m_dThetaBar;

		// Set the initial density based on the Exner pressure
		double dExnerP =
			- dG / (dCp * m_dThetaBar) * dZp + 1.0;
		double dRho =
			dP0 / (dRd * m_dThetaBar) *
			  pow(dExnerP, (dCv / dRd));

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
	    const double dG = phys.GetG();
		const double dCv = phys.GetCv();
		const double dCp = phys.GetCp();
		const double dRd = phys.GetR();
		const double dP0 = phys.GetP0();
		// Set the uniform U, V, W field for all time
		dState[0] = 0.0;
		dState[1] = 0.0;
		dState[3] = 0.0;

		// Set the initial potential temperature field
		dState[2] = m_dThetaBar + EvaluateTPrime(phys, dXp, dZp);

		// Set the initial density based on the Exner pressure
		double dExnerP =
			- dG / (dCp * m_dThetaBar) * dZp + 1.0;
		double dRho =
			dP0 / (dRd * m_dThetaBar) *
			  pow(dExnerP, (dCv / dRd));

		dState[4] = dRho;
	}
};

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	// Initialize Tempest
	TempestInitialize(&argc, &argv);

try {
	// Reference constant background pontential temperature
	double dThetaBar;

	// Parameter factor for temperature disturbance
	double dThetaC;

	// Parameter reference bubble radius
	double drC;

	// Parameter reference length x for temperature disturbance
	double dxC;

	// Parameter reference length z for temperature disturbance
	double dzC;

	// Parameter Archimede's Constant (essentially Pi but to some digits)
	double dpiC;

	// Parse the command line
	BeginTempestCommandLine("ThermalBubbleCartesianTest");
		SetDefaultResolutionX(36);
		SetDefaultResolutionY(1);
		SetDefaultLevels(72);
		SetDefaultOutputDeltaT("10s");
		SetDefaultDeltaT("10000u");
		SetDefaultEndTime("700s");
		SetDefaultHorizontalOrder(4);
		SetDefaultVerticalOrder(1);

		CommandLineDouble(dThetaBar, "ThetaBar", 300.0);
		CommandLineDouble(dThetaC, "ThetaC", 0.5);
		CommandLineDouble(drC, "rC", 250.0);
		CommandLineDouble(dxC, "xC", 500.0);
		CommandLineDouble(dzC, "zC", 350.0);
		CommandLineDouble(dpiC, "piC", 3.14159265);

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	// Create a new instance of the test
	ThermalBubbleCartesianTest * test =
		new ThermalBubbleCartesianTest(
			dThetaBar,
			dThetaC,
			drC,
			dxC,
			dzC,
			dpiC);

	// Setup the Model
	AnnounceBanner("MODEL SETUP");

	Model model(EquationSet::PrimitiveNonhydrostaticEquations);

	// Setup the cartesian model with dimensions and reference latitude
	TempestSetupCartesianModel(model, test->m_dGDim, 0.0, 
								test->m_iLatBC,  true);

	// Set the reference length to reduce diffusion relative to global scale
	const double XL = std::abs(test->m_dGDim[1] - test->m_dGDim[0]);
	const double oneDegScale = 110000.0;
	if (XL < oneDegScale) {
		model.GetGrid()->SetReferenceLength(XL);
	}
	else {
		model.GetGrid()->SetReferenceLength(oneDegScale);
	}

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

