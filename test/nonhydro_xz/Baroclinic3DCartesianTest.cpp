///////////////////////////////////////////////////////////////////////////////
///
///	\file    Baroclinic3DCartesianTest.cpp
///	\author  Paul Ullrich, Jorge Guerra
///	\version January 13, 2015
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
class Baroclinic3DCartesianTest : public TestCase {

public:
	/// <summary>
	///		Grid dimension array (FOR CARTESIAN GRIDS).
	///	</summary>
	double m_dGDim[6];

	/// <summary>
	///		Reference latitude for "large" domains
	///	</summary>
	double m_dRefLat;

private:

	///	<summary>
	///		Nondimensional vertical width parameter
	///	</summary>
	double m_dbC;

	///	<summary>
	///		Reference zonal U velocity (balanced jet).
	///	</summary>
	double m_dU0;

	///	<summary>
	///		Reference zonal wind perturbation.
	///	</summary>
	double m_dUp;

	///	<summary>
	///		Brunt-Vaisala frequency
	///	</summary>
	double m_dNbar;
	///	<summary>
	///		Reference constant background pontential temperature
	///	</summary>
	double m_dTheta0;

	///	<summary>
	///		Parameter reference length x for temperature disturbance
	///	</summary>
	double m_dxC;

	///	<summary>
	///		Parameter reference length z for temperature disturbance
	///	</summary>
	double m_dyC;

	///	<summary>
	///		Parameter reference width for perturtion gaussian
	///	</summary>
	double m_dLpC;

	///	<summary>
	///		Parameter Archimede's Constant (essentially Pi but to some digits)
	///	</summary>
	double m_dpiC;

public:
	///	<summary>
	///		Constructor. (with physical constants defined privately here)
	///	</summary>
	Baroclinic3DCartesianTest() :
		m_dbC(2.),
		m_dU0(35.),
		m_dUp(1.),
		m_dNbar(0.014),
		m_dTheta0(288.),
		m_dLpC(600000.),
		m_dxC(2000000.),
		m_dyC(2500000.),
		m_dpiC(3.14159265)
	{
		// Set the dimensions of the box
		m_dGDim[0] = 0.0;
		m_dGDim[1] = 40000000.0;
		m_dGDim[2] = 0.0;
		m_dGDim[3] = 6000000.0;
		m_dGDim[4] = 0.0;
		m_dGDim[5] = 30000.0;

		// Set the reference latitude
		m_dRefLat = 45.0 / 180.0 * m_dpiC;
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
		// Do nothing to the PhysicalConstants for global simulations
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
	///		Evaluate the zonal velocity field perturbation.
	///	</summary>
	double EvaluateUPrime(
		const PhysicalConstants & phys,
		double dXp,
		double dYp
	) const {

		// Gaussian perturbation for the zonal jet
		double xL2 = (dXp - m_dxC) * (dXp - m_dxC);
		double yL2 = (dYp - m_dyC) * (dYp - m_dyC);

		double dUpert = m_dUp * exp(-(xL2 + yL2) / (m_dLpC * m_dLpC));
 
		return dUpert;
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
		double dThetaBar = m_dTheta0 * exp(m_dNbar * m_dNbar / dG * dZp);

		// Set the uniform V, W field for all time
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
		double dEta = pow((dRho * dRd * dThetaBar / dP0),(dCp / dCv));
	
		// Set the balanced zonal jet from Ullrich, 2014
		double dUJet = -m_dU0 * sin(m_dpiC * dYp / m_dGDim[3]) * 
					sin(m_dpiC * dYp / m_dGDim[3]) * 
					log(dEta) * exp(-(log(dEta) / m_dbC) * (log(dEta) / m_dbC));
		dState[0] = dUJet;

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
		double dThetaBar = m_dTheta0 * exp(m_dNbar * m_dNbar / dG * dZp);

		// Set the uniform V, W field for all time
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
		double dEta = pow((dRho * dRd * dThetaBar / dP0),(dCp / dCv));
	
		// Set the balanced zonal jet from Ullrich, 2014
		double dUJet = -m_dU0 * sin(m_dpiC * dYp / m_dGDim[3]) * 
					sin(m_dpiC * dYp / m_dGDim[3]) * 
					log(dEta) * exp(-(log(dEta) / m_dbC) * (log(dEta) / m_dbC));
		dState[0] = dUJet  + EvaluateUPrime(phys, dXp, dYp);

		dState[4] = dRho;
	}
};

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	// Initialize Tempest
	TempestInitialize(&argc, &argv);

try {

	// Parse the command line
	BeginTempestCommandLine("Baroclinic3DCartesianTest");
		SetDefaultResolutionX(288);
		SetDefaultResolutionY(48);
		SetDefaultLevels(32);
		SetDefaultOutputDeltaT("3h");
		SetDefaultDeltaT("300s");
		SetDefaultEndTime("12d");
		SetDefaultHorizontalOrder(4);
		SetDefaultVerticalOrder(1);

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	// Create a new instance of the test
	Baroclinic3DCartesianTest * test =
		new Baroclinic3DCartesianTest();

	// Setup the Model
	AnnounceBanner("MODEL SETUP");

	Model model(EquationSet::PrimitiveNonhydrostaticEquations);
	
	// Setup the cartesian model with dimensions and reference latitude
	TempestSetupCartesianModel(model, test->m_dGDim, test->m_dRefLat);

	// Set the reference length to reduce diffusion (1100km)
	model.GetGrid()->SetReferenceLength(1100000.0);

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

