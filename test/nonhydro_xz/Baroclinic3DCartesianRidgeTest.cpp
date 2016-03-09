///////////////////////////////////////////////////////////////////////////////
///
///	\file    Baroclinic3DCartesianRidgeTest.cpp
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
#include "iomanip"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		MODIFIED Ullrich and Jablonowski (2012)
///
///		Baroclinic Wave in a 3D Channel with ridge
///	</summary>
class Baroclinic3DCartesianRidgeTest : public TestCase {

public:
        /// <summary>
	///		Lateral BC array (FOR CARTESIAN GRIDS).
	///	</summary>
	int m_iLatBC[4];

	/// <summary>
	///		Grid dimension array (FOR CARTESIAN GRIDS).
	///	</summary>
	double m_dGDim[6];

	/// <summary>
	///		Reference latitude for "large" domains
	///	</summary>
	double m_dRefLat;

	///	<summary>
	///		Parameter reference height for topography disturbance
	///	</summary>
	double m_dhC;

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
	///		Assumed lapse rate of absolute temperature
	///	</summary>
	double m_ddTdz;
	
	///	<summary>
	///		Reference constant surface absolute temperature
	///	</summary>
	double m_dT0;

	///	<summary>
	///		Parameter reference length a for temperature disturbance
	///	</summary>
	double m_daC;

	///	<summary>
	///		Parameter reference length x for temperature disturbance
	///	</summary>
	double m_dXC;

	///	<summary>
	///		Parameter reference length y for temperature disturbance
	///	</summary>
	double m_dYC;

	///	<summary>
	///		Parameter for the center of the y domain
	///	</summary>
	double m_dY0;

	///	<summary>
	///		Parameter reference width for perturtion gaussian
	///	</summary>
	double m_dLpC;

	///	<summary>
	///		Parameter Archimede's Constant (essentially Pi but to some digits)
	///	</summary>
	double m_dpiC;

	///	<summary>
	///		Flag indicating that Rayleigh friction is inactive.
	///	</summary>
	bool m_fNoRayleighFriction;

public:
	///	<summary>
	///		Constructor. (with physical constants defined privately here)
	///	</summary>
	Baroclinic3DCartesianRidgeTest(
		double dbC,
		double dU0,
		double dUp,
		double ddTdz,
		double dT0,
		double dLpC,
		double dhC,
		double daC,
		double dXC,
		double dYC,
		bool fNoRayleighFriction
	) :
		m_dbC(dbC),
		m_dU0(dU0),
		m_dUp(dUp),
		m_ddTdz(ddTdz),
		m_dT0(dT0),
		m_dLpC(dLpC),
		m_dhC(dhC),
		m_daC(daC),
		m_dXC(dXC),
		m_dYC(dYC),
		m_fNoRayleighFriction(fNoRayleighFriction)
	{
		m_dpiC = M_PI;

		// Set the dimensions of the box
		m_dGDim[0] = 0.0;
		m_dGDim[1] = 30000000.0;
		m_dGDim[2] = 0.0;
		m_dGDim[3] = 6000000.0;
		m_dGDim[4] = 0.0;
		m_dGDim[5] = 30000.0;

		// Set the reference latitude
		m_dRefLat = 45.0 / 180.0 * m_dpiC;

		// Set the center of the domain in Y
		m_dY0 = 0.5 * (m_dGDim[3] - m_dGDim[2]);

		// Set the boundary conditions for this test (no-flux in Y)
		m_iLatBC[0] = Grid::BoundaryCondition_Periodic;
		m_iLatBC[1] = Grid::BoundaryCondition_NoFlux;
		m_iLatBC[2] = Grid::BoundaryCondition_Periodic;
		m_iLatBC[3] = Grid::BoundaryCondition_NoFlux;
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
		const PhysicalConstants & phys,
		double dXp,
		double dYp
	) const {
		// Specify the ridge to a factor of 5 downstream from the perturbation
		double xLoc = 2.0 * m_dXC;
		double hsm = m_dhC / (1.0 + exp(((dXp - xLoc)/m_daC) *
					((dXp - xLoc)/m_daC)));
		//double yFac = 1.0 / (1.0 + exp(pow((dYp - m_dY0)/(m_daC),4)));
		double yFac = 1.0;
		hsm *= yFac;
		//std::cout << hsm << ' ' << xLoc << "\n";
		return hsm;
	}

	///	<summary>
	///		Flag indicating whether or not Rayleigh friction strength is given.
	///	</summary>
	virtual bool HasRayleighFriction() const {
		return !m_fNoRayleighFriction;
	}

	///	<summary>
	///		Evaluate the Rayleigh friction strength at the given point.
	///	</summary>
	virtual double EvaluateRayleighStrength(
		double dZ,
		double dXp,
		double dYp
	) const {
		const double dRayleighStrength = 1.0E-3;
		const double dRayleighDepth = 5000.0;
		const double dRayleighWidth = 2.0E6;

		double dNuDepth = 0.0;
		double dNuRight = 0.0;
		double dNuLeft  = 0.0;

		if (dZ > m_dGDim[5] - dRayleighDepth) {
			double dNormZ = (m_dGDim[5] - dZ) / dRayleighDepth;
			dNuDepth = 0.5 * dRayleighStrength * (1.0 + cos(M_PI * dNormZ));
			//dNuDepth = 0.0;
		}
		if (dXp > m_dGDim[1] - dRayleighWidth) {
			double dNormX = (m_dGDim[1] - dXp) / dRayleighWidth;
			dNuRight = 0.5 * dRayleighStrength * (1.0 + cos(M_PI * dNormX));
		}
		if (dXp < m_dGDim[0] + dRayleighWidth) {
			double dNormX = (dXp - m_dGDim[0]) / dRayleighWidth;
			dNuLeft = 0.5 * dRayleighStrength * (1.0 + cos(M_PI * dNormX));
		}

		if ((dNuDepth >= dNuRight) && (dNuDepth >= dNuLeft)) {
			return dNuDepth;
		}
		if (dNuRight >= dNuLeft) {
			return dNuRight;
		}
		return dNuLeft;
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
		double xL2 = (dXp - m_dXC) * (dXp - m_dXC);
		double yL2 = (dYp - m_dYC) * (dYp - m_dYC);

		double dUpert = m_dUp * exp(-(xL2 + yL2) / (m_dLpC * m_dLpC));
 
		return dUpert;
	}

	///	<summary>
	///		Calculate the geopotential and temperature at the given point.
	///	</summary>
	void CalculateGeopotentialTemperature(
		const PhysicalConstants & phys,
		double dEta,
		double dXp,
		double dYp,
		double dZp,
		double & dGeopotential,
		double & dTemperature
	) const {
		// Get some constants
		const double dG = phys.GetG();
		const double dCv = phys.GetCv();
		const double dCp = phys.GetCp();
		const double dRd = phys.GetR();
		const double dP0 = phys.GetP0();
		const double dae = phys.GetEarthRadius();
		const double df0 = 2 * phys.GetOmega() * sin(m_dRefLat);
		//const double df0 = 0.0;
		const double dbeta0 = 2 * phys.GetOmega() * cos(m_dRefLat) / (dae + dZp);
		//const double dbeta0 = 0.0;
		const double dLy = m_dGDim[3] - m_dGDim[2];

		// Horizontally averaged temperature
		double dAvgTemperature =
			m_dT0 * pow(dEta, dRd * m_ddTdz / dG);

		// Horizontally averaged geopotential
		double dAvgGeopotential =
			m_dT0 * dG / m_ddTdz * 
			(1.0 - pow(dEta, dRd * m_ddTdz / dG));

		// Horizontal variation geopotential function
		double dXYGeopotential = 0.5 * m_dU0 * 
			((df0 - dbeta0 * m_dY0) * (dYp - m_dY0 - 
			m_dY0 / m_dpiC * sin(2 * m_dpiC * dYp / dLy)) + 
			0.5 * dbeta0 * 
			(dYp * dYp - dLy * dYp / m_dpiC * sin(2 * m_dpiC * dYp / dLy) - 
			0.5 * dLy * dLy / (m_dpiC * m_dpiC) * cos(2 * m_dpiC * dYp / dLy) -
			dLy * dLy / 3 - 0.5 * dLy * dLy / (m_dpiC * m_dpiC)));

		double dExpDecay = exp(-(log(dEta) / m_dbC) * (log(dEta) / m_dbC));
		double dRefProfile1 = log(dEta);
		double dRefProfile2 = 2 / (m_dbC * m_dbC) * log(dEta) * log(dEta) - 1.0;

		// Total geopotential distribution
		dGeopotential = dAvgGeopotential + dXYGeopotential*
			dRefProfile1 * dExpDecay;
		
		// Total temperature distribution
		dTemperature = dAvgTemperature + dXYGeopotential / dRd *
			dRefProfile2 * dExpDecay;
	}

	///	<summary>
	///		Calculate eta at the given point via Newton iteration.  The
	///		geopotential and temperature at this point are also returned via
	///		command-line parameters.
	///	</summary>
	double EtaFromRLL(
		const PhysicalConstants &phys,
		double dZp,
		double dXp,
		double dYp,
		double & dGeopotential,
		double & dTemperature
	) const {
		const int MaxIterations  = 100;
		const double InitialEta  = 1.0e-5;
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
				phys, dEta, dXp, dYp, dZp, dGeopotential, dTemperature);

			dF     = - phys.GetG() * dZp + dGeopotential;
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
		double dZp,
		double dXp,
		double dYp,
		double * dState
	) const {
		const double dLy = m_dGDim[3] - m_dGDim[2];

		// Pressure coordinate
		double dGeopotential;
		double dTemperature;

		double dEta = EtaFromRLL(
			phys, dZp, dXp, dYp, dGeopotential, dTemperature);

		// Calculate zonal velocity and set other velocity components
		double dExpDecay = exp(-(log(dEta) / m_dbC) * (log(dEta) / m_dbC));
		double dUlon =
			-m_dU0 * sin(m_dpiC * dYp / dLy) * sin(m_dpiC * dYp / dLy) *
			 log(dEta) * dExpDecay;

		dState[0] = dUlon;
		dState[1] = 0.0;
		dState[3] = 0.0;

		// Calculate rho and theta
		double dPressure = phys.GetP0() * dEta;
		//std::cout << std::setprecision(16) << "Z = " << dZp << " Eta = " << dEta << "\n";

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
		double dZp,
		double dXp,
		double dYp,
		double * dState,
		double * dTracer
	) const {

		// Evaluate the reference state at this point
		EvaluateReferenceState(phys, dZp, dXp, dYp, dState);

		// Add perturbation in zonal velocity
		dState[0] += 0.0;
		//dState[0] += EvaluateUPrime(phys, dXp, dYp);
	}
};

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	// Initialize Tempest
	TempestInitialize(&argc, &argv);

try {
	// Nondimensional vertical width parameter
	double dbC;
	
	// Uniform zonal velocity
	double dU0;

	// Magnitude of the zonal wind perturbation
	double dUp;

	// Lapse rate
	double ddTdz;

	// Reference absolute temperature
	double dT0;

	// Width parameter for the perturbation
	double dLpC;

	// Amplitude parameter of the ridge
	double daC;

	// Parameter reference height for temperature disturbance
	double dhC;

	// Center position of the perturbation
	double dXC;

	// Center position of the perturbation
	double dYC;

	// No Rayleigh friction
	bool fNoRayleighFriction;

	// Parse the command line
	BeginTempestCommandLine("Baroclinic3DCartesianRidgeTest");
		SetDefaultResolutionX(288);
		SetDefaultResolutionY(48);
		SetDefaultLevels(32);
		SetDefaultOutputDeltaT("3h");
		SetDefaultDeltaT("300s");
		SetDefaultEndTime("12d");
		SetDefaultHorizontalOrder(4);
		SetDefaultVerticalOrder(4);

		CommandLineDouble(dbC, "b", 2.0);
		CommandLineDouble(dU0, "u0", 35.0);
		CommandLineDouble(dUp, "up", 1.0);
		CommandLineDouble(ddTdz, "gamma", 0.005);
		CommandLineDouble(dT0, "T0", 288.0);
		CommandLineDouble(dLpC, "Lp", 600000.0);
		CommandLineDouble(dhC, "hC", 400.0);
		CommandLineDouble(daC, "aC", 1000000.0);
		CommandLineDouble(dXC, "Xc", 2000000.0);
		CommandLineDouble(dYC, "Yc", 3000000.0);
		CommandLineBool(fNoRayleighFriction, "norayleigh");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	// Create a new instance of the test
	Baroclinic3DCartesianRidgeTest * test =
		new Baroclinic3DCartesianRidgeTest(dbC,
			dU0,
			dUp,
			ddTdz,
			dT0,
			dLpC,
			dhC,
			daC,
			dXC,
			dYC,
			fNoRayleighFriction);

	// Setup the Model
	AnnounceBanner("MODEL SETUP");

	Model model(EquationSet::PrimitiveNonhydrostaticEquations);
	
	// Setup the cartesian model with dimensions and reference latitude
	TempestSetupCartesianModel(model, test->m_dGDim, test->m_dRefLat, 
								test->m_iLatBC, false);

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

