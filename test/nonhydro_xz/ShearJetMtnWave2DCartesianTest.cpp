///////////////////////////////////////////////////////////////////////////////
///
///	\file    ShearJetMtnWave2DCartesianTest.cpp
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
///		Giraldo et al. (2007)
///
///		Thermal rising bubble test case.
///	</summary>
class ShearJetMtnWave2DCartesianTest : public TestCase {

public:
	/// <summary>
	///		Lateral BC array (FOR CARTESIAN GRIDS).
	///	</summary>
	int m_iLatBC[4];

	/// <summary>
	///		Grid dimension array (FOR CARTESIAN GRIDS).
	///	</summary>
	double m_dGDim[6];

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
	double m_dUj;

	///	<summary>
	///		Assumed lapse rate of absolute temperature
	///	</summary>
	double m_ddTdz;

	///	<summary>
	///		Assumed lapse rate of absolute temperature (stratosphere)
	///	</summary>
	double m_ddTdzSTR;
	
	///	<summary>
	///		Reference constant surface absolute temperature
	///	</summary>
	double m_dT0;

	///	<summary>
	///		Parameter reference length a for temperature disturbance
	///	</summary>
	double m_daC;

	///	<summary>
	///		Parameter reference length for mountain profile
	///	</summary>
	double m_dlC;

	///	<summary>
	///		Parameter for the center of the y domain
	///	</summary>
	double m_dY0;

	///	<summary>
	///		Parameter Archimede's Constant (essentially Pi but to some digits)
	///	</summary>
	double m_dpiC;

	///	<summary>
	///		Flag indicating that Rayleigh friction is inactive.
	///	</summary>
	bool m_fNoRayleighFriction;

	///	<summary>
	///		Elevation at the tropopause (m) USER SPECIFIED.
	///	</summary>
	double m_dTPHeight;

	///	<summary>
	///		Thickness of the mixed layer above tropopause (m) USER SPECIFIED.
	///	</summary>
	double m_dTPMixedLayerH;

	///	<summary>
	///		Geopotential at the tropopause (m) DERIVED.
	///	</summary>
	double m_dTPPhi1;

	///	<summary>
	///		Geopotential at the top of the mixed layer (m) DERIVED.
	///	</summary>
	double m_dTPPhi2;

	///	<summary>
	///		Temperature at the tropopause (K) DERIVED.
	///	</summary>
	double m_dTPTemp1;

	///	<summary>
	///		Temperature at the top of the mixed layer (K) DERIVED.
	///	</summary>
	double m_dTPTemp2;

	///	<summary>
	///		Sigma coordinate value at the tropopause DERIVED.
	///	</summary>
	double m_dTPEta1;

	///	<summary>
	///		Sigma coordinate value at the top of the mixed layer DERIVED.
	///	</summary>
	double m_dTPEta2;

	///<summary>
	///		Uniform diffusion coefficient for scalars.
	///</summary>
	double m_dUCoeffS;

	///<summary>
	///		Uniform diffusion coefficient for vectors.
	///	</summary>
	double m_dUCoeffV;

public:
	///	<summary>
	///		Constructor. (with physical constants defined privately here)
	///	</summary>
	ShearJetMtnWave2DCartesianTest(
		const PhysicalConstants & phys,
		double dbC,
		double dU0,
		double dUj,
		double ddTdz,
		double ddTdzSTR,
		double dT0,
		double dhC,
		double daC,
		double dlC,
		double dUCoeffS,
		double dUCoeffV,
		bool fNoRayleighFriction
	) :
		m_dbC(dbC),
		m_dU0(dU0),
		m_dUj(dUj),
		m_ddTdz(ddTdz),
		m_ddTdzSTR(ddTdzSTR),
		m_dT0(dT0),
		m_dhC(dhC),
		m_daC(daC),
		m_dlC(dlC),
		m_dUCoeffS(dUCoeffS),
		m_dUCoeffV(dUCoeffV),
		m_fNoRayleighFriction(fNoRayleighFriction)
	{
		m_dpiC = M_PI;

		// Set the dimensions of the box
		m_dGDim[0] = -40000.0;
		m_dGDim[1] = 40000.0;
		m_dGDim[2] = -500.0;
		m_dGDim[3] = 500.0;
		m_dGDim[4] = 0.0;
		m_dGDim[5] = 30000.0;

		// Set the center of the domain in Y
		m_dY0 = 0.5 * (m_dGDim[3] - m_dGDim[2]);

		// Set the tropopause elevation and mixed layer depth in meters
		m_dTPHeight = 12000.0;
		m_dTPMixedLayerH = 3000.0;

		// Find and set the tropopause temperature here before initializing the test
		double dGeopotential;
		double dTemperature;

		// Get the temperature at the tropopause
		double dEta = EtaFromRLL(
			phys, m_dTPHeight, 0.0, 0.0, dGeopotential, dTemperature);
		m_dTPTemp1 = dTemperature;
		m_dTPEta1 = dEta;
		m_dTPPhi1 = dGeopotential;

		// Get the pressure level at the top of the mixed layer
		dEta = EtaFromRLL(
			phys, m_dTPHeight + m_dTPMixedLayerH, 0.0, 0.0, 
			dGeopotential, dTemperature);
		m_dTPTemp2 = dTemperature;
		m_dTPEta2 = dEta;
		m_dTPPhi2 = dGeopotential;

		// Set the boundary conditions for this test
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
		// Do nothing to the PhysicalConstants for global simulations
	}

	///	<summary>
	///		Strength of the uniform diffusion (m^2/s)
	///	</summary>
	virtual void GetUniformDiffusionCoeffs(
		double & dScalarUniformDiffusionCoeff,
		double & dVectorUniformDiffusionCoeff
	) const {
		dScalarUniformDiffusionCoeff = m_dUCoeffS;
		dVectorUniformDiffusionCoeff = m_dUCoeffV;
	}
//
	///	<summary>
	///		Evaluate the topography at the given point. (cartesian version)
	///	</summary>
	virtual double EvaluateTopography(
		const PhysicalConstants & phys,
		double dXp,
		double dYp
	) const {
		// Specify the Schar Mountain (Test case 5 from Giraldo et al. 2008)
		double hsm = m_dhC * exp(-dXp/m_daC * dXp/m_daC) *
					 cos(M_PI * dXp / m_dlC) * cos(M_PI * dXp / m_dlC);
		//std::cout << hsm << "\n";

		return hsm;
	}
//
/*
	///	<summary>
	///		Evaluate the topography at the given point. (cartesian version)
	///	</summary>
	virtual double EvaluateTopography(
		const PhysicalConstants & phys,
		double dXp,
		double dYp
	) const {
		// Specify the Linear Mountain (case 6 from Giraldo et al. 2008)
		double hsm = m_dhC / (1.0 + ((dXp - 0.0)/m_daC) * ((dXp - 0.0)/m_daC) *
									((dXp - 0.0)/m_daC) * ((dXp - 0.0)/m_daC));
		//std::cout << hsm << "\n";
		return hsm;
	}
*/
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
		const double dRayleighStrengthZ = 1.0E-2;//8.0e-3;
		const double dRayleighStrengthX = 1.0 * dRayleighStrengthZ;
		const double dRayleighDepth = 5000.0;
		const double dRayleighWidth = 5000.0;

		double dNuDepth = 0.0;
		double dNuRight = 0.0;
		double dNuLeft  = 0.0;

		if (dZ > m_dGDim[5] - dRayleighDepth) {
			double dNormZ = (m_dGDim[5] - dZ) / dRayleighDepth;
			dNuDepth = 0.5 * dRayleighStrengthZ * (1.0 + cos(M_PI * dNormZ));
		}
		if (dXp > m_dGDim[1] - dRayleighWidth) {
			double dNormX = (m_dGDim[1] - dXp) / dRayleighWidth;
			dNuRight = 0.5 * dRayleighStrengthX * (1.0 + cos(M_PI * dNormX));
		}
		if (dXp < m_dGDim[0] + dRayleighWidth) {
			double dNormX = (dXp - m_dGDim[0]) / dRayleighWidth;
			dNuLeft = 0.5 * dRayleighStrengthX * (1.0 + cos(M_PI * dNormX));
		}
		
		//std::cout << dXp << ' ' << dZ << ' ' << dNuDepth << std::endl;
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
		return 0.0;
	}

	///	<summary>
	///		Calculate the geopotential and temperature at the given point.
	///	</summary>
	void CalculateGeopotentialTemperature(
		const PhysicalConstants & phys,
		double dEta,
		double dZp,
		double dXp,
		double dYp,
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
		const double df0 = 0.0;
		const double dbeta0 = 0.0;
		const double dLy = m_dGDim[3] - m_dGDim[2];

		// Horizontally averaged temperature profile (piecewise continuous)
		// Horizontally averaged geopotential
		double dAvgGeopotential = 0.0;
		double dAvgTemperature = 0.0;
		if (dZp <= m_dTPHeight) {
			dAvgTemperature = m_dT0 * pow(dEta, dRd * m_ddTdz / dG);
			dAvgGeopotential =
					m_dT0 * dG / m_ddTdz * 
					(1.0 - pow(dEta, dRd * m_ddTdz / dG));
		}
		else if ((dZp > m_dTPHeight)&&(dZp <= m_dTPHeight + m_dTPMixedLayerH)) {
			dAvgTemperature = m_dTPTemp1;
			dAvgGeopotential = -dRd * m_dTPTemp1 * log(dEta) + 
								dRd * m_dTPTemp1 * log(m_dTPEta1) + m_dTPPhi1;
		}
		else if (dZp > m_dTPHeight + m_dTPMixedLayerH) {
			dAvgTemperature = m_dTPTemp1 * 
					pow((dEta / m_dTPEta2), dRd * m_ddTdzSTR / dG);
			dAvgGeopotential =
					m_dTPTemp1 * dG / m_ddTdzSTR * 
					(1.0 - pow((dEta / m_dTPEta2), dRd * m_ddTdzSTR / dG)) 
					+ m_dTPPhi2;
		}

		// Horizontal variation geopotential function
		double dXYGeopotential = 0.0;

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
		const int MaxIterations  = 200;
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
				phys, dEta, dZp, dXp, dYp, dGeopotential, dTemperature);

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
		double dUlon = -m_dUj * 0.5 * log(dEta) * dExpDecay;

		dState[0] = dUlon + m_dU0;
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

	// Magnitude of the zonal wind jet
	double dUj;

	// Lapse rate troposphere
	double ddTdz;

	// Lapse rate stratosphere
	double ddTdzSTR;

	// Reference absolute temperature
	double dT0;

	// Parameter reference height for temperature disturbance
	double dhC;

	// Parameter reference length a for temperature disturbance
	double daC;

	// Parameter reference length for mountain profile
	double dlC;

	// Uniform diffusion coefficient scalars
	double dUCoeffS;

	// Uniform diffusion coefficient vectors
	double dUCoeffV;

	// No Rayleigh friction
	bool fNoRayleighFriction;

	// Parse the command line
	BeginTempestCommandLine("ShearJetMtnWave2DCartesianTest");
		SetDefaultResolutionX(288);
		SetDefaultResolutionY(48);
		SetDefaultLevels(32);
		SetDefaultOutputDeltaT("3h");
		SetDefaultDeltaT("300s");
		SetDefaultEndTime("12d");
		SetDefaultHorizontalOrder(4);
		SetDefaultVerticalOrder(1);

		CommandLineDouble(dbC, "b", 2.0);
		CommandLineDouble(dU0, "u0", 10.0);
		CommandLineDouble(dUj, "uj", 5.0);
		CommandLineDouble(ddTdz, "gamma", 0.0065);
		CommandLineDouble(ddTdzSTR, "gamma_str", -0.002);
		CommandLineDouble(dT0, "T0", 280.0);
		CommandLineDouble(dhC, "hC", 250.0);
		CommandLineDouble(daC, "aC", 5000.0);
		CommandLineDouble(dlC, "lC", 4000.0);
		CommandLineDouble(dUCoeffS, "nuDiffS", 0.0);
		CommandLineDouble(dUCoeffV, "nuDiffV", 0.0);
		CommandLineBool(fNoRayleighFriction, "norayleigh");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	// Setup the Model
	AnnounceBanner("MODEL SETUP");

	Model model(EquationSet::PrimitiveNonhydrostaticEquations);

	// Physical constants
	const PhysicalConstants & phys = model.GetPhysicalConstants();

	// Create a new instance of the test
	ShearJetMtnWave2DCartesianTest * test =
		new ShearJetMtnWave2DCartesianTest(phys, dbC,
				dU0,
				dUj,
				ddTdz,
				ddTdzSTR,
				dT0,
				dhC,
				daC,
				dlC,
				dUCoeffS,
				dUCoeffV,
				fNoRayleighFriction);

	// Setup the cartesian model with dimensions and reference latitude
	TempestSetupCartesianModel(model, test->m_dGDim, 0.0, 
								test->m_iLatBC, true);

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
	std::cout << "Try/catch block in the main program!" << std::endl;
}

	// Deinitialize Tempest
	TempestDeinitialize();
}

///////////////////////////////////////////////////////////////////////////////

