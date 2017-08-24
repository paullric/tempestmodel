///////////////////////////////////////////////////////////////////////////////
///
///	\file    ShearJetMtnWave2DCartesianTestCBVF.cpp
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
class ShearJetMtnWave2DCartesianTestCBVF : public TestCase {

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
	///		Reference constant surface absolute temperature
	///	</summary>
	double m_dTheta0;

	///	<summary>
	///		Constant Brunt-Vaisala frequency
	///	</summary>
	double m_dNbar;

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
	ShearJetMtnWave2DCartesianTestCBVF(
		const PhysicalConstants & phys,
		double dbC,
		double dU0,
		double dUj,
		double dT0,
		double dNbar,
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
		m_dTheta0(dT0),
		m_dNbar(dNbar),
		m_dhC(dhC),
		m_daC(daC),
		m_dlC(dlC),
		m_dUCoeffS(dUCoeffS),
		m_dUCoeffV(dUCoeffV),
		m_fNoRayleighFriction(fNoRayleighFriction)
	{
		m_dpiC = M_PI;

		// Set the dimensions of the box
		m_dGDim[0] = -60000.0;
		m_dGDim[1] = 60000.0;
		m_dGDim[2] = -100.0;
		m_dGDim[3] = 100.0;
		m_dGDim[4] = 0.0;
		m_dGDim[5] = 35000.0;

		// Set the center of the domain in Y
		m_dY0 = 0.5 * (m_dGDim[3] - m_dGDim[2]);

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
		const double dRayleighDepth = 10000.0;
		const double dRayleighWidthR = 20000.0;
		const double dRayleighWidthL = 20000.0;
		const double dRayDepthXi = dRayleighDepth / m_dGDim[5];

		double dNuDepth = 0.0;
		double dNuRight = 0.0;
		double dNuLeft  = 0.0;

		double dLayerZ = 1.0 - dRayDepthXi;
		//double dLayerZ = m_dGDim[5] - dRayleighDepth;
 		double dLayerR = m_dGDim[1] - dRayleighWidthR;
 		double dLayerL = m_dGDim[0] + dRayleighWidthL;

		if (dZ > dLayerZ) {
			//double dNormZ = (m_dGDim[5] - dZ) / dRayleighDepth;
			double dNormZ = (1.0 - dZ) / dRayDepthXi;
			dNuDepth = 0.5 * dRayleighStrengthZ * (1.0 + cos(M_PI * dNormZ));
			//dNuDepth = dRayleighStrengthZ * pow(cos(0.5 * M_PI * dNormZ),4);
		}

		if (dXp > dLayerR) {
			double dNormX = (m_dGDim[1] - dXp) / dRayleighWidthR;
			dNuRight = 0.5 * dRayleighStrengthX * (1.0 + cos(M_PI * dNormX));
			//dNuRight = dRayleighStrengthX * pow(cos(0.5 * M_PI * dNormX),4);
		}
		if (dXp < dLayerL) {
			double dNormX = (dXp - m_dGDim[0]) / dRayleighWidthL;
			dNuLeft = 0.5 * dRayleighStrengthX * (1.0 + cos(M_PI * dNormX));
			//dNuLeft = dRayleighStrengthX * pow(cos(0.5 * M_PI * dNormX),4);
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
	///		Evaluate the reference state at the given point.
	///	</summary>
	virtual void EvaluateReferenceState(
		const PhysicalConstants & phys,
		double dXi,
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
		const double dLy = m_dGDim[3] - m_dGDim[2];

		double dExnerP = 0.0;
		double dRho = 0.0;
		double dThetaBar = m_dTheta0 * exp(m_dNbar * m_dNbar / dG * dZp);
		// Zero gravity case
		if (dG == 0.0) {
			static double dT0 = 300.0;

			dState[2] = dT0 * pow(dP0, (dRd / dCp));
			dState[4] = dP0 / (dRd * dT0);

		// Stratification with gravity
		} else {
			// Set the initial density based on the Exner pressure
			dExnerP = (dG * dG) / (dCp * m_dTheta0 * m_dNbar * m_dNbar);
			dExnerP *= (exp(-m_dNbar * m_dNbar / dG * dZp) - 1.0);
			dExnerP += 1.0;
			dRho = dP0 / (dRd * dThetaBar) * pow(dExnerP,(dCv / dRd));
			dState[4] = dRho;

			// Set the initial potential temperature field
			dState[2] = dThetaBar;
		}

		// Compute parameters for the jet (on terrain following grid)
		double dZcomp = m_dGDim[5] * dXi;
		dThetaBar = m_dTheta0 * exp(m_dNbar * m_dNbar / dG * dZcomp);
		// Zero gravity case
		if (dG == 0.0) {
			static double dT0 = 300.0;

			dThetaBar = dT0 * pow(dP0, (dRd / dCp));
			dRho = dP0 / (dRd * dT0);

		// Stratification with gravity
		} else {
			// Set the initial density based on the Exner pressure
			dExnerP = (dG * dG) / (dCp * m_dTheta0 * m_dNbar * m_dNbar);
			dExnerP *= (exp(-m_dNbar * m_dNbar / dG * dZcomp) - 1.0);
			dExnerP += 1.0;
			dRho = dP0 / (dRd * dThetaBar) * pow(dExnerP,(dCv / dRd));
		}

		// Calculate zonal velocity and set other velocity components
		double dEtaComp = phys.PressureFromRhoTheta(dThetaBar * dRho)
					/ dP0;
		double dExpDecay = exp(-(log(dEtaComp) / m_dbC) * (log(dEtaComp) / m_dbC));
		double dUlon = -m_dUj * log(dEtaComp) * dExpDecay;

		dState[0] = dUlon + m_dU0;
		dState[1] = 0.0;
		dState[3] = 0.0;
	}

	///	<summary>
	///		Evaluate the state vector at the given point.
	///	</summary>
	virtual void EvaluatePointwiseState(
		const PhysicalConstants & phys,
		const Time & time,
		double dXi,
		double dZp,
		double dXp,
		double dYp,
		double * dState,
		double * dTracer
	) const {

		// Evaluate the reference state at this point
		EvaluateReferenceState(phys, dXi, dZp, dXp, dYp, dState);

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

	// Reference absolute temperature
	double dT0;

	// Constant Brunt-Vaisala frequency
	double dNbar;

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
	BeginTempestCommandLine("ShearJetMtnWave2DCartesianTestCBVF");
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
		CommandLineDouble(dT0, "T0", 292.15);
		CommandLineDouble(dNbar, "N0", 0.01);
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
	ShearJetMtnWave2DCartesianTestCBVF * test =
		new ShearJetMtnWave2DCartesianTestCBVF(phys, dbC,
				dU0,
				dUj,
				dT0,
				dNbar,
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
