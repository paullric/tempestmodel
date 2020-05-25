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

	///     <summary>
        ///             Height of the Rayleigh damped layer.
        ///     </summary>
        double m_dZh;

	///     <summary>
        ///             Rayleigh friction time scale.
        ///     </summary>
        double m_dTau0;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	HeldSuarezTest(
		double dZh,
		double dTau0,
		double dZtop
	) :
		ParamEarthRadiusScaling(1.0),
		ParamT0(280.0),

		m_dZh(dZh),
		m_dTau0(dTau0),
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
        
        	// Random perturbation on U and V to induce asymmetry
        	dState[0] += 1.0E-3 * static_cast<double>(rand()) / 
                        static_cast<double>(RAND_MAX);
        	dState[1] += 1.0E-3 * static_cast<double>(rand()) / 
                        static_cast<double>(RAND_MAX);
	}
	
    	///	<summary>
	///		Evaluate the perturbed state vector at the given point upon restart.
	///	</summary>
	virtual void EvaluatePointwisePerturbation(
		const PhysicalConstants & phys,
		const Time & time,
		double dZ,
		double dLon,
		double dLat,
		double * dState,
		double * dTracer
	) const {
		const int fMode = 2;

		double dH = phys.GetR() * ParamT0 / phys.GetG();

		double dP = phys.GetP0() * exp(- dZ / dH);

		double dRho = dP / phys.GetG() / dH;
        
        double dAE = phys.GetEarthRadius();

		// Store the state perturbation heating magnitude and scale parameters
        	double dPert = 1.0;
        	double dXLS = 5.0E6;
        	double dYLS = 1.2E6;
        	double dZLS = 5.0E3;
        	double dLonShift = 0.0;
        
        	// Apply a periodic shift so as to specify the perturbation function
        	if (dLon > M_PI) {
            		dLonShift = dLon - 2.0 * M_PI;
        	} else {
            		dLonShift = dLon;
       	 	}

		double dXL = dLonShift * dAE * std::cos(dLat) / dXLS;
		double dYL = dLat * dAE / dYLS;
        
        // Create zero mean double Gaussian perturbation in U-Theta
        double dZHeat = 0.2;
        double dPow = 1.0 / dZHeat - 1.0;
        double dAp = 1.0 / dZHeat * std::pow(1.0 - dZHeat, -dPow);
        double dXi = dZ / m_dZtop;

		double dkC = 0.0;
		double dFX_m1 = 0.0;
		double dFX_m2 = 0.0;
		double dDFXDX_m1 = 0.0;
		double dDFXDX_m2 = 0.0;
		
		double dGX = std::exp(-0.5 * dXL * dXL);
		double dGY = std::exp(-0.5 * dYL * dYL);
		
		double dVxi = dAp * std::pow(1.0 - dXi, dPow) * dXi;
		
		double dIntVxi = dAp / 30.0 * (1.0 - std::pow((1.0 - dXi), 5.0) * 
			(1.0 + 5.0 * dXi));

		double dInt2Vxi = dAp / 30.0 * (2.0 / 7.0 * dXi + std::pow(1.0 - dXi, 6.0) - 
				5.0 / 7.0 * std::pow(1.0 - dXi, 7.0) - 2.0 / 7.0);

		// Assuming c = 50 m/s for shallow waver gravity wave speed
		double dbetap = 2.0 * phys.GetOmega() * std::cos(dLat);
        
	        double dRTscale = dPert / ParamT0 *
                        (1.0 - phys.GetKappa()) / phys.GetR();
		
		double dUscale = dAE * (m_dZtop / dYLS) * phys.GetG() * 
                        (dPert / ParamT0) * (1.0 / dbetap);
		
		double dWscale = (m_dZtop / dXLS) * dUscale;

		//if (dLat <= 1.0E-6) {
		//	printf("%.16E %.16E \n",dUscale, dWscale);
		//}
		
		if (fMode == 1) {
			// SPECIFIES MODE 1 - with average removed where dkC is a 
            // constant that approximates the Gaussian
            dkC = 0.5 / (dXLS * dXLS) * std::sqrt(2.25);
            dFX_m1 = 1.0 - std::tanh(dkC * dLonShift) *
                    std::tanh(dkC * dLonShift) -
                    1.0 / (M_PI * dkC) * std::tanh(M_PI * dkC);

			dState[0] = dUscale * dIntVxi * dFX_m1 * dGY;
            dState[2] = dPert * dVxi * dFX_m1 * dGY;
			
			dState[3] = 0.0;
		} else if (fMode == 2) {
			// SPECIFIES MODE 2
            dFX_m2 = - std::pow(std::exp(1.0), 0.5) * dXL *
        	                std::exp(-0.5 * dXL * dXL);

			dState[0] = dUscale * dIntVxi * dFX_m2 * dGY;
			dState[2] = dRTscale * dP * dVxi * dFX_m2 * dGY;
			//if (dLat <= 1.0E-6) {
			//	printf("%.16E %.16E \n",dState[0], dState[2]);
			//}

			dDFXDX_m2 = dGX * (1.0 - (dXL * dXL));
			dState[3] = 0.0; //dWscale * dInt2Vxi * dDFXDX_m2 * dGY;
		} else {
			AnnounceBanner("PERTURBATION ASSUMED ZONALLY CONSTANT...");
			dState[0] = 0.0;
                        dState[2] = 0.0;
			dState[3] = 0.0;
		}
      
		dState[1] = 0.0;
		dState[4] = 0.0;
        	//std::cout << dState[2] << std::endl;
	}

	///     <summary>
        ///             Flag indicating whether or not Rayleigh friction strength is given.
        ///     </summary>
        virtual bool HasRayleighFriction() const {
                return true;
        }

        ///     <summary>
        ///             Evaluate the Rayleigh friction strength at the given point.
        ///     </summary>
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

                return (0.0 * dNuDepth / m_dTau0);
        }
};

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	// Initialize Tempest
	TempestInitialize(&argc, &argv);

try {
	// Model height cap
	double dZtop;
	
	// Height of the Rayleigh damped layer.
        double dZh;

        // Rayleigh friction time scale.
        double dTau0;
	

	// Parse the command line
	BeginTempestCommandLine("HeldSuarezTest");
		SetDefaultResolution(20);
		SetDefaultLevels(10);
		SetDefaultOutputDeltaT("200s");
		SetDefaultDeltaT("200s");
		SetDefaultEndTime("200s");
		SetDefaultHorizontalOrder(4);
		SetDefaultVerticalOrder(1);

		CommandLineDouble(dZtop, "ztop", 35000.0);
		CommandLineDouble(dZh, "zh", 30000.0);
                CommandLineDouble(dTau0, "tau0", 25.0);

		ParseCommandLine(argc, argv);
	EndTempestCommandLine(argv)

	// Setup the Model
	AnnounceBanner("MODEL SETUP");

	Model model(EquationSet::PrimitiveNonhydrostaticEquations);

	TempestSetupCubedSphereModel(model);

	// Set the test case for the model
	AnnounceStartBlock("Initializing test case");

	model.SetTestCase(new HeldSuarezTest(dZh, dTau0, dZtop));

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
