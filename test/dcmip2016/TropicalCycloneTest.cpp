///////////////////////////////////////////////////////////////////////////////
///
///	\file    TropicalCycloneTest.cpp
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

///////////////////////////////////////////////////////////////////////////////

 extern "C" {
	void tc_tropical_(
		double * dLon,
		double * dLat,
		double * dP,
        double * dZ,
        double * dU,
        double * dV,
        double * dT,
        double * dPhis,
        double * dPs,
        double * dRho,
        double * dQ
	);
} 

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		DCMIP 2016: Tropical Cyclone Test
///	</summary>
class TropicalCycloneTest : public TestCase {

protected:
	///	<summary>
	///		Model cap.
	///	</summary>
	double m_dZtop;

	///	<summary>
	///		Earth radius scaling parameter.
	///	</summary>
	double m_dEarthScaling;

	///	<summary>
	///		Earth rotation rate parameter.
	///	</summary>
	double m_dOmega;

	///	<summary>
	///		Background wind speed.
	///	</summary>
	double m_dU0;

	///	<summary>
	///		Background Brunt-Vaisala frequency.
	///	</summary>
	double m_dN;

	///	<summary>
	///		Surface temperature at the equator.
	///	</summary>
	double m_dTeq;

	///	<summary>
	///		Potential temperature perturbation width parameter (m).
	///	</summary>
	double m_dPertWidth;

	///	<summary>
	///		Longitudinal centerpoint of the potential temperature pert.
	///	</summary>
	double m_dPertLonC;

	///	<summary>
	///		Latitudinal centerpoint of the potential temperature pert.
	///	</summary>
	double m_dPertLatC;

	///	<summary>
	///		Magnitude of the potential temperature perturbation.
	///	</summary>
	double m_dPertMagnitude;

	///	<summary>
	///		Vertical wavelength of the potential temperature perturbation.
	///	</summary>
	double m_dPertVerticalWavelength;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	TropicalCycloneTest(
		double dZtop,
		double dEarthScaling,
		double dOmega,
		double dU0,
		double dN,
		double dTeq,
		double dPertWidth,
		double dPertLonC,
		double dPertLatC,
		double dPertMagnitude,
		double dPertVerticalWavelength
	) :
		m_dZtop(dZtop),
		m_dEarthScaling(dEarthScaling),
		m_dOmega(dOmega),
		m_dU0(dU0),
		m_dN(dN),
		m_dTeq(dTeq),
		m_dPertWidth(dPertWidth),
		m_dPertLonC(dPertLonC * M_PI / 180.0),
		m_dPertLatC(dPertLatC * M_PI / 180.0),
		m_dPertMagnitude(dPertMagnitude),
		m_dPertVerticalWavelength(dPertVerticalWavelength)
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
		double dZ,
		double dLon,
		double dLat,
		double * dState
	) const {

		// Store the state
		// State 0 = Zonal velocity (m/s)
		// State 1 = Meridional velocity (m/s)
		// State 2 = Theta (K)
		// State 3 = Vertical velocity (m/s)
		// State 4 = Density (kg/m^3)
		dState[0] = 0.0;
		dState[1] = 0.0;
		dState[2] = 0.5;
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
        
        
        double dRho;
        double dU;
        double dV;
        double dP;
        double dT;
        double dPhis;
        double dPs;
        double dQ;
        
        
		// Calculate the reference state
		//EvaluateReferenceState(phys, dZ, dLon, dLat, dState);
       tc_tropical_(&dLon, &dLat,&dP,&dZ,&dU,&dV,&dT,&dPhis,&dPs,&dRho,&dQ);
        
     
        dState[0] = dU;
		dState[1] = dV;
		dState[2] = dT*pow((phys.GetP0()/dP),(phys.GetR()/phys.GetCp()));
		dState[3] = 0.0;
		dState[4] = dRho;

     // std::cout << "z: " << dZ<< std::endl;
        
		// Add in the potential temperature perturbation
		double dR = phys.GetEarthRadius() * acos(
			sin(m_dPertLatC) * sin(dLat)
			+ cos(m_dPertLatC) * cos(dLat) * cos(dLon - m_dPertLonC));

		double dS = m_dPertWidth * m_dPertWidth /
			(m_dPertWidth * m_dPertWidth + dR * dR);

		dState[2] += m_dPertMagnitude
			* dS * sin(2.0 * M_PI * dZ / m_dPertVerticalWavelength);
	}
};

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	// Initialize Tempest
	TempestInitialize(&argc, &argv);
    

try {
	// Model cap.
	double dZtop;

	// Earth radius scaling parameter.
	double dEarthScaling;

	// Earth rotation rate parameter.
	double dOmega;

	// Background wind speed.
	double dU0;

	// Background Brunt-Vaisala frequency.
	double dN;

	// Surface temperature at the equator.
	double dTeq;

	// Potential temperature perturbation width parameter (m).
	double dPertWidth;

	// Longitudinal centerpoint of the potential temperature pert.
	double dPertLonC;

	// Latitudinal centerpoint of the potential temperature pert.
	double dPertLatC;

	// Magnitude of the potential temperature perturbation.
	double dPertMagnitude;

	// Vertical wavelength of the potential temperature perturbation.
	double dPertVerticalWavelength;

	// Parse the command line
	BeginTempestCommandLine("TropicalCycloneTest");
		SetDefaultResolution(20);
		SetDefaultLevels(10);
		SetDefaultOutputDeltaT("1500000u");
		SetDefaultDeltaT("1500000u");
		SetDefaultEndTime("1500000u");
		SetDefaultHorizontalOrder(4);
		SetDefaultVerticalOrder(1);

		CommandLineDouble(dZtop, "ztop", 10000.0);
		CommandLineDouble(dEarthScaling, "X", 125.0);
		CommandLineDouble(dOmega, "omega", 0.0);
		CommandLineDouble(dU0, "u0", 20.0);
		CommandLineDouble(dN, "N", 0.01);
		CommandLineDouble(dTeq, "Teq", 310.0);
		CommandLineDouble(dPertWidth, "d", 1.0/6.0);
		CommandLineDouble(dPertLonC, "lon_c", 3.14159265358979/9.0 );
		CommandLineDouble(dPertLatC, "lat_c", 2.0*3.14159265358979/9.0);
		CommandLineDouble(dPertMagnitude, "dtheta", 1.0);
		CommandLineDouble(dPertVerticalWavelength, "Lz", 20000.0);

		ParseCommandLine(argc, argv);
	EndTempestCommandLine(argv)

	// Setup the Model
	AnnounceBanner("MODEL SETUP");

	Model model(EquationSet::PrimitiveNonhydrostaticEquations);

	TempestSetupCubedSphereModel(model);
    

	// Set the test case for the model
	AnnounceStartBlock("Initializing test case");
	model.SetTestCase(
            new TropicalCycloneTest(
			dZtop,
			dEarthScaling,
			dOmega,
			dU0,
			dN,
			dTeq,
			dPertWidth,
			dPertLonC,
			dPertLatC,
			dPertMagnitude,
			dPertVerticalWavelength)
                      );
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

