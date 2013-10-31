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

#include "CommandLine.h"
#include "Announce.h"
#include "STLStringHelper.h"

#include "Model.h"
#include "TimestepSchemeARK4.h"
#include "HorizontalDynamicsFEM.h"
#include "VerticalDynamicsStub.h"

#include "PhysicalConstants.h"
#include "TestCase.h"
#include "OutputManagerComposite.h"
#include "OutputManagerReference.h"
#include "OutputManagerChecksum.h"
#include "GridData4D.h"
#include "EquationSet.h"

#include "GridCartesianGLL.h"

#include "mpi.h"

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
	///		Reference surface pressure.
	///	</summary>
	double m_dP0;

    ///	<summary>
	///		Specific heat of air at constant pressure.
	///	</summary>
	double m_dCp;

    ///	<summary>
	///		Specific heat of air at constant volume.
	///	</summary>
	double m_dCv;

    ///	<summary>
	///		Gas constant of air.
	///	</summary>
	double m_dR;

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

public:
	///	<summary>
	///		Constructor. (with physical constants defined privately here)
	///	</summary>
	InertialGravityCartesianXZTest() :
		m_dH0(10000.),
		m_dU0(20.),
        m_dP0(1.E5),
        m_dCp(1005.0),
        m_dCv(718.0),
        m_dR(287.058),
		m_dNbar(0.01),
		m_dTheta0(300.),
		m_dThetaC(0.01),
		m_dhC(10000.),
		m_daC(5000.),
		m_dxC(100.E+3),
		m_dpiC(3.14159265)
	{
        m_dGDim[0] = 0.0; m_dGDim[1] = 300000.0;
        m_dGDim[2] = -1000.0; m_dGDim[3] = 1000.0;
        m_dGDim[4] = 0.0; m_dGDim[5] = 10000.0;
    }

public:
	///	<summary>
	///		Number of tracers used in this test.
	///	</summary>
	virtual int GetTracerCount() const {
		return 0;
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
	///		Evaluate the perturbed potential temperature field.
	///	</summary>
	double EvaluateTPrime(
        const PhysicalConstants phys,
		double dxP,
		double dzP
	) const {
        double gsi = phys.GetG();
        // Base potential temperature field
		double dThetaBar = m_dTheta0 * exp(pow(m_dNbar,2.0)/gsi * dzP);
        // Potential temperature perturbation
        double dThetaHat = m_dThetaC * sin(m_dpiC * dxP / m_dhC)
                                     / (1.0 + pow((dxP - m_dxC)/m_daC,2.0));
        //std::cout << dThetaHat + dThetaBar << '\n';

		return dThetaHat + dThetaBar;
	}

	///	<summary>
	///		Evaluate the state vector at the given point.
	///	</summary>
	virtual void EvaluatePointwiseState(
		const PhysicalConstants & phys,
		double dTime,
		double dzP,
        double dxP,
        double dyP,
		double *dState,
		double *dTracer
	) const {

        // Set the uniform U, V, W field for all time
        dState[0] = m_dU0;
        dState[1] = 0.0;
        dState[3] = 0.0;

        // Set the initial potential temperature field
        dState[2] = EvaluateTPrime(phys, dxP, dzP);

        // Set the initial density based on the Exner pressure
        double gsi = phys.GetG();
        double dExnerP = pow(gsi,2.0) / (m_dCp * m_dTheta0 * pow(m_dNbar,2.0));
        dExnerP *= exp(-pow(m_dNbar,2.0)/gsi * dzP) - 1.0;
        dExnerP += 1.0;
        //std::cout << dExnerP << '\n';
        double dRho = m_dP0 / (m_dR * dState[2]) * pow(dExnerP,(m_dCv / m_dR));
        //std::cout << dRho << '\n';
        dState[4] = dRho;
	}
};

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	// Initialize MPI
	MPI_Init(&argc, &argv);

try {
	// Output directory
	std::string strOutputDir;

	// Output file prefix
	std::string strOutputPrefix;

	// Number of outputs per reference file
	int nOutputsPerFile;

	// Resolution
	int nResolution;

	// Order
	int nOrder;

	// Grid rotation angle
	double dAlpha;

	// Output time
	double dOutputDeltaT;

	// Numerical method
	std::string strHorizontalDynamics;

	// Use hyperdiffusion
	bool fNoHyperviscosity;

	// Model parameters
	ModelParameters params;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strOutputDir, "output_dir",
			"outInertialGravityCartesianXZTest");
		CommandLineString(strOutputPrefix, "output_prefix", "out");
		CommandLineInt(nOutputsPerFile, "output_perfile", -1);
		CommandLineInt(nResolution, "resolution", 20);
		CommandLineInt(nOrder, "order", 4);
		CommandLineDouble(dAlpha, "alpha", 0.0);
		CommandLineDouble(params.m_dDeltaT, "dt", 1.0);
		CommandLineDouble(params.m_dEndTime, "endtime", 2.0);
		CommandLineDouble(dOutputDeltaT, "outputtime", 1.0);
		CommandLineStringD(strHorizontalDynamics, "method", "DG", "(SE | DG)");
		CommandLineBool(fNoHyperviscosity, "nohypervis");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner("INITIALIZATION");

	// Create a new test case object
	AnnounceStartBlock("Creating test case");

	AnnounceEndBlock("Done");

	// Construct a model
	AnnounceStartBlock("Creating model");
	Model model(EquationSet::PrimitiveNonhydrostaticEquations);
	AnnounceEndBlock("Done");

	// Generate a new cartesian GLL grid (20 x 20 x 20 for now)
    // Initialize the test case here (to have grid dimensions available)
    InertialGravityCartesianXZTest test;
	AnnounceStartBlock("Creating grid");
	GridCartesianGLL grid(
		model,
		nResolution,
		1,
		nOrder,
        nOrder,
		nResolution,
        test.m_dGDim);

	AnnounceEndBlock("Done");

	// Set the parameters for the model
	AnnounceStartBlock("Initializing parameters");
	model.SetParameters(&params);
	AnnounceEndBlock("Done");

	// Set the timestep scheme
	TimestepSchemeARK4 timestep(model);
	AnnounceStartBlock("Initializing timestep scheme");
	model.SetTimestepScheme(&timestep);
	AnnounceEndBlock("Done");

	// Set the horizontal dynamics
	HorizontalDynamicsFEM::Type eHorizontalDynamicsType;
	STLStringHelper::ToLower(strHorizontalDynamics);
	if (strHorizontalDynamics == "se") {
		eHorizontalDynamicsType = HorizontalDynamicsFEM::SpectralElement;
	} else if (strHorizontalDynamics == "dg") {
		eHorizontalDynamicsType = HorizontalDynamicsFEM::DiscontinuousGalerkin;
	} else {
		_EXCEPTIONT("Invalid method: Expected \"SE\" or \"DG\"");
	}

	HorizontalDynamicsFEM hdyn(
		model, nOrder, eHorizontalDynamicsType, fNoHyperviscosity);
	AnnounceStartBlock("Initializing horizontal dynamics");
	model.SetHorizontalDynamics(&hdyn);
	AnnounceEndBlock("Done");

	// Set the vertical dynamics
	VerticalDynamicsStub vdyn(model);
	AnnounceStartBlock("Initializing vertical dynamics");
	model.SetVerticalDynamics(&vdyn);
	AnnounceEndBlock("Done");

	// Set the grid for the model
	AnnounceStartBlock("Initializing grid");
	model.SetGrid(&grid);
	AnnounceEndBlock("Done");

	// Set the test case for the model
    InertialGravityCartesianXZTest ctest;
	AnnounceStartBlock("Initializing test case");
	model.SetTestCase(&ctest);
	AnnounceEndBlock("Done");

	// Set the reference output manager for the model
	AnnounceStartBlock("Creating reference output manager");
	OutputManagerReference outmanRef(
		grid,
		dOutputDeltaT,
		strOutputDir,
		strOutputPrefix,
		nOutputsPerFile,
		20, 20);
	outmanRef.OutputVorticity();
	outmanRef.OutputDivergence();
	model.AttachOutputManager(&outmanRef);
	AnnounceEndBlock("Done");

	// Set the composite output manager for the model
	AnnounceStartBlock("Creating composite output manager");
	OutputManagerComposite outmanComp(
		grid, dOutputDeltaT, strOutputDir, strOutputPrefix);
	model.AttachOutputManager(&outmanComp);
	AnnounceEndBlock("Done");

	// Set the checksum output manager for the model
	AnnounceStartBlock("Creating checksum output manager");
	OutputManagerChecksum outmanChecksum(grid, dOutputDeltaT);
	model.AttachOutputManager(&outmanChecksum);
	AnnounceEndBlock("Done");

	// Begin execution
	AnnounceBanner("SIMULATION");
	model.Go();

	// Execution complete
	//AnnounceBanner("Execution completed successfully");

	// Compute error norms
	AnnounceBanner("RESULTS");
	model.ComputeErrorNorms();
	AnnounceBanner();

} catch(Exception & e) {
	std::cout << e.ToString() << std::endl;
}

	// Finalize MPI
	MPI_Finalize();
}

///////////////////////////////////////////////////////////////////////////////

