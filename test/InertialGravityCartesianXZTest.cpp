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
#include "VerticalDynamicsFEM.h"

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
		m_dP0(1.0E5),
		m_dCp(1005.0),
		m_dCv(718.0),
		m_dR(287.058),
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
		double dG = phys.GetG();
		double dThetaBar = m_dTheta0 * exp(m_dNbar * m_dNbar / dG * dZp);

		// Set the uniform U, V, W field for all time
		dState[0] = m_dU0;
		dState[1] = 0.0;
		dState[3] = 0.0;

		// Set the initial potential temperature field
		dState[2] = dThetaBar;

		// Set the initial density based on the Exner pressure
		double gsi = phys.GetG();
		double dExnerP = pow(gsi,2.0) / (m_dCp * m_dTheta0 *
                                        (m_dNbar * m_dNbar));
		dExnerP *= (exp(-pow(m_dNbar,2.0)/gsi * dZp) - 1.0);
		dExnerP += 1.0;
		double dRho = m_dP0 / (m_dR * dThetaBar) * pow(dExnerP,(m_dCv / m_dR));
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
		double dG = phys.GetG();
		double dThetaBar = m_dTheta0 * exp(m_dNbar * m_dNbar / dG * dZp);

		// Set the uniform U, V, W field for all time
		dState[0] = m_dU0;
		dState[1] = 0.0;
		dState[3] = 0.0;

		// Set the initial potential temperature field
		dState[2] = dThetaBar + EvaluateTPrime(phys, dXp, dZp);

		// Set the initial density based on the Exner pressure
		double gsi = phys.GetG();
		double dExnerP = (gsi * gsi) / (m_dCp * m_dTheta0 *
                                        (m_dNbar * m_dNbar));
		dExnerP *= (exp(-(m_dNbar * m_dNbar)/gsi * dZp) - 1.0);
		dExnerP += 1.0;
		double dRho = m_dP0 / (m_dR * dThetaBar) * pow(dExnerP,(m_dCv / m_dR));
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

	// Levels
	int nLevels;

	// Horizontal Order
	int nHorizontalOrder;

	// Vertical Order
	int nVerticalOrder;

	// Grid rotation angle
	double dAlpha;

	// Output time
	double dOutputDeltaT;

	// Numerical method
	std::string strHorizontalDynamics;

	// Use hyperdiffusion
	bool fNoHyperviscosity;

	// No reference state
	bool fNoReferenceState;

	// Vertical hyperdiffusion
	int nVerticalHyperdiffOrder;

	// Solve vertical using a fully explicit method
	bool fFullyExplicitVertical;

	// Store Exner pressure on edges (advanced)
	bool fExnerPressureOnREdges;

	// Store mass flux on levels (advanced)
	bool fMassFluxOnLevels;

	// Model parameters
	ModelParameters params;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strOutputDir, "output_dir",
			"outInertialGravityCartesianXZTest");
		CommandLineString(strOutputPrefix, "output_prefix", "out");
		CommandLineInt(nOutputsPerFile, "output_perfile", -1);
		CommandLineString(params.m_strRestartFile, "restart_file", "");
		CommandLineInt(nResolution, "resolution", 40);
		CommandLineInt(nLevels, "levels", 40);
		CommandLineInt(nHorizontalOrder, "order", 4);
		CommandLineInt(nVerticalOrder, "vertorder", 4);
		CommandLineDouble(dAlpha, "alpha", 0.0);
		CommandLineDouble(params.m_dDeltaT, "dt", 5.0);
		CommandLineDouble(params.m_dEndTime, "endtime", 3000.0);
		CommandLineDouble(dOutputDeltaT, "outputtime", 500.0);
		CommandLineStringD(strHorizontalDynamics, "method", "SE", "(SE | DG)");
		CommandLineBool(fNoHyperviscosity, "nohypervis");
		CommandLineBool(fNoReferenceState, "norefstate");
		CommandLineInt(nVerticalHyperdiffOrder, "verticaldifforder", 4);
		CommandLineBool(fFullyExplicitVertical, "explicitvertical");
		CommandLineBool(fExnerPressureOnREdges, "exneredges");
		CommandLineBool(fMassFluxOnLevels, "massfluxlevels");

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
	   model, nHorizontalOrder, eHorizontalDynamicsType, fNoHyperviscosity);
	AnnounceStartBlock("Initializing horizontal dynamics");
	model.SetHorizontalDynamics(&hdyn);
	AnnounceEndBlock("Done");

	// Set the vertical dynamics
	VerticalDynamicsFEM vdyn(
		model,
		nHorizontalOrder,
		nVerticalOrder,
		nVerticalHyperdiffOrder,
		fFullyExplicitVertical,
		!fExnerPressureOnREdges,
		fMassFluxOnLevels);

	AnnounceStartBlock("Initializing vertical dynamics");
	model.SetVerticalDynamics(&vdyn);
	AnnounceEndBlock("Done");

	// Generate a new cartesian GLL grid
	// Initialize the test case here (to have grid dimensions available)
	InertialGravityCartesianXZTest test(fNoReferenceState);

	AnnounceStartBlock("Constructing grid");
	GridCartesianGLL grid(
		model,
		nResolution,
		1,
		1,
		nHorizontalOrder,
		nVerticalOrder,
		nLevels,
		test.m_dGDim);

	model.SetGrid(&grid);
	AnnounceEndBlock("Done");

	// Set the reference output manager for the model
	AnnounceStartBlock("Creating reference output manager");
	OutputManagerReference outmanRef(
		grid,
		dOutputDeltaT,
		strOutputDir,
		strOutputPrefix,
		nOutputsPerFile,
		nResolution * (nHorizontalOrder - 1),
		2);
	model.AttachOutputManager(&outmanRef);
	AnnounceEndBlock("Done");

	// Set the composite output manager for the model
	AnnounceStartBlock("Creating composite output manager");
	OutputManagerComposite outmanComp(
		grid,
		dOutputDeltaT,
		strOutputDir,
		strOutputPrefix);
	model.AttachOutputManager(&outmanComp);
	AnnounceEndBlock("Done");

	// Set the checksum output manager for the model
	AnnounceStartBlock("Creating checksum output manager");
	OutputManagerChecksum outmanChecksum(grid, dOutputDeltaT);
	model.AttachOutputManager(&outmanChecksum);
	AnnounceEndBlock("Done");

	// Set the test case for the model
	InertialGravityCartesianXZTest ctest(fNoReferenceState);
	AnnounceStartBlock("Initializing test case");
	model.SetTestCase(&ctest);
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

