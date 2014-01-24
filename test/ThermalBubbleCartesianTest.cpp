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

#include "CommandLine.h"
#include "Announce.h"
#include "STLStringHelper.h"

#include "Model.h"
#include "TimestepSchemeStrang.h"
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
///		Thermal rising bubble test case.
///	</summary>
class ThermalBubbleCartesianTest : public TestCase {

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

	///	<summary>
	///		Flag indicating that the reference profile should be used.
	///	</summary>
	bool m_fNoReferenceState;

public:
	///	<summary>
	///		Constructor. (with physical constants defined privately here)
	///	</summary>
	ThermalBubbleCartesianTest(
		bool fNoReferenceState
	) :
		m_dH0(10000.),
		m_dThetaBar(300.0),
		m_dThetaC(0.5),
		m_drC(250.),
		m_dxC(500.),
		m_dzC(350.),
		m_dpiC(3.14159265),
		m_fNoReferenceState(fNoReferenceState)
	{
		// Set the dimensions of the box
		m_dGDim[0] = 0.0;
		m_dGDim[1] = 1000.0;
		m_dGDim[2] = -1000.0;
		m_dGDim[3] = 1000.0;
		m_dGDim[4] = 0.0;
		m_dGDim[5] = 1000.0;
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
		dState[2] = m_dThetaBar + EvaluateTPrime(phys, dXp, dZp);;

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

	// Off-centering
	double dOffCentering;

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
			"outThermalBubbleCartesianTest");
		CommandLineString(strOutputPrefix, "output_prefix", "out");
		CommandLineInt(nOutputsPerFile, "output_perfile", -1);
		CommandLineString(params.m_strRestartFile, "restart_file", "");
		CommandLineInt(nResolution, "resolution", 36);
		CommandLineInt(nLevels, "levels", 72);
		CommandLineInt(nHorizontalOrder, "order", 4);
		CommandLineInt(nVerticalOrder, "vertorder", 2);
		CommandLineDouble(dOffCentering, "offcentering", 0.0);
		CommandLineDouble(params.m_dDeltaT, "dt", 0.01);
		CommandLineDouble(params.m_dEndTime, "endtime", 700.0);
		CommandLineDouble(dOutputDeltaT, "outputtime", 50.0);
		CommandLineStringD(strHorizontalDynamics, "method", "SE", "(SE | DG)");
		CommandLineBool(fNoHyperviscosity, "nohypervis");
		CommandLineBool(fNoReferenceState, "norefstate");
		CommandLineInt(nVerticalHyperdiffOrder, "verticaldifforder", 0);
		CommandLineBool(fFullyExplicitVertical, "explicitvertical");
		CommandLineBool(fExnerPressureOnREdges, "exneredges");
		CommandLineBool(fMassFluxOnLevels, "massfluxlevels");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner("INITIALIZATION");

	// Construct a model
	AnnounceStartBlock("Creating model");
	Model model(EquationSet::PrimitiveNonhydrostaticEquations);
	AnnounceEndBlock("Done");

	// Set the parameters for the model
	AnnounceStartBlock("Initializing parameters");
	model.SetParameters(&params);
	AnnounceEndBlock("Done");

	// Set the timestep scheme
	TimestepSchemeStrang timestep(model, dOffCentering);
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
		true,
		!fExnerPressureOnREdges,
		fMassFluxOnLevels);

	AnnounceStartBlock("Initializing vertical dynamics");
	model.SetVerticalDynamics(&vdyn);
	AnnounceEndBlock("Done");

	// Generate a new cartesian GLL grid
	// Initialize the test case here (to have grid dimensions available)
	ThermalBubbleCartesianTest test(fNoReferenceState);

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

	grid.SetReferenceLength(1100000.0);

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
		1,
		false,  // Output variables in natural locations
		true);  // Remove reference profile in output

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
	ThermalBubbleCartesianTest ctest(fNoReferenceState);
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

