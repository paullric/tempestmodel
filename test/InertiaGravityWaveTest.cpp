///////////////////////////////////////////////////////////////////////////////
///
///	\file    InertiaGravityWaveTest.cpp
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

#include "Defines.h"
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

#include "GridCSGLL.h"

#ifdef USE_PETSC
#include <petscsnes.h>
#else
#include "mpi.h"
#endif

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Ullrich, Melvin, Jablonowski and Staniforth (2013) Baroclinic wave test
///	</summary>
class InertiaGravityWaveTest : public TestCase {

protected:
	///	<summary>
	///		Background temperature.
	///	</summary>
	double ParamT0;

	///	<summary>
	///		Model cap.
	///	</summary>
	double ParamZtop;

	///	<summary>
	///		Potential temperature perturbation.
	///	</summary>
	double ParamPertMagnitude;

protected:
	///	<summary>
	///		Earth radius scaling parameter.
	///	</summary>
	double m_dEarthRadiusScaling;

	///	<summary>
	///		Model height cap.
	///	</summary>
	double m_dZtop;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	InertiaGravityWaveTest(
		double dEarthRadiusScaling
	) :
		ParamT0(300.0),
		ParamZtop(10000.0),
		ParamPertMagnitude(1.0),

		m_dEarthRadiusScaling(dEarthRadiusScaling)
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
		return ParamZtop;
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
		phys.SetOmega(0.0);
		phys.SetEarthRadius(phys.GetEarthRadius() / m_dEarthRadiusScaling);
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

		// Calculate the isothermal pressure
		double dPressure =
			phys.GetP0() * exp(- phys.GetG() * dZ / phys.GetR() / ParamT0);

		// Calculate exact density
		double dRho = dPressure / (phys.GetR() * ParamT0);

		// Store the state
		dState[0] = 0.0;
		dState[1] = 0.0;
		dState[2] = phys.RhoThetaFromPressure(dPressure) / dRho;
		dState[3] = 0.0;
		dState[4] = dRho;
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

		// Calculate the isothermal pressure
		double dPressure =
			phys.GetP0() * exp(- phys.GetG() * dZ / (phys.GetR() * ParamT0));

		// Calculate exact density
		double dRho = dPressure / (phys.GetR() * ParamT0);

		// Calculate exact temperature
		double dTemperature = dPressure / (dRho * phys.GetR());

		dTemperature += ParamPertMagnitude
			* exp(-100.0 * dLat * dLat)
			* sin(M_PI * dZ / ParamZtop);

		// Recalculate Rho
		dRho = dPressure / (phys.GetR() * dTemperature);

/*
		// Add the perturbation
		double dTb = ParamPertMagnitude
			* exp(100.0 * (sin(dLat) - 1.0))
			* sin(M_PI * dZ / ParamZtop);

		double dRhob =
			phys.GetP0() / (phys.GetR() * ParamT0) * (-dTb / ParamT0);

		dRho += exp(- 0.5 * phys.GetG() * dZ / phys.GetR() / ParamT0) * dRhob;
*/
		// Calculate the potential temperature
		double dTheta = phys.RhoThetaFromPressure(dPressure) / dRho;
/*
		// Add the perturbation
		dTheta += ParamPertMagnitude
			* exp(-100.0 * dLat * dLat)
			* sin(M_PI * dZ / ParamZtop);
*/
		// Store the state
		dState[0] = 0.0;
		dState[1] = 0.0;
		dState[2] = dTheta;
		dState[3] = 0.0;
		dState[4] = dRho;
	}
};

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

#ifdef USE_PETSC
	// Initialize PetSc
	PetscInitialize(&argc,&argv, NULL, NULL);
#else
	// Initialize MPI
	MPI_Init(&argc, &argv);
#endif

try {
	// Output directory
	std::string strOutputDir;

	// Output file prefix
	std::string strOutputPrefix;

	// Number of outputs per reference file
	int nOutputsPerFile;

	// Resolution
	int nResolution;

	// Vertical resolution
	int nLevels;

	// Order
	int nHorizontalOrder;

	// Vertical Order
	int nVerticalOrder;

	// Use reference state flag
	bool fNoReferenceState;

	// Solve vertical using a fully explicit method
	bool fFullyExplicitVertical;

	// Earth radius scaling
	double dEarthRadiusScaling;

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
			"outInertiaGravityWaveTest");
		CommandLineString(strOutputPrefix, "output_prefix", "out");
		CommandLineInt(nOutputsPerFile, "output_perfile", -1);
		CommandLineInt(nResolution, "resolution", 20);
		CommandLineInt(nLevels, "levels", 10);
		CommandLineInt(nHorizontalOrder, "order", 4);
		CommandLineInt(nVerticalOrder, "vertorder", 1);
		CommandLineBool(fNoReferenceState, "norefstate");
		CommandLineDouble(dEarthRadiusScaling, "radiusscale", 1.0);
		CommandLineDouble(params.m_dDeltaT, "dt", 200.0);
		CommandLineDouble(params.m_dEndTime, "endtime", 200.0);
		CommandLineDouble(dOutputDeltaT, "outputtime", 21600.0);
		CommandLineStringD(strHorizontalDynamics, "method", "SE", "(SE | DG)");
		CommandLineBool(fNoHyperviscosity, "nohypervis");
		CommandLineBool(fFullyExplicitVertical, "explicitvertical");

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
	TimestepSchemeStrang timestep(model);
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

	AnnounceStartBlock("Initializing horizontal dynamics");
	HorizontalDynamicsFEM hdyn(
		model,
		nHorizontalOrder,
		eHorizontalDynamicsType,
		fNoHyperviscosity);

	model.SetHorizontalDynamics(&hdyn);
	AnnounceEndBlock("Done");

	// Set the vertical dynamics
	VerticalDynamicsFEM vdyn(
		model,
		nHorizontalOrder,
		nVerticalOrder,
		0,
		fFullyExplicitVertical,
		!fNoReferenceState);

	AnnounceStartBlock("Initializing vertical dynamics");
	model.SetVerticalDynamics(&vdyn);
	AnnounceEndBlock("Done");

	// Construct the cubed-sphere grid for the model
	AnnounceStartBlock("Constructing grid");
	GridCSGLL grid(
		model,
		nResolution,
		4,
		nHorizontalOrder,
		nVerticalOrder,
		nLevels);

	grid.SetReferenceLength(0.5 * M_PI / 30.0 * dEarthRadiusScaling);

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
		360, 180,
		false,   // Output variables in natural locations
		true);   // Remove reference profile in output

	outmanRef.OutputVorticity();
	outmanRef.OutputDivergence();
	outmanRef.OutputTemperature();
	model.AttachOutputManager(&outmanRef);
	AnnounceEndBlock("Done");

	// Set the checksum output manager for the model
	AnnounceStartBlock("Creating checksum output manager");
	OutputManagerChecksum outmanChecksum(grid, dOutputDeltaT);
	model.AttachOutputManager(&outmanChecksum);
	AnnounceEndBlock("Done");

	// Set the test case for the model
	InertiaGravityWaveTest test(dEarthRadiusScaling);

	AnnounceStartBlock("Initializing test case");
	model.SetTestCase(&test);
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

#ifdef USE_PETSC
	// Finalize PetSc
	PetscFinalize();
#else
	// Finalize MPI
	MPI_Finalize();
#endif
}

///////////////////////////////////////////////////////////////////////////////

