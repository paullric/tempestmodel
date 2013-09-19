///////////////////////////////////////////////////////////////////////////////
///
///	\file    SWTest2.cpp
///	\author  Paul Ullrich
///	\version February 24, 2013
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

#include "GridCSGLL.h"

#include "mpi.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Williamson et al. (1994) Test case 2
///
///		Geostrophically balanced flow with a shallow-water model.
///	</summary>
class ShallowWaterTestCase2 : public TestCase {

private:
	///	<summary>
	///		Background height field.
	///	</summary>
	double m_dH0;

	///	<summary>
	///		Maximum velocity.
	///	</summary>
	double m_dU0;

	///	<summary>
	///		Grid inclination.
	///	</summary>
	double m_dAlpha;

	///	<summary>
	///		True if a cosine bell tracer should be included.
	///	</summary>
	bool m_fTracerOn;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	ShallowWaterTestCase2(
		double dH0,
		double dU0,
		double dAlpha,
		bool fTracerOn
	) :
		m_dH0(dH0),
		m_dU0(dU0),
		m_dAlpha(dAlpha * M_PI / 180.0),
		m_fTracerOn(fTracerOn)
	{ }

public:
	///	<summary>
	///		Number of tracers used in this test.
	///	</summary>
	virtual int GetTracerCount() const {
		return (m_fTracerOn)?(1):(0);
	}

	///	<summary>
	///		Obtain test case specific physical constants.
	///	</summary>
	virtual void EvaluatePhysicalConstants(
		PhysicalConstants & phys
	) const {
		// Set the alpha parameter
		phys.SetAlpha(m_dAlpha);
	}

	///	<summary>
	///		Evaluate the topography at the given point.
	///	</summary>
	virtual double EvaluateTopography(
		const PhysicalConstants & phys,
		double dLon,
		double dLat
	) const {
		return 0.0;
	}

	///	<summary>
	///		Evaluate the state vector at the given point.
	///	</summary>
	virtual void EvaluatePointwiseState(
		const PhysicalConstants & phys,
		double dTime,
		double dZ,
		double dLon,
		double dLat,
		double * dState,
		double * dTracer
	) const {
		if (fabs(dLat - 0.5 * M_PI) < 1.0e-12) {
			dLat -= 1.0e-12;
		}
		if (fabs(dLat + 0.5 * M_PI) < 1.0e-12) {
			dLat += 1.0e-12;
		}

		// Pointwise zonal component of velocity
		dState[0] = m_dU0 * cos(dLat)
			* (cos(m_dAlpha) + cos(dLon) * tan(dLat) * sin(m_dAlpha));

		// Pointwise meridional component of velocity
		dState[1] = - m_dU0 * sin(dLon) * sin(m_dAlpha);

		// Height field
		double dHTrig =
			- cos(dLon) * cos(dLat) * sin(m_dAlpha)
			+ sin(dLat) * cos(m_dAlpha);

		dState[2] = m_dH0
			- (phys.GetEarthRadius() * phys.GetOmega() + 0.5 * m_dU0)
				* m_dU0 * dHTrig * dHTrig / phys.GetG();
/*
		dState[0] = 100.0;

		double dXX = cos(dLat) * cos(dLon) - 1.0;
		double dYY = cos(dLat) * sin(dLon);
		double dZZ = sin(dLat);

		double dR = (dXX * dXX + dYY * dYY + dZZ * dZZ);

		dState[0] += 100.0 * exp(-2.0 * dR);
*/

		// Cosine bell tracer field
		if (m_fTracerOn) {
			const double ParamLonC = 1.5 * M_PI;
			const double ParamLatC = 0.0;
			const double ParamH0 = 1000.0;
			const double ParamR = phys.GetEarthRadius() / 3.0;

			double dR = phys.GetEarthRadius()
				* acos(sin(ParamLatC) * sin(dLat)
					+ cos(ParamLatC) * cos(dLat) * cos(dLon - ParamLonC));

			if (dR < ParamR) {
				dTracer[0] = (ParamH0 / 2.0) * (1.0 + cos(M_PI * dR / ParamR));
			} else {
				dTracer[0] = 0.0;
			}
		}
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

	// Resolution
	int nResolution;

	// Order
	int nOrder;

	// Background velocity
	double dU0;

	// Background height field
	double dH0;

	// Grid rotation angle
	double dAlpha;

	// Include tracer field
	bool fTracersOn;

	// Output time
	double dOutputDeltaT;

	// Numerical method
	std::string strHorizontalDynamics;

	// Hyperviscosity
	bool fNoHyperviscosity;

	// Model parameters
	ModelParameters params;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strOutputDir, "output_dir", "out");
		CommandLineString(strOutputPrefix, "output_prefix", "out");
		CommandLineInt(nResolution, "resolution", 40);
		CommandLineInt(nOrder, "order", 4);
		CommandLineDouble(dU0, "u0", 38.61068277);
		CommandLineDouble(dH0, "h0", 2998.104995);
		CommandLineDouble(dAlpha, "alpha", 0.0);
		CommandLineBool(fTracersOn, "with_tracer");
		CommandLineDouble(params.m_dDeltaT, "dt", 200.0);
		CommandLineDouble(params.m_dEndTime, "endtime", 200.0);//86400.0 * 5.0);
		CommandLineDouble(dOutputDeltaT, "outputtime", 86400.0);
		CommandLineStringD(strHorizontalDynamics, "method", "SE", "(SE | DG)");
		CommandLineBool(fNoHyperviscosity, "nohypervis");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner("INITIALIZATION");

	// Create a new test case object
	AnnounceStartBlock("Creating test case");
	
	AnnounceEndBlock("Done");

	// Construct a model
	AnnounceStartBlock("Creating model");
	Model model(EquationSet::ShallowWaterEquations);
	AnnounceEndBlock("Done");

	// Generate a new cubed-sphere GLL grid
	AnnounceStartBlock("Creating grid");
	GridCSGLL grid(
		model,
		nResolution,
		4,
		nOrder,
		1,
		1);

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
	ShallowWaterTestCase2 test(dH0, dU0, dAlpha, fTracersOn);
	AnnounceStartBlock("Initializing test case");
	model.SetTestCase(&test);
	AnnounceEndBlock("Done");

	// Set the reference output manager for the model
	AnnounceStartBlock("Creating reference output manager");
	OutputManagerReference outmanRef(
		grid, dOutputDeltaT, strOutputDir, strOutputPrefix,
		720, 360, false, false, true, true);
	model.AttachOutputManager(&outmanRef);

	outmanRef.InitializeNcOutput("ref.nc");
	AnnounceEndBlock("Done");

	// Set the composite output manager for the model
	AnnounceStartBlock("Creating composite output manager");
	OutputManagerComposite outmanComp(
		grid, dOutputDeltaT, strOutputDir, strOutputPrefix);
	model.AttachOutputManager(&outmanComp);

	outmanComp.InitializeNcOutput("comp.nc");
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

