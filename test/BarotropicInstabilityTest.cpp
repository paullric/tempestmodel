///////////////////////////////////////////////////////////////////////////////
///
///	\file    BarotropicInstability.cpp
///	\author  Paul Ullrich
///	\version August 13, 2013
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
///		Galewsky et al. (2004)
///
///		Barotropic instability test case.
///	</summary>
class BarotropicInstabilityTestCase : public TestCase {

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
	///		Southern latitude of barotropic jet.
	///	</summary>
	double m_dTheta0;

	///	<summary>
	///		Northern latitude of barotropic jet.
	///	</summary>
	double m_dTheta1;

	///	<summary>
	///		Stretching parameter.
	///	</summary>
	double m_dXE;

	///	<summary>
	///		Parameter used for computing height field (??)
	///	</summary>
	double m_dHHat;

	///	<summary>
	///		Parameter used for computing height field (??)
	///	</summary>
	double m_dHPhi2;

	///	<summary>
	///		Parameter used for computing height field (??)
	///	</summary>
	double m_dHAlpha;

	///	<summary>
	///		Parameter used for computing height field (??)
	///	</summary>
	double m_dHBeta;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	BarotropicInstabilityTestCase(
		double dAlpha
	) :
		m_dH0(10158.18617045463179),
		m_dU0(80.0),
		m_dAlpha(dAlpha),
		m_dTheta0(M_PI / 7.0),
		m_dTheta1(M_PI / 2.0 - m_dTheta0),
		m_dXE(0.3),
		m_dHHat(120.0),
		m_dHPhi2(M_PI / 4.0),
		m_dHAlpha(1.0 / 3.0),
		m_dHBeta(1.0 / 15.0)
	{ }

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
	///		Determined the corresponding point in rotated spherical coordinates.
	///	</summary>
	void CalculateRLLPrime(
		double dLon,
		double dLat,
		double & dLonP,
		double & dLatP
	) const {
		if (m_dAlpha == 0.0) {
			dLonP = dLon;
			dLatP = dLat;
			return;
		}

		// Calculate latitude and longitude in rotated coordinates
		dLatP = asin(sin(dLat) * cos(m_dAlpha)
			- cos(dLat) * cos(dLon) * sin(m_dAlpha));

		dLonP = asin(sin(dLon) * cos(dLat) / cos(dLatP));

		double dTemp =
			cos(m_dAlpha) * cos(dLon) * cos(dLat)
			+ sin(m_dAlpha) * sin(dLat);

		if (dTemp < 0.0) {
			dLonP = M_PI - dLonP;
		}
		if (dLonP < 0.0) {
			dLonP += 2.0 * M_PI;
		}
	}

	///	<summary>
	///		Evaluate the perturbed velocity field.
	///	</summary>
	double EvaluateUPrime(
		double dLonP,
		double dLatP
	) const {
		if (dLatP < m_dTheta0) {
			return 0.0;

		} else if (dLatP > m_dTheta1) {
			return 0.0;
		}

		double dNormalizer =
			exp(- 4.0 / (m_dTheta1 - m_dTheta0) / (m_dTheta1 - m_dTheta0));

		double dUp =
			exp(1.0 / (dLatP - m_dTheta0) / (dLatP - m_dTheta1));

		return (m_dU0 / dNormalizer * dUp);
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

		// Calculate rotated RLL coordinates of this point
		double dLonP;
		double dLatP;

		CalculateRLLPrime(dLon, dLat, dLonP, dLatP);

		// Calculate height field via numerical integration
		int nIntervals =
			static_cast<unsigned int>((dLatP + 0.5 * M_PI) / (1.0e-2));

		if (nIntervals < 1) {
			nIntervals = 1;
		}

		// Numerically integrate H using 4th order Gaussian quadrature
		DataVector<double> dLatX;
		dLatX.Initialize(nIntervals+1);

		for (int i = 0; i <= nIntervals; i++) {
			dLatX[i] = - 0.5 * M_PI +
				((dLatP + 0.5 * M_PI) / static_cast<double>(nIntervals))
					* static_cast<double>(i);
		}

		double dH = 0.0;

		double dXeval;
		double dU;

		for (int i = 0; i < nIntervals; i++) {
		for (int m = -1; m <= 1; m += 2) {
			dXeval =
				0.5 * (dLatX[i+1] + dLatX[i])
				+ static_cast<double>(m)
					* sqrt(1.0 / 3.0) * 0.5 * (dLatX[i+1] - dLatX[i]);

			dU = EvaluateUPrime(dLonP, dXeval);

			dH += (2.0 * phys.GetEarthRadius() * phys.GetOmega() * sin(dXeval)
				+ dU * tan(dXeval)) * dU;
		}
		}

		dH *= 0.5 * (dLatX[1] - dLatX[0]);

		dState[2] = m_dH0 - dH / phys.GetG();

		// Add perturbation
		if (dLon > M_PI) {
			dLon = dLon - 2.0 * M_PI;
		}
		if ((dLon < -M_PI) || (dLon > M_PI)) {
			_EXCEPTIONT("Invalid value of longitude.");
		}

		dState[2] += m_dHHat * cos(dLat)
			* exp(-((dLon * dLon) / (m_dHAlpha * m_dHAlpha)))
			* exp(-((m_dHPhi2 - dLat) * (m_dHPhi2 - dLat)
				/ (m_dHBeta * m_dHBeta)));

		// Evaluate the velocity field
		double dUP = EvaluateUPrime(dLonP, dLatP);

		double dUlat =
			(- dUP * sin(m_dAlpha) * sin(dLonP)) / cos(dLat);

		double dUlon;
		if (fabs(cos(dLon)) < 1.0e-13) {
			if (fabs(m_dAlpha) > 1.0e-13) {
				if (cos(dLon) > 0.0) {
					dUlon = - dUlat * cos(dLat) / tan(m_dAlpha);
				} else {
					dUlon = dUlat * cos(dLat) / tan(m_dAlpha);
				}
			} else {
				dUlon = dUP;
			}

		} else {
			dUlon = (dUlat * sin(dLat) * sin(dLon) + dUP * cos(dLonP))
				/ cos(dLon);
		}

		if (dUlon > 80.0) {
			printf("%1.5e %1.5e\n", dUlon, dUlat);
		}

		dState[0] = dUlon;
		dState[1] = dUlat;
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
			"outBarotropicInstabilityTest");
		CommandLineString(strOutputPrefix, "output_prefix", "out");
		CommandLineInt(nOutputsPerFile, "output_perfile", -1);
		CommandLineInt(nResolution, "resolution", 20);
		CommandLineInt(nOrder, "order", 4);
		CommandLineDouble(dAlpha, "alpha", 0.0);
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

	// Construct the cubed-sphere grid for the model
	AnnounceStartBlock("Constructing grid");
	GridCSGLL grid(
		model,
		nResolution,
		4,
		nOrder,
		1,
		1);

	grid.AddDefaultPatches();
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
		720, 360);
	outmanRef.OutputVorticity();
	outmanRef.OutputDivergence();
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
	BarotropicInstabilityTestCase test(dAlpha);
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

	// Finalize MPI
	MPI_Finalize();
}

///////////////////////////////////////////////////////////////////////////////

