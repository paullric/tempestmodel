///////////////////////////////////////////////////////////////////////////////
///
///	\file    TempestInitialize.h
///	\author  Paul Ullrich
///	\version June 4, 2014
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

#ifndef _TEMPESTINITIALIZE_H_
#define _TEMPESTINITIALIZE_H_

#include "Defines.h"
#include "Announce.h"
#include "Model.h"
#include "TimestepSchemeStrang.h"
#include "TimestepSchemeERK.h"
#include "TimestepSchemeARS222.h"
#include "TimestepSchemeARS232.h"
#include "TimestepSchemeARK232.h"
#include "TimestepSchemeGARK2.h"
#include "TimestepSchemeARS343.h"
#include "TimestepSchemeARS343b.h"
#include "TimestepSchemeARS443.h"
#include "TimestepSchemeSSP3332.h"
#include "TimestepSchemeSplitExp.h"
#include "HorizontalDynamicsStub.h"
#include "HorizontalDynamicsFEM.h"
#include "HorizontalDynamicsFEMV2.h"
#include "SplitExplicitDynamics.h"
#include "HighSpeedDynamics.h"
#include "VerticalDynamicsStub.h"
#include "VerticalDynamicsFEM.h"
#include "VerticalDynamicsFEMV2.h"
#include "VerticalDynamicsSchur.h"
#include "OutputManagerComposite.h"
#include "OutputManagerReference.h"
#include "OutputManagerChecksum.h"
#include "GridCSGLL.h"
#include "GridCartesianGLL.h"
#include "VerticalStretch.h"

#include "TimeObj.h"
#include "Announce.h"
#include "CommandLine.h"
#include "STLStringHelper.h"

#ifdef TEMPEST_MPIOMP
#include <mpi.h>
#endif

#include <string>

#ifdef TEMPEST_PETSC
#include <petscsnes.h>
#endif

///////////////////////////////////////////////////////////////////////////////

struct _TempestCommandLineVariables {
	bool fNoOutput;
	std::string strOutputDir;
	std::string strOutputPrefix;
	std::string strRestartFile;
	int nOutputsPerFile;
	Time timeOutputDeltaT;
	Time timeOutputRestartDeltaT;
	Time timeDeltaT;
	Time timeEndTime;
	int nOutputResX;
	int nOutputResY;
	int nOutputResZ;
	bool fOutputVorticity;
	bool fOutputDivergence;
	bool fOutputTemperature;
	bool fOutputSurfacePressure;
	bool fOutputRichardson;
	bool fOutputConvective;
	bool fOutputZonalForce;
	bool fNoReferenceState;
	bool fNoTracers;
	bool fNoHyperviscosity;
	int nHyperviscosityOrder;
	double dNuScalar;
	double dNuDiv;
	double dNuVort;
	double dInstepNuDiv;
	bool fExplicitVertical;
	std::string strVerticalStaggering;
	bool fForceMassFluxOnLevels;
	std::string strVerticalStretch;
	std::string strVerticalDiscretization;
	int nVerticalHyperdiffOrder;
	std::string strTimestepScheme;
	std::string strHorizontalDynamics;
	std::string strVerticalDynamics;
	int nResolutionX;
	int nResolutionY;
	int nLevels;
	int nHorizontalOrder;
	int nVerticalOrder;
};

///////////////////////////////////////////////////////////////////////////////

#define _TempestDefineCommandLineDefault(TestCaseName) \
	CommandLineBool(_tempestvars.fNoOutput, "output_none"); \
	CommandLineString(_tempestvars.strOutputDir, "output_dir", "out" TestCaseName); \
	CommandLineString(_tempestvars.strOutputPrefix, "output_prefix", "out"); \
	CommandLineString(_tempestvars.strRestartFile, "restart_file", ""); \
	CommandLineInt(_tempestvars.nOutputsPerFile, "output_perfile", -1); \
	CommandLineDeltaTime(_tempestvars.timeOutputRestartDeltaT, "output_restart_dt", ""); \
	CommandLineInt(_tempestvars.nOutputResX, "output_x", 360); \
	CommandLineInt(_tempestvars.nOutputResY, "output_y", 180); \
	CommandLineInt(_tempestvars.nOutputResZ, "output_z", 0); \
	CommandLineBool(_tempestvars.fOutputVorticity, "output_vort"); \
	CommandLineBool(_tempestvars.fOutputDivergence, "output_div"); \
	CommandLineBool(_tempestvars.fOutputTemperature, "output_temp"); \
	CommandLineBool(_tempestvars.fOutputSurfacePressure, "output_ps"); \
	CommandLineBool(_tempestvars.fOutputRichardson, "output_Ri"); \
	CommandLineBool(_tempestvars.fOutputConvective, "output_convs"); \
	CommandLineBool(_tempestvars.fOutputZonalForce, "output_drag"); \
	CommandLineBool(_tempestvars.fNoReferenceState, "norefstate"); \
	CommandLineBool(_tempestvars.fNoTracers, "notracers"); \
	CommandLineBool(_tempestvars.fNoHyperviscosity, "nohypervis"); \
	CommandLineInt(_tempestvars.nHyperviscosityOrder, "hypervisorder", 4); \
	CommandLineDouble(_tempestvars.dNuScalar, "nu", 1.0e15); \
	CommandLineDouble(_tempestvars.dNuDiv, "nud", 1.0e15); \
	CommandLineDouble(_tempestvars.dNuVort, "nuv", 1.0e15); \
	CommandLineDouble(_tempestvars.dInstepNuDiv, "inud", 0.0); \
	CommandLineBool(_tempestvars.fExplicitVertical, "explicitvertical"); \
	CommandLineStringD(_tempestvars.strVerticalStaggering, "vstagger", "LOR", "(LEV | INT | LOR | CPH)"); \
	CommandLineStringD(_tempestvars.strVerticalDiscretization, "vdisc", "FE", "(FE | FV)"); \
	CommandLineBool(_tempestvars.fForceMassFluxOnLevels, "vmassfluxlevels"); \
	CommandLineString(_tempestvars.strVerticalStretch, "vstretch", "uniform"); \
	CommandLineInt(_tempestvars.nVerticalHyperdiffOrder, "vhypervisorder", 0); \
	CommandLineString(_tempestvars.strTimestepScheme, "timescheme", "strang"); \
	CommandLineStringD(_tempestvars.strHorizontalDynamics, "hmethod", "V1", "(V1 | V2 | SPEX)"); \
	CommandLineStringD(_tempestvars.strVerticalDynamics, "vmethod", "V1", "(V1 | V2 | SCHUR | NONE)");

///////////////////////////////////////////////////////////////////////////////

#define BeginTempestCommandLine(TestCaseName) \
	_TempestCommandLineVariables _tempestvars; \
	BeginCommandLine() \
	_TempestDefineCommandLineDefault(TestCaseName)

#define EndTempestCommandLine(argv) \
	EndCommandLine(argv)

#define SetDefaultResolution(default_resolution) \
	CommandLineInt(_tempestvars.nResolutionX, "resolution", default_resolution);

#define SetDefaultResolutionX(default_resolution_x) \
	CommandLineInt(_tempestvars.nResolutionX, "resx", default_resolution_x);

#define SetDefaultResolutionY(default_resolution_y) \
	CommandLineInt(_tempestvars.nResolutionY, "resy", default_resolution_y);

#define SetDefaultLevels(default_levels) \
	CommandLineInt(_tempestvars.nLevels, "levels", default_levels);

#define SetDefaultOutputDeltaT(default_output_dt) \
	CommandLineDeltaTime(_tempestvars.timeOutputDeltaT, "outputtime", default_output_dt);

#define SetDefaultDeltaT(default_dt) \
	CommandLineDeltaTime(_tempestvars.timeDeltaT, "dt", default_dt);

#define SetDefaultEndTime(default_endtime) \
	CommandLineFixedTime(_tempestvars.timeEndTime, "endtime", default_endtime);

#define SetDefaultHorizontalOrder(default_order) \
	CommandLineInt(_tempestvars.nHorizontalOrder, "order", default_order)

#define SetDefaultVerticalOrder(default_vertorder) \
	CommandLineInt(_tempestvars.nVerticalOrder, "vertorder", default_vertorder)

///////////////////////////////////////////////////////////////////////////////

void _TempestSetupMethodOfLines(
	Model & model,
	_TempestCommandLineVariables & vars
) {
	// Set the timestep scheme
	AnnounceStartBlock("Initializing time scheme");

	STLStringHelper::ToLower(vars.strTimestepScheme);
	if (vars.strTimestepScheme == "strang") {
		model.SetTimestepScheme(
			new TimestepSchemeStrang(model));

	} else if (vars.strTimestepScheme == "strang/fe") {
		model.SetTimestepScheme(
			new TimestepSchemeStrang(
				model, 0.0, TimestepSchemeStrang::ForwardEuler));

	} else if (vars.strTimestepScheme == "strang/rk4") {
		model.SetTimestepScheme(
			new TimestepSchemeStrang(
				model, 0.0, TimestepSchemeStrang::RungeKutta4));

	} else if (vars.strTimestepScheme == "strang/rk3") {
		model.SetTimestepScheme(
			new TimestepSchemeStrang(
				model, 0.0, TimestepSchemeStrang::RungeKuttaSSP3));

	} else if (vars.strTimestepScheme == "strang/kgu35") {
		model.SetTimestepScheme(
			new TimestepSchemeStrang(
				model, 0.0, TimestepSchemeStrang::KinnmarkGrayUllrich35));

	} else if (vars.strTimestepScheme == "strang/ssprk53") {
		model.SetTimestepScheme(
			new TimestepSchemeStrang(
				model, 0.0, TimestepSchemeStrang::RungeKuttaSSPRK53));

	} else if (vars.strTimestepScheme == "erk") {
		model.SetTimestepScheme(
			new TimestepSchemeERK(model));

	} else if (vars.strTimestepScheme == "erk/fe") {
		model.SetTimestepScheme(
			new TimestepSchemeERK(
				model, TimestepSchemeERK::ForwardEuler));

	} else if (vars.strTimestepScheme == "erk/rk4") {
		model.SetTimestepScheme(
			new TimestepSchemeERK(
				model, TimestepSchemeERK::RungeKutta4));

	} else if (vars.strTimestepScheme == "erk/rk3") {
		model.SetTimestepScheme(
			new TimestepSchemeERK(
				model, TimestepSchemeERK::RungeKuttaSSP3));

	} else if (vars.strTimestepScheme == "erk/kgu35") {
		model.SetTimestepScheme(
			new TimestepSchemeERK(
				model, TimestepSchemeERK::KinnmarkGrayUllrich35));

	} else if (vars.strTimestepScheme == "erk/ssprk53") {
		model.SetTimestepScheme(
			new TimestepSchemeERK(
				model, TimestepSchemeERK::RungeKuttaSSPRK53));

	} else if (vars.strTimestepScheme == "ars222") {
		model.SetTimestepScheme(
			new TimestepSchemeARS222(model));

	} else if (vars.strTimestepScheme == "ars232") {
		model.SetTimestepScheme(
			new TimestepSchemeARS232(model));

	} else if (vars.strTimestepScheme == "ark232") {
		model.SetTimestepScheme(
			new TimestepSchemeARK232(model));

		} else if (vars.strTimestepScheme == "gark2") {
		model.SetTimestepScheme(
			new TimestepSchemeGARK2(model));

	} else if (vars.strTimestepScheme == "ars343") {
		model.SetTimestepScheme(
			new TimestepSchemeARS343(model));

	} else if (vars.strTimestepScheme == "ars343b") {
		model.SetTimestepScheme(
			new TimestepSchemeARS343b(model));

	} else if (vars.strTimestepScheme == "ars443") {
		model.SetTimestepScheme(
			new TimestepSchemeARS443(model));

	} else if (vars.strTimestepScheme == "ssp3_332") {
		model.SetTimestepScheme(
			new TimestepSchemeSSP3332(model));

	} else if (vars.strTimestepScheme == "spex") {
		model.SetTimestepScheme(
			new TimestepSchemeSplitExp(model));

	} else {
		_EXCEPTIONT("Invalid timescheme: Expected "
			"\"Strang\", \"ERK\", \"ARS222\", \"ARS232\", \"ARK232\", "
			"\"ARS343\", \"ARS443\", \"SSP3_332\"");
	}
	AnnounceEndBlock("Done");

	// Set the horizontal dynamics
	AnnounceStartBlock("Initializing horizontal dynamics");

	if (vars.fNoHyperviscosity) {
		vars.dNuScalar = 0.0;
		vars.dNuDiv = 0.0;
		vars.dNuVort = 0.0;
	}

	STLStringHelper::ToLower(vars.strHorizontalDynamics);

	if (vars.strHorizontalDynamics == "v1") {
		model.SetHorizontalDynamics(
			new HorizontalDynamicsFEM(
				model,
				vars.nHorizontalOrder,
				vars.nHyperviscosityOrder,
				vars.dNuScalar,
				vars.dNuDiv,
				vars.dNuVort,
				vars.dInstepNuDiv));

	} else if (vars.strHorizontalDynamics == "v2") {
		model.SetHorizontalDynamics(
			new HorizontalDynamicsFEMV2(
				model,
				vars.nHorizontalOrder,
				vars.nHyperviscosityOrder,
				vars.dNuScalar,
				vars.dNuDiv,
				vars.dNuVort,
				vars.dInstepNuDiv));

	} else if (vars.strHorizontalDynamics == "spex") {
		model.SetHorizontalDynamics(
			new SplitExplicitDynamics(
				model,
				vars.nHorizontalOrder,
				vars.nHyperviscosityOrder,
				vars.dNuScalar,
				vars.dNuDiv,
				vars.dNuVort,
				vars.dInstepNuDiv));

	} else if (vars.strHorizontalDynamics == "hs") {
		model.SetHorizontalDynamics(
			new HighSpeedDynamics(
				model,
				vars.nHorizontalOrder,
				vars.nHyperviscosityOrder,
				vars.dNuScalar,
				vars.dNuDiv,
				vars.dNuVort,
				vars.dInstepNuDiv));

	} else {
		_EXCEPTIONT("Invalid --hmethod");
	}

	AnnounceEndBlock("Done");

	// Vertical staggering
	STLStringHelper::ToLower(vars.strVerticalStaggering);

	// Vertical method
	STLStringHelper::ToLower(vars.strVerticalDynamics);

	AnnounceStartBlock("Initializing vertical dynamics");
	if (vars.nLevels == 1) {
		model.SetVerticalDynamics(
			new VerticalDynamicsStub(model));

	} else if (vars.strVerticalDynamics == "v1") {
		model.SetVerticalDynamics(
			new VerticalDynamicsFEM(
				model,
				vars.nHorizontalOrder,
				vars.nVerticalOrder,
				vars.nVerticalHyperdiffOrder,
				vars.fExplicitVertical,
				!vars.fNoReferenceState,
				vars.fForceMassFluxOnLevels));

	} else if (vars.strVerticalDynamics == "v2") {
		model.SetVerticalDynamics(
			new VerticalDynamicsFEMV2(
				model,
				vars.nHorizontalOrder,
				vars.nVerticalOrder,
				vars.nVerticalHyperdiffOrder,
				vars.fExplicitVertical,
				!vars.fNoReferenceState,
				vars.fForceMassFluxOnLevels));

	} else if (vars.strVerticalDynamics == "schur") {
		model.SetVerticalDynamics(
			new VerticalDynamicsSchur(
				model,
				vars.nHorizontalOrder,
				vars.nVerticalOrder,
				vars.nVerticalHyperdiffOrder,
				vars.fExplicitVertical,
				!vars.fNoReferenceState,
				vars.fForceMassFluxOnLevels));

	} else if (vars.strVerticalDynamics == "none") {
		model.SetVerticalDynamics(
			new VerticalDynamicsStub(
				model));

	} else {
		_EXCEPTIONT("Invalid --vmethod");
	}

	AnnounceEndBlock("Done");
}

///////////////////////////////////////////////////////////////////////////////

void _TempestSetupOutputManagers(
	Model & model,
	_TempestCommandLineVariables & vars
) {
	// Set the reference output manager for the model
	if (!vars.fNoOutput) {
		AnnounceStartBlock("Creating reference output manager");
		OutputManagerReference * pOutmanRef =
			new OutputManagerReference(
				*(model.GetGrid()),
				vars.timeOutputDeltaT,
				vars.strOutputDir,
				vars.strOutputPrefix,
				vars.nOutputsPerFile,
				vars.nOutputResX,
				vars.nOutputResY,
				vars.nOutputResZ,
				false,
				false);

		if (vars.fOutputVorticity) {
			pOutmanRef->OutputVorticity();
		}
		if (vars.fOutputDivergence) {
			pOutmanRef->OutputDivergence();
		}
		if (vars.fOutputTemperature) {
			pOutmanRef->OutputTemperature();
		}
		if (vars.fOutputSurfacePressure) {
			pOutmanRef->OutputSurfacePressure();
		}
		if (vars.fOutputRichardson) {
			pOutmanRef->OutputRichardson();
		}
		if (vars.fOutputConvective) {
			pOutmanRef->OutputConvective();
		}
		if (vars.fOutputZonalForce) {
                        pOutmanRef->OutputZonalForce();
		}

		model.AttachOutputManager(pOutmanRef);
		AnnounceEndBlock("Done");
	}

	// Set the composite output manager for the model
	if (!vars.timeOutputRestartDeltaT.IsZero()) {
		AnnounceStartBlock("Creating composite output manager");
		model.AttachOutputManager(
			new OutputManagerComposite(
				*(model.GetGrid()),
				vars.timeOutputRestartDeltaT,
				vars.strOutputDir,
				vars.strOutputPrefix));
		AnnounceEndBlock("Done");
	}

	// Set the checksum output manager for the model
	AnnounceStartBlock("Creating checksum output manager");
	model.AttachOutputManager(
		new OutputManagerChecksum(
			*(model.GetGrid()),
			vars.timeOutputDeltaT));
	AnnounceEndBlock("Done");
}

///////////////////////////////////////////////////////////////////////////////

void _TempestSetupCubedSphereModel(
	Model & model,
	_TempestCommandLineVariables & vars
) {
	// Set the time step size and end time
	model.SetDeltaT(vars.timeDeltaT);
	model.SetEndTime(vars.timeEndTime);

	// Setup Method of Lines
	_TempestSetupMethodOfLines(model, vars);

	// Get the vertical discretization from --vdisc
	Grid::VerticalDiscretization eVerticalDiscretization;
	STLStringHelper::ToLower(vars.strVerticalDiscretization);
	if (vars.strVerticalDiscretization == "fe") {
		eVerticalDiscretization = Grid::VerticalDiscretization_FiniteElement;

	} else if (vars.strVerticalDiscretization == "fv") {
		eVerticalDiscretization = Grid::VerticalDiscretization_FiniteVolume;

	} else {
		_EXCEPTIONT("Invalid value for --vdisc");
	}

	// Get the vertical staggering from --vstagger
	Grid::VerticalStaggering eVerticalStaggering;
	STLStringHelper::ToLower(vars.strVerticalStaggering);
	if (vars.strVerticalStaggering == "lev") {
		eVerticalStaggering = Grid::VerticalStaggering_Levels;

	} else if (vars.strVerticalStaggering == "int") {
		eVerticalStaggering = Grid::VerticalStaggering_Interfaces;

	} else if (vars.strVerticalStaggering == "lor") {
		eVerticalStaggering = Grid::VerticalStaggering_Lorenz;

	} else if (vars.strVerticalStaggering == "cph") {
		eVerticalStaggering = Grid::VerticalStaggering_CharneyPhillips;

	} else {
		_EXCEPTIONT("Invalid value for --vstagger");
	}

	// Construct the Grid
	if (vars.strRestartFile == "") {
		AnnounceStartBlock("Constructing grid");

		// Maximum number of patches currently equals communicator size
		int nCommSize = 1;
#ifdef TEMPEST_MPIOMP
		MPI_Comm_size(MPI_COMM_WORLD, &nCommSize);
#endif

		if (nCommSize < 6) {
			nCommSize = 6;
		}

		GridCSGLL * pGrid = new GridCSGLL(model);

		pGrid->DefineParameters();

		pGrid->SetParameters(
			vars.nLevels,
			nCommSize,
			vars.nResolutionX,
			4,
			vars.nHorizontalOrder,
			vars.nVerticalOrder,
			eVerticalDiscretization,
			eVerticalStaggering);

		pGrid->InitializeDataLocal();

		// Set the vertical stretching function
		STLStringHelper::ToLower(vars.strVerticalStretch);
		if (vars.strVerticalStretch == "uniform") {

		} else if (vars.strVerticalStretch == "cubic") {
			pGrid->SetVerticalStretchFunction(
				new VerticalStretchCubic);

		} else if (vars.strVerticalStretch == "pwlinear") {
			pGrid->SetVerticalStretchFunction(
				new VerticalStretchPiecewiseLinear);

		} else {
			_EXCEPTIONT("Invalid value for --vstretch");
		}

		// Vertical type
		STLStringHelper::ToLower(vars.strVerticalDiscretization);

		// Set the Model Grid
		model.SetGrid(pGrid);

	// Set the Grid from Restart file
	} else {
		AnnounceStartBlock("Constructing Grid from restart file");

		Grid * pGrid = new GridCSGLL(model);

		pGrid->DefineParameters();

		model.SetGridFromRestartFile(pGrid, vars.strRestartFile);
	}

	AnnounceEndBlock("Done");

	// Setup OutputManagers
	_TempestSetupOutputManagers(model, vars);
}

///////////////////////////////////////////////////////////////////////////////

void _TempestSetupCartesianModel(
	Model & model,
	double dGDim[],
	double dRefLat,
	int iLatBC[],
	bool fCartesianXZ,
	_TempestCommandLineVariables & vars
) {
	// Set the time step size and end time
	model.SetDeltaT(vars.timeDeltaT);
	model.SetEndTime(vars.timeEndTime);

	// Setup Method of Lines
	_TempestSetupMethodOfLines(model, vars);

	// Get the vertical discretization from --vdisc
	Grid::VerticalDiscretization eVerticalDiscretization;
	STLStringHelper::ToLower(vars.strVerticalDiscretization);
	if (vars.strVerticalDiscretization == "fe") {
		eVerticalDiscretization = Grid::VerticalDiscretization_FiniteElement;

	} else if (vars.strVerticalDiscretization == "fv") {
		eVerticalDiscretization = Grid::VerticalDiscretization_FiniteVolume;

	} else {
		_EXCEPTIONT("Invalid value for --vdisc");
	}

	// Get the vertical staggering from --vstagger
	Grid::VerticalStaggering eVerticalStaggering;
	STLStringHelper::ToLower(vars.strVerticalStaggering);
	if (vars.strVerticalStaggering == "lev") {
		eVerticalStaggering = Grid::VerticalStaggering_Levels;

	} else if (vars.strVerticalStaggering == "int") {
		eVerticalStaggering = Grid::VerticalStaggering_Interfaces;

	} else if (vars.strVerticalStaggering == "lor") {
		eVerticalStaggering = Grid::VerticalStaggering_Lorenz;

	} else if (vars.strVerticalStaggering == "cph") {
		eVerticalStaggering = Grid::VerticalStaggering_CharneyPhillips;

	} else {
		_EXCEPTIONT("Invalid value for --vstagger");
	}

	// Construct the Grid
	if (vars.strRestartFile == "") {
		AnnounceStartBlock("Constructing grid");

		// Maximum number of patches currently equals communicator size
		int nCommSize = 1;
#ifdef TEMPEST_MPIOMP
		MPI_Comm_size(MPI_COMM_WORLD, &nCommSize);
#endif

		GridCartesianGLL * pGrid = new GridCartesianGLL(model);

		pGrid->DefineParameters();

		pGrid->SetParameters(
			vars.nLevels,
			nCommSize,
			vars.nResolutionX,
			vars.nResolutionY,
			4,
			vars.nHorizontalOrder,
			vars.nVerticalOrder,
			dGDim,
			dRefLat,
			iLatBC,
			fCartesianXZ,
			eVerticalDiscretization,
			eVerticalStaggering);

		pGrid->InitializeDataLocal();

		// Set the vertical stretching function
		STLStringHelper::ToLower(vars.strHorizontalDynamics);
		if (vars.strVerticalStretch == "uniform") {

		} else if (vars.strVerticalStretch == "cubic") {
			pGrid->SetVerticalStretchFunction(
				new VerticalStretchCubic);

		} else if (vars.strVerticalStretch == "pwlinear") {
			pGrid->SetVerticalStretchFunction(
				new VerticalStretchPiecewiseLinear);

		} else {
			_EXCEPTIONT("Invalid value for --vstretch");

		}

		// Set the Model Grid
		model.SetGrid(pGrid);

	// Set the Grid from Restart file
	} else {
		AnnounceStartBlock("Constructing Grid from restart file");

		Grid * pGrid = new GridCartesianGLL(model);

		pGrid->DefineParameters();

		model.SetGridFromRestartFile(pGrid, vars.strRestartFile);
	}

	AnnounceEndBlock("Done");

	// Setup OutputManagers
	_TempestSetupOutputManagers(model, vars);
}

///////////////////////////////////////////////////////////////////////////////

#define TempestSetupCubedSphereModel(model) \
	_TempestSetupCubedSphereModel(model, _tempestvars);

#define TempestSetupCartesianModel(model, dimensions, latitude, lateralBC, CartesianXZ) \
	_TempestSetupCartesianModel(model, dimensions, latitude, lateralBC, CartesianXZ, _tempestvars);

///////////////////////////////////////////////////////////////////////////////

void TempestInitialize(int * argc, char*** argv) {

#ifdef TEMPEST_PETSC
	// Initialize PetSc
	PetscInitialize(argc, argv, NULL, NULL);
#endif
#ifdef TEMPEST_MPIOMP
	// Initialize MPI
	MPI_Init(argc, argv);
#endif

	AnnounceOnlyOutputOnRankZero();
}

///////////////////////////////////////////////////////////////////////////////

void TempestAbort() {

#ifdef TEMPEST_MPIOMP
	// Abort
	MPI_Abort(MPI_COMM_WORLD, 1);
#endif

}

///////////////////////////////////////////////////////////////////////////////

void TempestDeinitialize() {

#ifdef TEMPEST_PETSC
	// Finalize PetSc
	PetscFinalize();
#endif
#ifdef TEMPEST_MPIOMP
	// Finalize MPI
	MPI_Finalize();
#endif

}

///////////////////////////////////////////////////////////////////////////////

#endif
