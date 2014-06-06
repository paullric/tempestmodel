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
#include "Model.h"
#include "TimestepSchemeStrang.h"
#include "HorizontalDynamicsFEM.h"
#include "HorizontalDynamicsDG.h"
#include "VerticalDynamicsStub.h"
#include "VerticalDynamicsFEM.h"
#include "OutputManagerComposite.h"
#include "OutputManagerReference.h"
#include "OutputManagerChecksum.h"
#include "GridCSGLL.h"
#include "GridCartesianGLL.h"

#include "Announce.h"
#include "CommandLine.h"
#include "STLStringHelper.h"
#include <string>
#include "mpi.h"

#ifdef USE_PETSC
#include <petscsnes.h>
#endif

///////////////////////////////////////////////////////////////////////////////

struct _TempestCommandLineVariables {
	ModelParameters param;
	std::string strOutputDir;
	std::string strOutputPrefix;
	int nOutputsPerFile;
	double dOutputDeltaT;
	double dOutputRestartDeltaT;
	int nOutputResX;
	int nOutputResY;
	bool fOutputVorticity;
	bool fOutputDivergence;
	bool fOutputTemperature;
	bool fNoReferenceState;
	bool fNoTracers;
	bool fNoHyperviscosity;
	int nVerticalHyperdiffOrder;
	std::string strTimestepScheme;
	std::string strHorizontalDynamics;
	int nResolutionX;
	int nResolutionY;
	int nLevels;
	int nHorizontalOrder;
	int nVerticalOrder;
};

///////////////////////////////////////////////////////////////////////////////

#define _TempestDefineCommandLineDefault(TestCaseName) \
	CommandLineString(_tempestvars.strOutputDir, "output_dir", "out" TestCaseName); \
	CommandLineString(_tempestvars.strOutputPrefix, "output_prefix", "out"); \
	CommandLineString(_tempestvars.param.m_strRestartFile, "restart_file", ""); \
	CommandLineInt(_tempestvars.nOutputsPerFile, "output_perfile", -1); \
	CommandLineDouble(_tempestvars.dOutputRestartDeltaT, "output_restart_dt", 0.0); \
	CommandLineInt(_tempestvars.nOutputResX, "output_x", 360); \
	CommandLineInt(_tempestvars.nOutputResY, "output_y", 180); \
	CommandLineBool(_tempestvars.fOutputVorticity, "output_vort"); \
	CommandLineBool(_tempestvars.fOutputDivergence, "output_div"); \
	CommandLineBool(_tempestvars.fOutputTemperature, "output_temp"); \
	CommandLineBool(_tempestvars.fNoReferenceState, "norefstate"); \
	CommandLineBool(_tempestvars.fNoTracers, "notracers"); \
	CommandLineBool(_tempestvars.fNoHyperviscosity, "nohypervis"); \
	CommandLineInt(_tempestvars.nVerticalHyperdiffOrder, "verthypervisorder", 0); \
	CommandLineString(_tempestvars.strTimestepScheme, "timescheme", "strang"); \
	CommandLineStringD(_tempestvars.strHorizontalDynamics, "method", "SE", "(SE | DG)");

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

#define SetDefaultOutputTime(default_outputtime) \
	CommandLineDouble(_tempestvars.dOutputDeltaT, "outputtime", default_outputtime);

#define SetDefaultDeltaT(default_deltat) \
	CommandLineDouble(_tempestvars.param.m_dDeltaT, "dt", default_deltat);

#define SetDefaultEndTime(default_endtime) \
	CommandLineDouble(_tempestvars.param.m_dEndTime, "endtime", default_endtime);

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

	} else {
		_EXCEPTIONT("Invalid timescheme: Expected \"Strang\"");
	}
	AnnounceEndBlock("Done");

	// Set the horizontal dynamics
	AnnounceStartBlock("Initializing horizontal dynamics");

	STLStringHelper::ToLower(vars.strHorizontalDynamics);
	if (vars.strHorizontalDynamics == "se") {
		model.SetHorizontalDynamics(
			new HorizontalDynamicsFEM(
				model, vars.nHorizontalOrder, vars.fNoHyperviscosity));

	} else if (vars.strHorizontalDynamics == "dg") {
		model.SetHorizontalDynamics(
			new HorizontalDynamicsDG(
				model, vars.nHorizontalOrder, vars.fNoHyperviscosity));

	} else {
		_EXCEPTIONT("Invalid method: Expected \"SE\" or \"DG\"");
	}
	AnnounceEndBlock("Done");

	// Set the vertical dynamics
	AnnounceStartBlock("Initializing vertical dynamics");
	if (vars.nLevels == 1) {
		model.SetVerticalDynamics(
			new VerticalDynamicsStub(model));

	} else {
		model.SetVerticalDynamics(
			new VerticalDynamicsFEM(
				model,
				vars.nHorizontalOrder,
				vars.nVerticalOrder,
				vars.nVerticalHyperdiffOrder,
				false, // Implicit vertical
				!vars.fNoReferenceState));
	}

	AnnounceEndBlock("Done");
}

///////////////////////////////////////////////////////////////////////////////

void _TempestSetupOutputManagers(
	Model & model,
	_TempestCommandLineVariables & vars
) {
	// Set the reference output manager for the model
	AnnounceStartBlock("Creating reference output manager");
	OutputManagerReference * pOutmanRef =
		new OutputManagerReference(
			*(model.GetGrid()),
			vars.dOutputDeltaT,
			vars.strOutputDir,
			vars.strOutputPrefix,
			vars.nOutputsPerFile,
			vars.nOutputResX,
			vars.nOutputResY);

	if (vars.fOutputVorticity) {
		pOutmanRef->OutputVorticity();
	}
	if (vars.fOutputDivergence) {
		pOutmanRef->OutputDivergence();
	}
	if (vars.fOutputTemperature) {
		pOutmanRef->OutputTemperature();
	}

	model.AttachOutputManager(pOutmanRef);
	AnnounceEndBlock("Done");

	// Set the composite output manager for the model
	if (vars.dOutputRestartDeltaT != 0.0) {
		AnnounceStartBlock("Creating composite output manager");
		model.AttachOutputManager(
			new OutputManagerComposite(
				*(model.GetGrid()),
				vars.dOutputRestartDeltaT,
				vars.strOutputDir,
				vars.strOutputPrefix));
		AnnounceEndBlock("Done");
	}

	// Set the checksum output manager for the model
	AnnounceStartBlock("Creating checksum output manager");
	model.AttachOutputManager(
		new OutputManagerChecksum(
			*(model.GetGrid()),
			vars.dOutputDeltaT));
	AnnounceEndBlock("Done");
}

///////////////////////////////////////////////////////////////////////////////

void _TempestSetupCubedSphereModel(
	Model & model,
	_TempestCommandLineVariables & vars
) {
	// Set the parameters
	model.SetParameters(vars.param);

	// Setup Method of Lines
	_TempestSetupMethodOfLines(model, vars);

	// Construct the Grid
	AnnounceStartBlock("Constructing grid");
	model.SetGrid(
		new GridCSGLL(
			model,
			vars.nResolutionX,
			4,
			vars.nHorizontalOrder,
			vars.nVerticalOrder,
			vars.nLevels));

	AnnounceEndBlock("Done");

	// Setup OutputManagers
	_TempestSetupOutputManagers(model, vars);
}

///////////////////////////////////////////////////////////////////////////////

void _TempestSetupCartesianModel(
	Model & model,
	double dGDim[],
	_TempestCommandLineVariables & vars
) {
	// Set the parameters
	model.SetParameters(vars.param);

	// Setup Method of Lines
	_TempestSetupMethodOfLines(model, vars);

	// Set the model grid
	AnnounceStartBlock("Constructing grid");
	model.SetGrid(
		new GridCartesianGLL(
			model,
			vars.nResolutionX,
			vars.nResolutionY,
			4,
			vars.nHorizontalOrder,
			vars.nVerticalOrder,
			vars.nLevels,
			dGDim));

	AnnounceEndBlock("Done");

	// Setup OutputManagers
	_TempestSetupOutputManagers(model, vars);
}

///////////////////////////////////////////////////////////////////////////////

#define TempestSetupCubedSphereModel(model) \
	_TempestSetupCubedSphereModel(model, _tempestvars);

#define TempestSetupCartesianModel(model, dimensions) \
	_TempestSetupCartesianModel(model, dimensions, _tempestvars);

///////////////////////////////////////////////////////////////////////////////

void TempestInitialize(int * argc, char*** argv) {

#ifdef USE_PETSC
	// Initialize PetSc
	PetscInitialize(argc, argv, NULL, NULL);
#else
	// Initialize MPI
	MPI_Init(argc, argv);
#endif

}

///////////////////////////////////////////////////////////////////////////////

void TempestDeinitialize() {

#ifdef USE_PETSC
	// Finalize PetSc
	PetscFinalize();
#else
	// Finalize MPI
	MPI_Finalize();
#endif

}

///////////////////////////////////////////////////////////////////////////////

#endif

