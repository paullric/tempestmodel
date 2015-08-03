///////////////////////////////////////////////////////////////////////////////
///
///	\file    OutputManagerComposite.cpp
///	\author  Paul Ullrich
///	\version August 11, 2010
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

#include "OutputManagerComposite.h"

#include "Model.h"
#include "Grid.h"
#include "ConsolidationStatus.h"

#include "TimeObj.h"
#include "Announce.h"

#include "mpi.h"

#include <iostream>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <sys/stat.h>

///////////////////////////////////////////////////////////////////////////////

OutputManagerComposite::OutputManagerComposite(
	Grid & grid,
	const Time & timeOutputFrequency,
	std::string strOutputDir,
	std::string strOutputFormat,
	std::string strRestartFile
) :
	OutputManager(
		grid,
		timeOutputFrequency,
		strOutputDir,
		strOutputFormat,
		1),
	m_pActiveNcOutput(NULL),
	m_varZs(NULL),
	m_varRayleighStrengthNode(NULL),
	m_varRayleighStrengthREdge(NULL),
	m_strRestartFile(strRestartFile)
{
}

///////////////////////////////////////////////////////////////////////////////

OutputManagerComposite::~OutputManagerComposite() {
	CloseFile();
}

///////////////////////////////////////////////////////////////////////////////

bool OutputManagerComposite::OpenFile(
	const std::string & strFileName
) {
	// Determine processor rank
	int nRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

	// The active model
	const Model & model = m_grid.GetModel();

	// Open new NetCDF file
	NcVar * varZs;
	NcVar * varRayleighStrengthNode;
	NcVar * varRayleighStrengthREdge;

	if (nRank == 0) {

		// Allocate receive buffer
		m_vecRecvBuffer.Allocate(m_grid.GetMaxDegreesOfFreedom());

		// Check for existing NetCDF file
		if (m_pActiveNcOutput != NULL) {
			_EXCEPTIONT("NetCDF file already open");
		}

		// Open new NetCDF file
		std::string strNcFileName = strFileName + ".restart.nc";
		m_pActiveNcOutput = new NcFile(strNcFileName.c_str(), NcFile::Replace);
		if (m_pActiveNcOutput == NULL) {
			_EXCEPTIONT("Error opening NetCDF file");
		}

		// Create nodal index dimension
		NcDim * dimNodeIndex2D =
			m_pActiveNcOutput->add_dim(
				"node_index_2d", m_grid.GetTotalNodeCount2D());

		NcDim * dimNodeIndex =
			m_pActiveNcOutput->add_dim(
				"node_index", m_grid.GetTotalNodeCount(DataLocation_Node));

		NcDim * dimREdgeIndex =
			m_pActiveNcOutput->add_dim(
				"redge_index", m_grid.GetTotalNodeCount(DataLocation_REdge));

		// Output start time and current time
		m_pActiveNcOutput->add_att("start_time",
			model.GetStartTime().ToLongString().c_str());

		// Output equation set
		const EquationSet & eqn = model.GetEquationSet();

		m_pActiveNcOutput->add_att("equation_set", eqn.GetName().c_str());

		// Output the grid
		m_grid.ToFile(*m_pActiveNcOutput);

		// Get the patch index dimension
		NcDim * dimPatchIndex = m_pActiveNcOutput->get_dim("patch_index");
		if (dimPatchIndex == NULL) {
			_EXCEPTIONT("Dimension \'patch_index\' not found");
		}

		// Create variables in NetCDF output file
		for (int c = 0; c < eqn.GetComponents(); c++) {

			// State variables (only store on relevant level)
			if (m_grid.GetVarLocation(c) == DataLocation_Node) {
				m_vecStateVar.push_back(
					m_pActiveNcOutput->add_var(
						eqn.GetComponentShortName(c).c_str(),
						ncDouble, dimNodeIndex));

			} else if (m_grid.GetVarLocation(c) == DataLocation_REdge) {
				m_vecStateVar.push_back(
					m_pActiveNcOutput->add_var(
						eqn.GetComponentShortName(c).c_str(),
						ncDouble, dimREdgeIndex));

			} else {
				_EXCEPTIONT("(UNIMPLEMENTED) Invalid DataLocation");
			}

			// Reference state variables (store on all levels)
			if (m_grid.HasReferenceState()) {
				std::string strRefVarNameNode =
					eqn.GetComponentShortName(c) + "_RefNode";

				m_vecRefStateVarNode.push_back(
					m_pActiveNcOutput->add_var(
						strRefVarNameNode.c_str(),
						ncDouble, dimNodeIndex));

				std::string strRefVarNameREdge =
					eqn.GetComponentShortName(c) + "_RefREdge";

				m_vecRefStateVarREdge.push_back(
					m_pActiveNcOutput->add_var(
						strRefVarNameREdge.c_str(),
						ncDouble, dimREdgeIndex));
			}
		}

		for (int c = 0; c < eqn.GetTracers(); c++) {
			m_vecTracersVar.push_back(
				m_pActiveNcOutput->add_var(
					eqn.GetTracerShortName(c).c_str(),
					ncDouble, dimNodeIndex));
		}

		// Output topography for each patch
		m_varZs = m_pActiveNcOutput->add_var("ZS", ncDouble, dimNodeIndex2D);

		// Output Rayleigh strength for each patch
		m_varRayleighStrengthNode =
			m_pActiveNcOutput->add_var(
				"Rayleigh_Node", ncDouble, dimNodeIndex);

		m_varRayleighStrengthREdge =
			m_pActiveNcOutput->add_var(
				"Rayleigh_REdge", ncDouble, dimREdgeIndex);
	}

	// Wait for all processes to complete
	MPI_Barrier(MPI_COMM_WORLD);

	return true;
}	

///////////////////////////////////////////////////////////////////////////////

void OutputManagerComposite::CloseFile() {
	if (m_pActiveNcOutput != NULL) {
		delete(m_pActiveNcOutput);
		m_pActiveNcOutput = NULL;

		m_vecStateVar.clear();
		m_vecRefStateVarNode.clear();
		m_vecRefStateVarREdge.clear();
		m_vecTracersVar.clear();

		m_vecRecvBuffer.Deallocate();
	}
}

///////////////////////////////////////////////////////////////////////////////

void OutputManagerComposite::Output(
	const Time & time
) {
	// Check for open file
	if (!IsFileOpen()) {
		_EXCEPTIONT("No file available for output");
	}

	// Verify that only one output has been performed
	if (m_ixOutputTime != 0) {
		_EXCEPTIONT("Only one Composite output allowed per file");
	}

	// Determine processor rank
	int nRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

	// Equation set
	const EquationSet & eqn = m_grid.GetModel().GetEquationSet();

	// Output start time and current time
	if (nRank == 0) {
		m_pActiveNcOutput->add_att(
			"current_time", time.ToLongString().c_str());
	}

	// Begin data consolidation
	std::vector<DataTypeLocationPair> vecDataTypes;
	vecDataTypes.push_back(DataType_Topography);
	vecDataTypes.push_back(
		DataTypeLocationPair(DataType_RayleighStrength, DataLocation_Node));
	vecDataTypes.push_back(
		DataTypeLocationPair(DataType_RayleighStrength, DataLocation_REdge));

	if (m_grid.GetVarsAtLocation(DataLocation_Node) != 0) {
		vecDataTypes.push_back(
			DataTypeLocationPair(DataType_State, DataLocation_Node));
	}
	if (m_grid.GetVarsAtLocation(DataLocation_REdge) != 0) {
		vecDataTypes.push_back(
			DataTypeLocationPair(DataType_State, DataLocation_REdge));
	}

	if (m_grid.HasReferenceState()) {
		vecDataTypes.push_back(
			DataTypeLocationPair(DataType_RefState, DataLocation_Node));
		vecDataTypes.push_back(
			DataTypeLocationPair(DataType_RefState, DataLocation_REdge));
	}

	if (eqn.GetTracers() != 0) {
		vecDataTypes.push_back(DataType_Tracers);
	}

	ConsolidationStatus status(m_grid, vecDataTypes);

	m_grid.ConsolidateDataToRoot(status);

	// Receive all data objects from neighbors
	while ((nRank == 0) && (!status.Done())) {
		int nRecvCount;
		int ixRecvPatch;
		DataType eRecvDataType;
		DataLocation eRecvDataLocation;

		m_grid.ConsolidateDataAtRoot(
			status,
			m_vecRecvBuffer,
			nRecvCount,
			ixRecvPatch,
			eRecvDataType,
			eRecvDataLocation);

		// Store topography data
		if (eRecvDataType == DataType_Topography) {
			int nPatchIx = m_grid.GetCumulativePatch2DNodeIndex(ixRecvPatch);
			m_varZs->set_cur(nPatchIx);
			m_varZs->put(&(m_vecRecvBuffer[0]), nRecvCount);

		// Store Rayleigh strength data at nodes and edges
		} else if (eRecvDataType == DataType_RayleighStrength) {
			if (eRecvDataLocation == DataLocation_Node) {
				if (nRecvCount % m_grid.GetRElements() != 0) {
					_EXCEPTIONT("Invalid message length");
				}

			} else if (eRecvDataLocation == DataLocation_REdge) {
				if (nRecvCount % (m_grid.GetRElements()+1) != 0) {
					_EXCEPTIONT("Invalid message length");
				}

			} else {
				_EXCEPTIONT("Invalid DataLocation");
			}

			int ixCumulative2DNode =
				m_grid.GetCumulativePatch2DNodeIndex(ixRecvPatch);

			int nComponentSize =
				m_grid.GetPatch(ixRecvPatch)->GetTotalNodeCount(
					eRecvDataLocation);

			if (eRecvDataLocation == DataLocation_Node) {
				m_varRayleighStrengthNode->set_cur(
					ixCumulative2DNode * m_grid.GetRElements());
				m_varRayleighStrengthNode->put(
					&(m_vecRecvBuffer[0]), nComponentSize);

			} else if (eRecvDataLocation == DataLocation_REdge) {
				m_varRayleighStrengthREdge->set_cur(
					ixCumulative2DNode * (m_grid.GetRElements()+1));
				m_varRayleighStrengthREdge->put(
					&(m_vecRecvBuffer[0]), nComponentSize);

			} else {
				_EXCEPTION();
			}

		// Store state variable data
		} else if (eRecvDataType == DataType_State) {
			if (eRecvDataLocation == DataLocation_Node) {
				if (nRecvCount % m_grid.GetRElements() != 0) {
					_EXCEPTIONT("Invalid message length");
				}

			} else if (eRecvDataLocation == DataLocation_REdge) {
				if (nRecvCount % (m_grid.GetRElements()+1) != 0) {
					_EXCEPTIONT("Invalid message length");
				}

			} else {
				_EXCEPTIONT("Invalid DataLocation");
			}

			int ixRecvPtr = 0;

			int ixCumulative2DNode =
				m_grid.GetCumulativePatch2DNodeIndex(ixRecvPatch);

			for (int c = 0; c < eqn.GetComponents(); c++) {

				// Advance the pointer in the receive buffer
				int nComponentSize =
					m_grid.GetPatch(ixRecvPatch)->GetTotalNodeCount(
						eRecvDataLocation);

				// Only output data at relevant location
				if (m_grid.GetVarLocation(c) == eRecvDataLocation) {
					int nRadialDegreesOfFreedom;
					if (eRecvDataLocation == DataLocation_Node) {
						nRadialDegreesOfFreedom = m_grid.GetRElements();
					} else if (eRecvDataLocation == DataLocation_REdge) {
						nRadialDegreesOfFreedom = m_grid.GetRElements()+1;
					} else {
						_EXCEPTION();
					}

					m_vecStateVar[c]->set_cur(
						ixCumulative2DNode * nRadialDegreesOfFreedom);
					m_vecStateVar[c]->put(
						&(m_vecRecvBuffer[ixRecvPtr]), nComponentSize);
				}

				ixRecvPtr += nComponentSize;
			}

		// Store reference state variable data
		} else if (eRecvDataType == DataType_RefState) {
			if (eRecvDataLocation == DataLocation_Node) {
				if (nRecvCount % m_grid.GetRElements() != 0) {
					_EXCEPTIONT("Invalid message length");
				}

			} else if (eRecvDataLocation == DataLocation_REdge) {
				if (nRecvCount % (m_grid.GetRElements()+1) != 0) {
					_EXCEPTIONT("Invalid message length");
				}

			} else {
				_EXCEPTIONT("Invalid DataLocation");
			}

			int ixRecvPtr = 0;

			int ixCumulative2DNode =
				m_grid.GetCumulativePatch2DNodeIndex(ixRecvPatch);

			for (int c = 0; c < eqn.GetComponents(); c++) {

				// Advance the pointer in the receive buffer
				int nComponentSize =
					m_grid.GetPatch(ixRecvPatch)->GetTotalNodeCount(
						eRecvDataLocation);

				// Write data to file
				int nRadialDegreesOfFreedom;
				if (eRecvDataLocation == DataLocation_Node) {
					nRadialDegreesOfFreedom = m_grid.GetRElements();

					m_vecRefStateVarNode[c]->set_cur(
						ixCumulative2DNode * nRadialDegreesOfFreedom);
					m_vecRefStateVarNode[c]->put(
						&(m_vecRecvBuffer[ixRecvPtr]), nComponentSize);

				} else if (eRecvDataLocation == DataLocation_REdge) {
					nRadialDegreesOfFreedom = m_grid.GetRElements()+1;

					m_vecRefStateVarREdge[c]->set_cur(
						ixCumulative2DNode * nRadialDegreesOfFreedom);
					m_vecRefStateVarREdge[c]->put(
						&(m_vecRecvBuffer[ixRecvPtr]), nComponentSize);

				} else {
					_EXCEPTION();
				}

				ixRecvPtr += nComponentSize;
			}

		// Store tracer variable data
		} else if (eRecvDataType == DataType_Tracers) {
			int nComponentSize =
				nRecvCount / eqn.GetTracers();

			if (nRecvCount % eqn.GetTracers() != 0) {
				_EXCEPTIONT("Invalid message length");
			}

			int nCumulative3DNodeIx =
				m_grid.GetCumulativePatch3DNodeIndex(ixRecvPatch);

			for (int c = 0; c < eqn.GetTracers(); c++) {
				m_vecTracersVar[c]->set_cur(nCumulative3DNodeIx);
				m_vecTracersVar[c]->put(
					&(m_vecRecvBuffer[nComponentSize * c]), nComponentSize);
			}
		}
	}

	// Barrier
	MPI_Barrier(MPI_COMM_WORLD);
}

///////////////////////////////////////////////////////////////////////////////

Time OutputManagerComposite::Input(
	const std::string & strFileName
) {
	// Set the flag indicating that output came from a restart file
	m_fFromRestartFile = true;

	// Determine processor rank
	int nRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

	// The active model
	const Model & model = m_grid.GetModel();

	// Open new NetCDF file
	NcFile * pNcFile = NULL;

	// Open NetCDF file
	pNcFile = new NcFile(strFileName.c_str(), NcFile::ReadOnly);
	if (pNcFile == NULL) {
		_EXCEPTIONT("Error opening NetCDF file");
	}
	if (! pNcFile->is_valid()) {
		_EXCEPTIONT("Error opening NetCDF file");
	}

	// Get the time
	NcAtt * attCurrentTime = pNcFile->get_att("current_time");
	if (attCurrentTime == NULL) {
		_EXCEPTIONT("Attribute \'current_time\' not found in restart file");
	}
	std::string strCurrentTime = attCurrentTime->as_string(0);

	Time timeCurrent(strCurrentTime);

	// Equation set
	const EquationSet & eqn = m_grid.GetModel().GetEquationSet();

	for (int n = 0; n < m_grid.GetActivePatchCount(); n++) {
		GridPatch * pPatch = m_grid.GetActivePatch(n);

		int nCumulative2DNodeIx =
			m_grid.GetCumulativePatch2DNodeIndex(pPatch->GetPatchIndex());

		int nPatchNodeCount = pPatch->GetTotalNodeCount2D();

		// Input topography here
		DataArray2D<double> & dataTopography =
			pPatch->GetTopography();

		NcVar * varTopography = pNcFile->get_var("ZS");
		if (varTopography == NULL) {
			_EXCEPTIONT("Cannot find variable \'ZS\' in file");
		}
		varTopography->set_cur(nCumulative2DNodeIx);

		varTopography->get(dataTopography[0], nPatchNodeCount);

		// Input Rayleigh strength here
		DataArray3D<double> & dataRayleighStrengthNode =
			pPatch->GetRayleighStrength(DataLocation_Node);
		DataArray3D<double> & dataRayleighStrengthREdge =
			pPatch->GetRayleighStrength(DataLocation_REdge);

		NcVar * varRayleighNode = pNcFile->get_var("Rayleigh_Node");
		NcVar * varRayleighREdge = pNcFile->get_var("Rayleigh_REdge");

		varRayleighNode->set_cur(
			nCumulative2DNodeIx * m_grid.GetRElements());
		varRayleighREdge->set_cur(
			nCumulative2DNodeIx * (m_grid.GetRElements()+1));

		varRayleighNode->get(
			dataRayleighStrengthNode[0][0],
			pPatch->GetTotalNodeCount(DataLocation_Node));

		varRayleighREdge->get(
			dataRayleighStrengthREdge[0][0],
			pPatch->GetTotalNodeCount(DataLocation_REdge));

		// Input state
		DataArray4D<double> & dataStateNode =
			pPatch->GetDataState(0, DataLocation_Node);
		DataArray4D<double> & dataStateREdge =
			pPatch->GetDataState(0, DataLocation_REdge);

		DataArray4D<double> & dataRefStateNode =
			pPatch->GetReferenceState(DataLocation_Node);
		DataArray4D<double> & dataRefStateREdge =
			pPatch->GetReferenceState(DataLocation_REdge);

		for (int c = 0; c < eqn.GetComponents(); c++) {

			DataLocation loc = m_grid.GetVarLocation(c);

			// Number of radial elements for this component
			int nRadialDegreesOfFreedom;
			if (loc == DataLocation_Node) {
				nRadialDegreesOfFreedom = m_grid.GetRElements();
			} else if (loc == DataLocation_REdge) {
				nRadialDegreesOfFreedom = m_grid.GetRElements()+1;
			} else {
				_EXCEPTION();
			}

			// Number of nodes on this patch for this component
			int nComponentSize = pPatch->GetTotalNodeCount(loc);

			// Load state variable
			std::string strComponentName = eqn.GetComponentShortName(c);

			NcVar * var = pNcFile->get_var(strComponentName.c_str());
			if (var == NULL) {
				_EXCEPTION1("Cannot find variable \'%s\' in file",
					strComponentName.c_str());
			}
			var->set_cur(nCumulative2DNodeIx * nRadialDegreesOfFreedom);

			if (loc == DataLocation_Node) {
				var->get(dataStateNode[c][0][0], nComponentSize);
			} else if (loc == DataLocation_REdge) {
				var->get(dataStateREdge[c][0][0], nComponentSize);
			} else {
				_EXCEPTION();
			}

			// Load reference state
			if (m_grid.HasReferenceState()) {

				// Load reference state on nodes
				std::string strRefComponentNameNode =
					strComponentName + "_RefNode";

				NcVar * varRefNode =
					pNcFile->get_var(strRefComponentNameNode.c_str());
				if (varRefNode == NULL) {
					_EXCEPTION1("Cannot find variable \'%s\' in file",
						strRefComponentNameNode.c_str());
				}
				varRefNode->set_cur(
					nCumulative2DNodeIx * m_grid.GetRElements());

				varRefNode->get(
					dataRefStateNode[c][0][0],
					pPatch->GetTotalNodeCount(DataLocation_Node));

				// Load reference state on radial edges
				std::string strRefComponentNameREdge =
					strComponentName + "_RefREdge";

				NcVar * varRefREdge =
					pNcFile->get_var(strRefComponentNameREdge.c_str());
				if (varRefREdge == NULL) {
					_EXCEPTION1("Cannot find variable \'%s\' in file",
						strRefComponentNameREdge.c_str());
				}
				varRefREdge->set_cur(
					nCumulative2DNodeIx * (m_grid.GetRElements()+1));

				varRefREdge->get(
					dataRefStateREdge[c][0][0],
					pPatch->GetTotalNodeCount(DataLocation_REdge));
			}
		}

		// Input tracers

	}

	// Close the file
	if (pNcFile != NULL) {
		delete pNcFile;
	}

	return timeCurrent;
}

///////////////////////////////////////////////////////////////////////////////

