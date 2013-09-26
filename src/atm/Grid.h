///////////////////////////////////////////////////////////////////////////////
///
///	\file    Grid.h
///	\author  Paul Ullrich
///	\version February 25, 2013
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

#ifndef _GRID_H_
#define _GRID_H_

#include "GridPatch.h"
#include "ChecksumType.h"

#include "mpi.h"

#include <string>
#include <vector>
#include <map>

///////////////////////////////////////////////////////////////////////////////

class Model;
class TestCase;
class ConsolidationStatus;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Atmospheric model grid data.  Container for GridPatch objects.
///	</summary>
class Grid {

friend class Model;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	Grid(
		const Model & model,
		int nABaseResolution,
		int nBBaseResolution,
		int nRefinementRatio,
		int nRElements
	);

	///	<summary>
	///		Virtual destructor.
	///	</summary>
	virtual ~Grid();

protected:
	///	<summary>
	///		Initialize grid patches.
	///	</summary>
	virtual void Initialize() {
		m_fInitialized = true;
	}

	///	<summary>
	///		Initialize the vertical coordinate.
	///	</summary>
	virtual void InitializeVerticalCoordinate(
		const GridSpacing & aGridSpacing
	);

	///	<summary>
	///		Evaluate the test case.
	///	</summary>
	void EvaluateTestCase(
		const TestCase & test,
		double dTime = 0.0,
		int iDataIndex = 0
	);

public:
	///	<summary>
	///		Perform checksum calculation on all state variables.
	///	</summary>
	void Checksum(
		DataType eDataType,
		DataVector<double> & dChecksums,
		int iDataIndex = 0,
		ChecksumType eChecksumType = ChecksumType_Sum
	) const;

public:
	///	<summary>
	///		Exchange data between processors.
	///	</summary>
	void Exchange(
		DataType eDataType,
		int iDataIndex
	);

public:
	///	<summary>
	///		Get the total number of patches on the grid.
	///	</summary>
	int GetPatchCount() const {
		return m_vecGridPatches.size();
	}

	///	<summary>
	///		Get the total number of active patches on this processor.
	///	</summary>
	int GetActivePatchCount() const {
		return m_vecActiveGridPatches.size();
	}

	///	<summary>
	///		Get the total number of nodes in the largest patch, returned to
	///		the root node.
	///	</summary>
	int GetLargestGridPatchNodes() const;

	///	<summary>
	///		Get the total count of nodes over the global grid, returned to
	///		the root node.
	///	</summary>
	int GetTotalNodeCount() const;

	///	<summary>
	///		Get the maximum number of degrees of freedom required for a
	///		consolidate operation, returned to the root node.
	///	</summary>
	int GetMaximumDegreesOfFreedom() const;

public:
	///	<summary>
	///		Receive a data object on the root process.
	///	</summary>
	void ConsolidateDataAtRoot(
		ConsolidationStatus & status,
		DataVector<double> & dataRecvBuffer,
		int & nRecvCount,
		int & ixRecvPatch,
		DataType & eRecvDataType
	) const;

	///	<summary>
	///		Send data to the root process and block on receipt.
	///	</summary>
	void ConsolidateDataToRoot(
		ConsolidationStatus & status
	) const;

	///	<summary>
	///		Compute vorticity and divergence on the grid.
	///	</summary>
	virtual void ComputeVorticityDivergence(
		int iDataIndex
	);

	///	<summary>
	///		Interpolate data vertically from Nodes to REdges.
	///	</summary>
	void InterpolateNodeToREdge(
		int iVar,
		int iDataIndex
	);

	///	<summary>
	///		Interpolate data vertically from REdges to Nodes.
	///	</summary>
	void InterpolateREdgeToNode(
		int iVar,
		int iDataIndex
	);

public:
	///	<summary>
	///		Get the bounds on the reference grid.
	///	</summary>
	virtual void GetReferenceGridBounds(
		double & dX0,
		double & dX1,
		double & dY0,
		double & dY1
	) {
		_EXCEPTIONT("Not implemented");
	}

	///	<summary>
	///		Perform interpolation on a node array and send data to root
	///		(generally used for serial output on reference grid)
	///	</summary>
	void ReduceInterpolate(
		const DataVector<double> & dAlpha,
		const DataVector<double> & dBeta,
		const DataVector<int> & iPatch,
		DataType eDataType,
		DataMatrix3D<double> & dInterpData,
		bool fIncludeReferenceState = true,
		bool fConvertToPrimitive = true
	) const;

	///	<summary>
	///		Convert an array of coordinate variables to coordinates on the
	///		reference grid (RLL on the sphere, Cartesian on the plane).
	///	</summary>
	virtual void ConvertReferenceToPatchCoord(
		const DataVector<double> & dXReference,
		const DataVector<double> & dYReference,
		DataVector<double> & dAlpha,
		DataVector<double> & dBeta,
		DataVector<int> & iPatch
	) const;

protected:
	///	<summary>
	///		Add a patch to the grid.
	///	</summary>
	GridPatch * AddPatch(
		GridPatch * pPatch
	);

	///	<summary>
	///		Distribute patches among processors and allocate local patches.
	///	</summary>
	void DistributePatches();

public:
	///	<summary>
	///		Copy data from one data index to another.
	///	</summary>
	void CopyData(
		int ixSource,
		int ixDest,
		DataType eDataType
	);

	///	<summary>
	///		Compute a linear combination of data and store at specified index.
	///	</summary>
	void LinearCombineData(
		const DataVector<double> & dCoeff,
		int ixDest,
		DataType eDataType
	);

	///	<summary>
	///		Set the state to zero.
	///	</summary>
	void ZeroData(
		int ixData,
		DataType eDataType
	);

	///	<summary>
	///		Add the reference state to the specified state data index.
	///	</summary>
	void AddReferenceState(
		int ix
	);

public:
	///	<summary>
	///		Get a reference to the model.
	///	</summary>
	const Model & GetModel() const {
		return m_model;
	}

	///	<summary>
	///		Get the base resolution of the grid in the alpha direction.
	///	</summary>
	int GetABaseResolution() const {
		return m_nABaseResolution;
	}

	///	<summary>
	///		Get the base resolution of the grid in the beta direction.
	///	</summary>
	int GetBBaseResolution() const {
		return m_nBBaseResolution;
	}

	///	<summary>
	///		Get the refinement ratio.
	///	</summary>
	int GetRefinementRatio() const {
		return m_nRefinementRatio;
	}

	///	<summary>
	///		Get the reference length scale.
	///	</summary>
	double GetReferenceLength() const {
		return m_dReferenceLength;
	}

	///	<summary>
	///		Get the number of radial elements.
	///	</summary>
	int GetRElements() const {
		return m_nRElements;
	}

	///	<summary>
	///		Get the altitude of the model cap.
	///	</summary>
	double GetZtop() const {
		return m_dZtop;
	}

	///	<summary>
	///		Get the radial coordinate at the given level.
	///	</summary>
	double GetREtaLevel(int ix) const {
		return m_dREtaLevels[ix];
	}

	///	<summary>
	///		Get the vector of radial element levels.
	///	</summary>
	const DataVector<double> & GetREtaLevels() const {
		return m_dREtaLevels;
	}

	///	<summary>
	///		Get the radial coordinate at the given interface.
	///	</summary>
	double GetREtaInterface(int ix) const {
		return m_dREtaInterfaces[ix];
	}

	///	<summary>
	///		Get the vector of radial element interfaces.
	///	</summary>
	const DataVector<double> & GetREtaInterfaces() const {
		return m_dREtaInterfaces;
	}

	///	<summary>
	///		Get the DataLocation of the specified equation set variable.
	///	</summary>
	DataLocation GetVarLocation(int ix) const {
		return m_vecVarLocation[ix];
	}

	///	<summary>
	///		Get the local variable index associated with this variable.
	///	</summary>
	int GetVarIndex(int ix) const {
		return m_vecVarIndex[ix];
	}

	///	<summary>
	///		Get the number of variables at the specified DataLocation.
	///	</summary>
	int GetVarsAtLocation(DataLocation loc) const {
		return m_vecVarsAtLocation[(int)loc];
	}

	///	<summary>
	///		Get the reference state availability flag.
	///	</summary>
	bool HasReferenceState() const {
		return m_fHasReferenceState;
	}

public:
	///	<summary>
	///		Get the specified cumulative patch 2D node index.
	///	</summary>
	int GetCumulativePatch2DNodeIndex(int ix) const {
		if ((ix < 0) || (ix >= m_vecCumulativePatch2DNodeIndex.size())) {
			_EXCEPTIONT("Invalid patch index");
		}
		return m_vecCumulativePatch2DNodeIndex[ix];
	}

	///	<summary>
	///		Get the specified cumulative patch 3D node index.
	///	</summary>
	int GetCumulativePatch3DNodeIndex(int ix) const {
		if ((ix < 0) || (ix >= m_vecCumulativePatch2DNodeIndex.size())) {
			_EXCEPTIONT("Invalid patch index");
		}
		return m_vecCumulativePatch2DNodeIndex[ix] * GetRElements();
	}

	///	<summary>
	///		Get the patch with the specified index.
	///	</summary>
	GridPatch * GetPatch(int ix) {
		if ((ix < 0) || (ix >= m_vecGridPatches.size())) {
			_EXCEPTIONT("Invalid patch index");
		}
		return m_vecGridPatches[ix];
	}

	///	<summary>
	///		Get the patch with the specified index.
	///	</summary>
	const GridPatch * GetPatch(int ix) const {
		if ((ix < 0) || (ix >= m_vecGridPatches.size())) {
			_EXCEPTIONT("Invalid patch index");
		}
		return m_vecGridPatches[ix];
	}

	///	<summary>
	///		Get the active patch with the specified index.
	///	</summary>
	GridPatch * GetActivePatch(int ix) {
		if ((ix < 0) || (ix >= m_vecActiveGridPatches.size())) {
			_EXCEPTIONT("Invalid active patch index");
		}
		return m_vecActiveGridPatches[ix];
	}

	///	<summary>
	///		Get the active patch with the specified index.
	///	</summary>
	const GridPatch * GetActivePatch(int ix) const {
		if ((ix < 0) || (ix >= m_vecActiveGridPatches.size())) {
			_EXCEPTIONT("Invalid active patch index");
		}
		return m_vecActiveGridPatches[ix];
	}

protected:
	///	<summary>
	///		Initialization flag.
	///	</summary>
	bool m_fInitialized;

	///	<summary>
	///		Reference to the model.
	///	</summary>
	const Model & m_model;

	///	<summary>
	///		Base resolution of the method in the alpha direction.
	///	</summary>
	int m_nABaseResolution;

	///	<summary>
	///		Base resolution of the method in the beta direction.
	///	</summary>
	int m_nBBaseResolution;

	///	<summary>
	///		Refinement ratio.
	///	</summary>
	int m_nRefinementRatio;

	///	<summary>
	///		Reference length scale.
	///	</summary>
	double m_dReferenceLength;

protected:
	///	<summary>
	///		Number of radial elements.
	///	</summary>
	int m_nRElements;

	///	<summary>
	///		Model height cap.
	///	</summary>
	double m_dZtop;

	///	<summary>
	///		REta coordinates of levels along radial axis.
	///	</summary>
	DataVector<double> m_dREtaLevels;

	///	<summary>
	///		REta coordinates of interfaces along radial axis.
	///	</summary>
	DataVector<double> m_dREtaInterfaces;

	///	<summary>
	///		Location of each equation set variable.
	///	</summary>
	std::vector<DataLocation> m_vecVarLocation;

	///	<summary>
	///		Map from equation set variable index to location-dependent
	///		structure.
	///	</summary>
	std::vector<int> m_vecVarIndex;

	///	<summary>
	///		Number of state variables at each location.
	///	</summary>
	std::vector<int> m_vecVarsAtLocation;

	///	<summary>
	///		Flag indicating whether or not a reference state is available.
	///	</summary>
	bool m_fHasReferenceState;

private:
	///	<summary>
	///		Vector of cumulative indices associated with 2D nodal values on
	///		each patch.
	///	</summary>
	std::vector<int> m_vecCumulativePatch2DNodeIndex;

	///	<summary>
	///		Vector of grid patches which are active locally.
	///	</summary>
	GridPatchVector m_vecActiveGridPatches;

	///	<summary>
	///		Vector of grid patches.
	///	</summary>
	GridPatchVector m_vecGridPatches;
};

///////////////////////////////////////////////////////////////////////////////

#endif

