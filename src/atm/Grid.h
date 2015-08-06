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
#include "MathHelper.h"

#include "mpi.h"

#include <string>
#include <vector>
#include <map>

///////////////////////////////////////////////////////////////////////////////

class Time;
class Model;
class TestCase;
class ConsolidationStatus;
class VerticalStretchFunction;

class NcFile;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Atmospheric model grid data.  Container for GridPatch objects.
///	</summary>
class Grid {

friend class Model;

public:
	///	<summary>
	///		Type of vertical staggering to be used.
	///	</summary>
	enum VerticalStaggering {
		VerticalStaggering_Levels,
		VerticalStaggering_Interfaces,
		VerticalStaggering_CharneyPhillips,
		VerticalStaggering_Lorenz
	};

public:
	///	<summary>
	///		Boundary condition type applied in each direction.
	///	</summary>
	enum BoundaryCondition {
		BoundaryCondition_Default = 0,
		BoundaryCondition_Periodic = BoundaryCondition_Default,
		BoundaryCondition_NoFlux = 1
	};

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	Grid(
		Model & model,
		int nABaseResolution,
		int nBBaseResolution,
		int nRefinementRatio,
		int nRElements,
		VerticalStaggering eVerticalStaggering
			= VerticalStaggering_CharneyPhillips
	);

	///	<summary>
	///		Virtual destructor.
	///	</summary>
	virtual ~Grid();

public:
	///	<summary>
	///		Get the boundary condition in the specified direction.
	///	</summary>
	BoundaryCondition GetBoundaryCondition(
		Direction eDir
	) {
		return m_eBoundaryCondition[static_cast<int>(eDir)];
	}

	///	<summary>
	///		Set the boundary condition in the specified direction.
	///	</summary>
	void SetBoundaryCondition(
		Direction eDir,
		BoundaryCondition eBoundaryCondition
	);

public:
	///	<summary>
	///		Set the flag indicating blocking of parallel exchanges.
	///	</summary>
	void SetBlockParallelExchange(bool fBlockParallelExchange) {
		m_fBlockParallelExchange = fBlockParallelExchange;
	}

public:
	///	<summary>
	///		Set the vertical stretching function.  Grid retains ownership
	///		of the pointer after assignment.
	///	</summary>
	void SetVerticalStretchFunction(
		VerticalStretchFunction * pVerticalStretchF
	);

	///	<summary>
	///		Evaluate the vertical stretching function.
	///	</summary>
	void EvaluateVerticalStretchF(
		double dREta,
		double & dREtaStretch,
		double & dDxREtaStretch
	);

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
	///		Initialize topography height/derivatives, state and tracer data
	///		from a TestCase.
	///	</summary>
	virtual void EvaluateTestCase(
		const TestCase & test,
		const Time & time,
		int iDataIndex = 0
	);

	///	<summary>
	///		Initialize state and tracer data from a TestCase.
	///	</summary>
	void EvaluateTestCase_StateOnly(
		const TestCase & test,
		const Time & time,
		int iDataIndex = 0
	);

	///	<summary>
	///		Evaluate relevant geometric terms.
	///	</summary>
	void EvaluateGeometricTerms();

public:
	///	<summary>
	///		Initialize state and tracer data from a TestCase.
	///	</summary>
	void ApplyBoundaryConditions(
		int iDataIndex = 0,
		DataType eDataType = DataType_State
	);

	///	<summary>
	///		Perform post-processing of variables on the grid after each
	///		TimeStep substage.
	///	</summary>
	virtual void PostProcessSubstage(
		int iDataUpdate,
		DataType eDataType = DataType_State
	) {
		_EXCEPTIONT("Unimplemented");
	}

public:
	///	<summary>
	///		Perform checksum calculation on all state variables.
	///	</summary>
	void Checksum(
		DataType eDataType,
		DataArray1D<double> & dChecksums,
		int iDataIndex = 0,
		ChecksumType eChecksumType = ChecksumType_Sum
	) const;

	///	<summary>
	///		Compute total energy on the grid.
	///	</summary>
	double ComputeTotalEnergy(
		int iDataIndex
	) const;

	///	<summary>
	///		Compute total energy on the grid.
	///	</summary>
	double ComputeTotalPotentialEnstrophy(
		int iDataIndex
	);

public:
	///	<summary>
	///		Exchange data between processors.
	///	</summary>
	void Exchange(
		DataType eDataType,
		int iDataIndex
	);

	///	<summary>
	///		Exchange connectivity buffers between processors.
	///	</summary>
	void ExchangeBuffers();

	///	<summary>
	///		Exchange connectivity buffers between processors.
	///	</summary>
	void ExchangeBuffersAndUnpack(
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
	///		Get the longest perimeter over all patches.
	///	</summary>
	int GetLongestActivePatchPerimeter() const;

	///	<summary>
	///		Get the maximum number of nodes in 2D over all patches.
	///	</summary>
	int GetMaxNodeCount2D() const;

	///	<summary>
	///		Get the total count of nodes in 2D (ignoring the number of vertical
	///		degrees of freedom).
	///	</summary>
	int GetTotalNodeCount2D() const;

	///	<summary>
	///		Get the maximum number of nodes over all patches.
	///	</summary>
	int GetMaxNodeCount(
		DataLocation loc = DataLocation_Node
	) const;

	///	<summary>
	///		Get the total count of nodes over the global grid, returned to
	///		the root node.
	///	</summary>
	int GetTotalNodeCount(
		DataLocation loc = DataLocation_Node
	) const;

	///	<summary>
	///		Get the maximum number of degrees of freedom required for a
	///		consolidate operation, returned to the root node.
	///	</summary>
	int GetMaxDegreesOfFreedom() const;

public:
	///	<summary>
	///		Receive a data object on the root process.
	///	</summary>
	void ConsolidateDataAtRoot(
		ConsolidationStatus & status,
		DataArray1D<double> & dataRecvBuffer,
		int & nRecvCount,
		int & ixRecvPatch,
		DataType & eRecvDataType,
		DataLocation & eRecvDataLocation
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
	///		Compute temperature on the grid.
	///	</summary>
	virtual void ComputeTemperature(
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
	///	<param name="eDataLocation">
	///		DataLocation_Node  = Interpolate all variables on nodes
	///		DataLocation_REdge = Interpolate all variables on redges
	///	</param>
	void ReduceInterpolate(
		const DataArray1D<double> & dAlpha,
		const DataArray1D<double> & dBeta,
		const DataArray1D<int> & iPatch,
		DataType eDataType,
		DataLocation eDataLocation,
		bool fInterpAllVariables,
		DataArray3D<double> & dInterpData,
		bool fIncludeReferenceState = true,
		bool fConvertToPrimitive = true
	) const;

	///	<summary>
	///		Convert an array of coordinate variables to coordinates on the
	///		reference grid (RLL on the sphere, Cartesian on the plane).
	///	</summary>
	virtual void ConvertReferenceToPatchCoord(
		const DataArray1D<double> & dXReference,
		const DataArray1D<double> & dYReference,
		DataArray1D<double> & dAlpha,
		DataArray1D<double> & dBeta,
		DataArray1D<int> & iPatch
	) const;

	///	<summary>
	///		Get the patch and coordinate index for the specified node.
	///	</summary>
	virtual void GetPatchFromCoordinateIndex(
		int iRefinementLevel,
		const DataArray1D<int> & vecIxA,
		const DataArray1D<int> & vecIxB,
		const DataArray1D<int> & vecPanel,
		DataArray1D<int> & vecPatchIndex,
		int nVectorLength = (-1)
	) {
		_EXCEPTIONT("Not implemented.");
	}

	///	<summary>
	///		Get the relation between coordinate vectors across panel
	///		boundaries.
	///	</summary>
	virtual void GetOpposingDirection(
		int ixPanelSrc,
		int ixPanelDest,
		Direction dir,
		Direction & dirOpposing,
		bool & fSwitchParallel,
		bool & fSwitchPerpendicular
	) const {
		_EXCEPTIONT("Not implemented.");
	}

public:
	///	<summary>
	///		Add the default set of patches.
	///	</summary>
	virtual void AddDefaultPatches() {
		_EXCEPTIONT("Not implemented");
	}

	///	<summary>
	///		Add a patch to the grid with the specified index and PatchBox.
	///	</summary>
	virtual GridPatch * AddPatch(
		int ixPatch,
		const PatchBox & box
	) {
		_EXCEPTIONT("Not implemented");
	}

	///	<summary>
	///		Add a patch to the grid.
	///	</summary>
	GridPatch * AddPatch(
		GridPatch * pPatch
	);

	///	<summary>
	///		Write a grid to a NetCDF file.
	///	</summary>
	void ToFile(
		NcFile & ncfile
	);

	///	<summary>
	///		Load a grid from a file.
	///	</summary>
	void FromFile(
		const std::string & strGridFile
	);

protected:
	///	<summary>
	///		Distribute patches among processors and allocate local patches.
	///	</summary>
	void DistributePatches();

	///	<summary>
	///		Initialize connectivity between patches.
	///	</summary>
	void InitializeConnectivity();

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
		const DataArray1D<double> & dCoeff,
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
	///		Get the grid stamp.
	///	</summary>
	int GetGridStamp() const {
		return m_iGridStamp;
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
	///		Get the base resolution at the specified refinement level.
	///	</summary>
	inline int GetABaseResolution(int iRefineLevel) {
		return (m_nABaseResolution
			* IntPow(
				m_nRefinementRatio,
				static_cast<unsigned int>(iRefineLevel)));
	}

	///	<summary>
	///		Get the base resolution at the specified refinement level.
	///	</summary>
	inline int GetBBaseResolution(int iRefineLevel) {
		return (m_nBBaseResolution
			* IntPow(
				m_nRefinementRatio,
				static_cast<unsigned int>(iRefineLevel)));
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
		return (m_vecCumulativePatch2DNodeIndex[ix]
			* m_nDegreesOfFreedomPerColumn);
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

public:
	///	<summary>
	///		Get the reference length scale.
	///	</summary>
	double GetReferenceLength() const {
		return m_dReferenceLength;
	}

	///	<summary>
	///		Set the reference length scale.
	///	</summary>
	void SetReferenceLength(double dReferenceLength) {
		m_dReferenceLength = dReferenceLength;
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
	const DataArray1D<double> & GetREtaLevels() const {
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
	const DataArray1D<double> & GetREtaInterfaces() const {
		return m_dREtaInterfaces;
	}

	///	<summary>
	///		Get the vector of normalized areas for radial element levels.
	///	</summary>
	const DataArray1D<double> & GetREtaLevelsNormArea() const {
		return m_dREtaLevelsNormArea;
	}

	///	<summary>
	///		Get the vector of normalized areas for radial element interfaces.
	///	</summary>
	const DataArray1D<double> & GetREtaInterfacesNormArea() const {
		return m_dREtaInterfacesNormArea;
	}

	///	<summary>
	///		Get the radial coordinate at the given level.
	///	</summary>
	double GetREtaStretchLevel(int ix) const {
		return m_dREtaStretchLevels[ix];
	}

	///	<summary>
	///		Get the vector of radial element levels.
	///	</summary>
	const DataArray1D<double> & GetREtaStretchLevels() const {
		return m_dREtaStretchLevels;
	}

	///	<summary>
	///		Get the radial coordinate at the given interface.
	///	</summary>
	double GetREtaStretchInterface(int ix) const {
		return m_dREtaStretchInterfaces[ix];
	}

	///	<summary>
	///		Get the vector of radial element interfaces.
	///	</summary>
	const DataArray1D<double> & GetREtaStretchInterfaces() const {
		return m_dREtaStretchInterfaces;
	}

	///	<summary>
	///		Get the type of vertical staggering.
	///	</summary>
	VerticalStaggering GetVerticalStaggering() const {
		return m_eVerticalStaggering;
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
	///		Get the number of degrees of freedom per column.
	///	</summary>
	int GetDegreesOfFreedomPerColumn() const {
		return m_nDegreesOfFreedomPerColumn;
	}

	///	<summary>
	///		Get the reference state availability flag.
	///	</summary>
	bool HasReferenceState() const {
		return m_fHasReferenceState;
	}

	///	<summary>
	///		Get the active Rayleigh friction flag.
	///	</summary>
	bool HasRayleighFriction() const {
		return m_fHasRayleighFriction;
	}

protected:
	///	<summary>
	///		Initialization flag.
	///	</summary>
	bool m_fInitialized;

	///	<summary>
	///		Reference to the model.
	///	</summary>
	Model & m_model;

	///	<summary>
	///		Block exchange operations between processors.
	///	</summary>
	bool m_fBlockParallelExchange;

	///	<summary>
	///		Pointer to the vertical stretching function.
	///	</summary>
	VerticalStretchFunction * m_pVerticalStretchF;

protected:
	///	<summary>
	///		Vector of cumulative indices associated with 2D nodal values on
	///		each patch.
	///	</summary>
	std::vector<int> m_vecCumulativePatch2DNodeIndex;

	///	<summary>
	///		Vector of grid patches.
	///	</summary>
	GridPatchVector m_vecGridPatches;

	///	<summary>
	///		Vector of grid patches which are active locally.
	///	</summary>
	GridPatchVector m_vecActiveGridPatches;

	///	<summary>
	///		Type of vertical stretching being applied.
	///	</summary>
	VerticalStaggering m_eVerticalStaggering;

protected:
	///	<summary>
	///		DataContainer for Grid data.
	///	</summary>
	DataContainer m_dcGridData;

protected:
	///	<summary>
	///		Boundary condition in each coordinate direction.
	///	</summary>
	DataArray1D<BoundaryCondition> m_eBoundaryCondition;

	///	<summary>
	///		Grid stamp.  This value is incremented whenever the grid changes.
	///	</summary>
	int m_iGridStamp;

	///	<summary>
	///		Base resolution in the alpha direction.
	///	</summary>
	int m_nABaseResolution;

	///	<summary>
	///		Base resolution in the beta direction.
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
	DataArray1D<double> m_dREtaLevels;

	///	<summary>
	///		REta coordinates of interfaces along radial axis.
	///	</summary>
	DataArray1D<double> m_dREtaInterfaces;

	///	<summary>
	///		Normalized area of REta elements along radial axis.
	///	</summary>
	DataArray1D<double> m_dREtaLevelsNormArea;

	///	<summary>
	///		Normalized area of REta interfaces along radial axis.
	///	</summary>
	DataArray1D<double> m_dREtaInterfacesNormArea;

	///	<summary>
	///		Stretched REta coordinates of levels along radial axis.
	///	</summary>
	DataArray1D<double> m_dREtaStretchLevels;

	///	<summary>
	///		Stretched REta coordinates of interfaces along radial axis.
	///	</summary>
	DataArray1D<double> m_dREtaStretchInterfaces;

	///	<summary>
	///		Location of each equation set variable.
	///	</summary>
	DataArray1D<DataLocation> m_vecVarLocation;

	///	<summary>
	///		Map from equation set variable index to location-dependent
	///		structure.
	///	</summary>
	DataArray1D<int> m_vecVarIndex;

	///	<summary>
	///		Number of state variables at each location.
	///	</summary>
	DataArray1D<int> m_vecVarsAtLocation;

	///	<summary>
	///		Number of degrees of freedom per column.
	///	</summary>
	int m_nDegreesOfFreedomPerColumn;

	///	<summary>
	///		Flag indicating whether or not a reference state is available.
	///	</summary>
	bool m_fHasReferenceState;

	///	<summary>
	///		Flag indicating whether or not Rayleigh friction is used.
	///	</summary>
	bool m_fHasRayleighFriction;
};

///////////////////////////////////////////////////////////////////////////////

#endif

