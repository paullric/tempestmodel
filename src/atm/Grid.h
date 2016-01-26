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
#include "Connectivity.h"
#include "DataStruct.h"

#ifdef USE_MPI
#include "mpi.h"
#endif

#include <string>
#include <vector>
#include <map>

///////////////////////////////////////////////////////////////////////////////

class Time;
class Model;
class TestCase;
class ConsolidationStatus;
class GridSpacing;
class VerticalStretchFunction;

#ifndef NO_NETCDF
class NcFile;
#else
typedef int NcFile;
#endif

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
	typedef int VerticalStaggering;
	static const int VerticalStaggering_Levels = 0;
	static const int VerticalStaggering_Interfaces = 1;
	static const int VerticalStaggering_CharneyPhillips = 2;
	static const int VerticalStaggering_Lorenz = 3;

public:
	///	<summary>
	///		Boundary condition type applied in each direction.
	///	</summary>
	typedef int BoundaryCondition;
	static const int BoundaryCondition_Default = 0;
	static const int BoundaryCondition_Periodic = 0;
	static const int BoundaryCondition_NoFlux = 1;
        static const int BoundaryCondition_NoSlip = 2;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	Grid(
		Model & model
	);

	///	<summary>
	///		Virtual destructor.
	///	</summary>
	virtual ~Grid();

	///	<summary>
	///		Define the parameters for the Grid.
	///	</summary>
	virtual void DefineParameters();

	///	<summary>
	///		Set the parameters for the Grid
	///	</summary>
	virtual void SetParameters(
		int nRElements,
		int nMaxPatchCount,
		int nABaseResolution,
		int nBBaseResolution,
		int nRefinementRatio,
		VerticalStaggering eVerticalStaggering
			= VerticalStaggering_CharneyPhillips
	);

	///	<summary>
	///		Allocate GridData.
	///	</summary>
	void InitializeDataLocal();

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

		// Initialize the vertical coordinate
		InitializeVerticalCoordinate();

		// Set the grid as initialized
		m_fInitialized = true;
	}

	///	<summary>
	///		Initialize the vertical coordinate.
	///	</summary>
	virtual void InitializeVerticalCoordinate();

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
		return m_nInitializedPatchBoxes;
	}

	///	<summary>
	///		Get the total number of active patches on this processor.
	///	</summary>
	int GetActivePatchCount() const {
		return m_vecActiveGridPatches.size();
	}

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

public:
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
	///		Build the default patch layout.
	///	</summary>
	virtual void ApplyDefaultPatchLayout(
		int nPatchCount
	) {
		_EXCEPTIONT("Not implemented");
	}

	///	<summary>
	///		Return a pointer to a new GridPatch.
	///	</summary>
	virtual GridPatch * NewPatch(
		int ixPatch
	) {
		_EXCEPTIONT("Not implemented");
	}

	///	<summary>
	///		Create a new empty GridPatch and activate it.
	///	</summary>
	GridPatch * ActivateEmptyPatch(
		int ixPatch
	);

	///	<summary>
	///		Deactivate the given patch.
	///	</summary>
	void DeactivatePatch(
		int ixPatch
	);

protected:
	///	<summary>
	///		Add a patch to the grid.
	///	</summary>
	GridPatch * AddPatch(
		GridPatch * pPatch
	);

public:
	///	<summary>
	///		Get a reference to the exchange buffer registry.
	///	</summary>
	const ExchangeBufferRegistry & GetExchangeBufferRegistry() const {
		return m_aExchangeBufferRegistry;
	}

	///	<summary>
	///		Get a reference to the exchange buffer registry.
	///	</summary>
	ExchangeBufferRegistry & GetExchangeBufferRegistry() {
		return m_aExchangeBufferRegistry;
	}

public:
	///	<summary>
	///		Distribute patches among processors and allocate local patches.
	///	</summary>
	void DistributePatches();

protected:
	///	<summary>
	///		Register an ExchangeBuffer.
	///	</summary>
	void RegisterExchangeBuffer(
		int ixSourcePatch,
		int ixTargetPatch,
		Direction dir
	);

	///	<summary>
	///		Register an ExchangeBuffer.
	///	</summary>
	void RegisterExchangeBuffer(
		int ixSourcePatch,
		int ixTargetPatch,
		Direction dir,
		int ixFirst,
		int ixSecond
	);

	///	<summary>
	///		Build exchange buffer information for the specified GridPatch.
	///	</summary>
	void InitializeExchangeBuffersFromPatch(
		int ixSourcePatch
	);

public:
	///	<summary>
	///		Build active exchange buffer information using the grid layout.
	///	</summary>
	void InitializeExchangeBuffersFromActivePatches();

	///	<summary>
	///		Build all exchange buffer information using the grid layout.
	///	</summary>
	void InitializeAllExchangeBuffers();

	///	<summary>
	///		Initialize connectivity between patches.
	///	</summary>
	void InitializeConnectivity(
		bool fAllocate = true
	);

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
	///		Get the DataContainer storing Grid parameters.
	///	</summary>
	const DataContainer & GetDataContainerParameters() const {
		return m_dcGridParameters;
	}

	///	<summary>
	///		Get the DataContainer storing Grid data.
	///	</summary>
	const DataContainer & GetDataContainerPatchData() const {
		return m_dcGridPatchData;
	}

public:
	///	<summary>
	///		Get the DataContainer storing Grid parameters.
	///	</summary>
	DataContainer & GetDataContainerParameters() {
		return m_dcGridParameters;
	}

	///	<summary>
	///		Get the DataContainer storing Grid data.
	///	</summary>
	DataContainer & GetDataContainerPatchData() {
		return m_dcGridPatchData;
	}

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

#ifdef USE_MPI
	///	<summary>
	///		Get the node that contains the specified GridPatch.
	///	</summary>
	int GetPatchProcessor(int ixPatch) const {
		if ((ixPatch < 0) || (ixPatch >= m_vecPatchProcessor.size())) {
			_EXCEPTIONT("Invalid active patch index");
		}
		return m_vecPatchProcessor[ixPatch];
	}
#endif

public:
	///	<summary>
	///		Get the PatchBox with the specified index.
	///	</summary>
	const PatchBox & GetPatchBox(int ix) const {
		if ((ix < 0) || (ix > m_aPatchBoxes.GetRows())) {
			_EXCEPTIONT("Invalid PatchBox index");
		}
		return m_aPatchBoxes[ix];
	}

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
	///		Get the type of vertical staggering.
	///	</summary>
	VerticalStaggering GetVerticalStaggering() const {
		return m_eVerticalStaggering;
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
	///		Get the order of accuracy of the method.
	///	</summary>
	int GetHorizontalOrder() const {
		return m_nHorizontalOrder;
	}

	///	<summary>
	///		Get the order of accuracy in the vertical.
	///	</summary>
	int GetVerticalOrder() const {
		return m_nVerticalOrder;
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
/*
	///	<summary>
	///		Vector of cumulative indices associated with 2D nodal values on
	///		each patch.
	///	</summary>
	std::vector<int> m_vecCumulativePatch2DNodeIndex;
*/
	///	<summary>
	///		Vector of grid patches which are active locally.
	///	</summary>
	GridPatchVector m_vecActiveGridPatches;

	///	<summary>
	///		Vector of grid patch ids which are active locally.
	///	</summary>
	std::vector<int> m_vecActiveGridPatchIndices;

#ifdef USE_MPI
	///	<summary>
	///		Vector of processors that contain the specified GridPatch.
	///	</summary>
	std::vector<int> m_vecPatchProcessor;
#endif

	///	<summary>
	///		Exchange buffer registry.
	///	</summary>
	ExchangeBufferRegistry m_aExchangeBufferRegistry;

protected:
	///	<summary>
	///		DataContainer for Grid parameters.
	///	</summary>
	DataContainer m_dcGridParameters;

protected:
	///	<summary>
	///		Maximum number of PatchBoxes.
	///	</summary>
	DataStruct<int> m_nMaxPatchCount;

	///	<summary>
	///		Number of radial elements.
	///	</summary>
	DataStruct<int> m_nRElements;

	///	<summary>
	///		Base resolution in the alpha direction.
	///	</summary>
	DataStruct<int> m_nABaseResolution;

	///	<summary>
	///		Base resolution in the beta direction.
	///	</summary>
	DataStruct<int> m_nBBaseResolution;

	///	<summary>
	///		Refinement ratio.
	///	</summary>
	DataStruct<int> m_nRefinementRatio;

	///	<summary>
	///		Type of vertical stretching being applied.
	///	</summary>
	DataStruct<VerticalStaggering> m_eVerticalStaggering;

	///	<summary>
	///		Order of accuracy in the horizontal.
	///	</summary>
	DataStruct<int> m_nHorizontalOrder;

	///	<summary>
	///		Order of accuracy in the vertical.
	///	</summary>
	DataStruct<int> m_nVerticalOrder;


protected:
	///	<summary>
	///		DataContainer for Grid data.
	///	</summary>
	DataContainer m_dcGridPatchData;

protected:
	///	<summary>
	///		Number of initialized PatchBoxes.
	///	</summary>
	DataStruct<int> m_nInitializedPatchBoxes;

	///	<summary>
	///		Array of PatchBoxes.
	///	</summary>
	DataArray1D<PatchBox> m_aPatchBoxes;

	///	<summary>
	///		Boundary condition in each coordinate direction.
	///	</summary>
	DataArray1D<BoundaryCondition> m_eBoundaryCondition;

	///	<summary>
	///		Grid stamp.  This value is incremented whenever the grid changes.
	///	</summary>
	DataStruct<int> m_iGridStamp;

	///	<summary>
	///		Reference length scale.
	///	</summary>
	DataStruct<double> m_dReferenceLength;

protected:
	///	<summary>
	///		Model height cap.
	///	</summary>
	DataStruct<double> m_dZtop;

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
	///		Flag indicating whether or not a reference state is available.
	///	</summary>
	DataStruct<bool> m_fHasReferenceState;

	///	<summary>
	///		Flag indicating whether or not Rayleigh friction is used.
	///	</summary>
	DataStruct<bool> m_fHasRayleighFriction;
};

///////////////////////////////////////////////////////////////////////////////

#endif

