///////////////////////////////////////////////////////////////////////////////
///
///	\file    HorizontalDynamicsFEM.h
///	\author  Paul Ullrich
///	\version June 18, 2013
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

#ifndef _HORIZONTALDYNAMICSFEM_H_
#define _HORIZONTALDYNAMICSFEM_H_

#include "HorizontalDynamics.h"
#include "DataVector.h"
#include "DataMatrix.h"
#include "DataMatrix3D.h"

///////////////////////////////////////////////////////////////////////////////

class GridData3D;
class GridData4D;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Finite-element method (FEM) based atmospheric horizontal dynamics.
///	</summary>
class HorizontalDynamicsFEM : public HorizontalDynamics {

public:
	///	<summary>
	///		Type of HorizontalDynamics.
	///	</summary>
	enum Type {
		SpectralElement,
		DiscontinuousGalerkin
	};

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	HorizontalDynamicsFEM(
		Model & model,
		int nHorizontalOrder,
		Type eHorizontalDynamicsType = SpectralElement,
		bool fUseHyperdiffusion = false
	);

public:
	///	<summary>
	///		Generate the hyperdiffusion matrix of the specified order.
	///	</summary>
	void GenerateHyperdiffusionMatrix();

public:
	///	<summary>
	///		Get the number of halo elements needed by the model.
	///	</summary>
	virtual int GetHaloElements() const {
		return 1;
	}
/*
protected:
	///	<summary>
	///		Apply the DSS procedure.
	///	</summary>
	void ApplyDSS(
		int iDataUpdate
	);
*/
public:
	///	<summary>
	///		Perform one Forward Euler step for the interior terms of the
	///		shallow water equations.
	///	</summary>
	void StepShallowWater(
		int iDataInitial,
		int iDataUpdate,
		double dTime,
		double dDeltaT
	);

	///	<summary>
	///		Perform one Forward Euler step for the element fluxes of the
	///		shallow water equations.
	///	</summary>
	void ElementFluxesShallowWater(
		int iDataInitial,
		int iDataUpdate,
		double dTime,
		double dDeltaT
	);

	///	<summary>
	///		Perform one horizontal Forward Euler step for the interior terms of
	///		the non-hydrostatic primitive equations.
	///	</summary>
	void StepNonhydrostaticPrimitive(
		int iDataInitial,
		int iDataUpdate,
		double dTime,
		double dDeltaT
	);

public:
	///	<summary>
	///		Perform one horizontal Forward Euler step.
	///	</summary>
	virtual void StepExplicit(
		int iDataInitial,
		int iDataUpdate,
		double dTime,
		double dDeltaT
	);

protected:
	///	<summary>
	///		Apply the scalar Laplacian operator.
	///	</summary>
	void ApplyScalarHyperdiffusion(
		int iDataInitial,
		int iDataUpdate,
		double dDeltaT,
		double dNu,
		bool fScaleNuLocally
	);

	///	<summary>
	///		Apply the vector Laplacian operator.
	///	</summary>
	void ApplyVectorHyperdiffusion(
		int iDataInitial,
		int iDataUpdate,
		double dDeltaT,
		double dNuDiff,
		double dNuVort,
		bool fScaleNuLocally
	);

public:
	///	<summary>
	///		Apply hyperdiffusion.
	///	</summar>
	virtual void StepAfterSubCycle(
		int iDataInitial,
		int iDataUpdate,
		double dTime,
		double dDeltaT
	);

private:
	///	<summary>
	///		Type of dynamics to use.
	///	</summary>
	Type m_eHorizontalDynamicsType;

	///	<summary>
	///		Spatial order of accuracy.
	///	</summary>
	int m_nHorizontalOrder;
/*
	///	<summary>
	///		Type of mass matrix employed by this method.
	///	</summary>
	MassMatrixType m_eMassMatrixType;

	///	<summary>
	///		Horizontal time integrator being used.
	///	</summary>
	TimeIntegrator m_eTimeIntegrator;
*/
	///	<summary>
	///		Derivatives of the basis functions at nodal points on the
	///		reference element.
	///	</summary>
	DataMatrix<double> m_dDxBasis1D;

	///	<summary>
	///		Components of the stiffness matrix.
	///	</summary>
	DataMatrix<double> m_dStiffness1D;

	///	<summary>
	///		Derivatives of the flux reconstruction function (used by
	///		discontinuous Galerkin dynamics).
	///	</summary>
	DataVector<double> m_dFluxDeriv1D;

	///	<summary>
	///		Nodal weights in the reference element.
	///	</summary>
	DataVector<double> m_dGLLWeight;

	///	<summary>
	///		Nodal alpha fluxes.
	///	</summary>
	DataMatrix<double> m_dAlphaFlux;

	///	<summary>
	///		Nodal beta fluxes.
	///	</summary>
	DataMatrix<double> m_dBetaFlux;

	///	<summary>
	///		Nodal pointwise pressures.
	///	</summary>
	DataMatrix<double> m_dPressure;

	///	<summary>
	///		Map from model interfaces to model levels.
	///	</summary>
	DataMatrix<double> m_matInterfaceToLevel;

	///	<summary>
	///		Map from model levels to interfaces.
	///	</summary>
	DataMatrix<double> m_matLevelToInterface;

private:
	///	<summary>
	///		Flag indicating whether or not hyperdiffusion should be used.
	///	</summary>
	bool m_fUseHyperdiffusion;

	///	<summary>
	///		Nodal pointwise gradient.
	///	</summary>
	DataMatrix3D<double> m_dGradient;

	///	<summary>
	///		Scalar hyperviscosity coefficient (at 1 degree resolution).
	///	</summary>
	double m_dNuScalar;

	///	<summary>
	///		Divergent hyperviscosity coefficient (at 1 degree resolution).
	///	</summary>
	double m_dNuDiv;

	///	<summary>
	///		Vortical hyperviscosity coefficient (at 1 degree resolution).
	///	</summary>
	double m_dNuVort;

};

///////////////////////////////////////////////////////////////////////////////

#endif

