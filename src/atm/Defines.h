///////////////////////////////////////////////////////////////////////////////
///
///	\file    Defines.h
///	\author  Paul Ullrich
///	\version August 5, 2013
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

#ifndef _DEFINES_H_
#define _DEFINES_H_

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Method to solve the vertical system.
///	</summary>
//#define USE_JFNK_PETSC
//#define USE_JFNK_GMRES
//#define USE_DIRECTSOLVE_APPROXJ
#define USE_DIRECTSOLVE

///	<summary>
///		Jacobian storage form
///	</summary>
//#define USE_JACOBIAN_DEBUG
#define USE_JACOBIAN_GENERAL
//#define USE_JACOBIAN_DIAGONAL

///	<summary>
///		Thermodynamic closure to use.
///	</summary>
//#define FORMULATION_PRESSURE
//#define FORMULATION_RHOTHETA_PI
//#define FORMULATION_RHOTHETA_P
#define FORMULATION_THETA
//#define FORMULATION_THETA_FLUX

///	<summary>
///		Compute residuals in TimestepScheme classes.
///	</summary>
//#define RESIDUAL_DIFFUSION

///	<summary>
///		Hardcoded horizontal order.
///	</summary>
//#define FIXED_HORIZONTAL_ORDER 4

///	<summary>
///		Hardcoded number of levels.
///	</summary>
//#define FIXED_RELEMENTS 20

///	<summary>
///		Prognostic contravariant velocities.
///	</summary>
//#define PROGNOSTIC_CONTRAVARIANT_MOMENTA

///	<summary>
///		Explicit vertical advection of vertical velocity variable.
///	</summary>
//#define EXPLICIT_VERTICAL_VELOCITY_ADVECTION

///	<summary>
///		Disable uniform diffusion checks.
///	</summary>
//#define DISABLE_UNIFORM_DIFFUSION_CHECKS

///	<summary>
///		When to apply Rayleigh damping.
///	</summary>
//#define APPLY_RAYLEIGH_WITH_HYPERVIS
#define APPLY_RAYLEIGH_WITH_VERTICALDYN

///	<summary>
///		Apply positive definite filter to tracers.
///	</summary>
#define POSITIVE_DEFINITE_FILTER_TRACERS

///////////////////////////////////////////////////////////////////////////////

#endif
