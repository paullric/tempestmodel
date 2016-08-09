///////////////////////////////////////////////////////////////////////////////
///
///	\file    DataType.h
///	\author  Paul Ullrich
///	\version March 24, 2013
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

#ifndef _DATATYPE_H_
#define _DATATYPE_H_

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		DataTypes that can be exchanged between processors.
///	</summary>
enum DataType {
	DataType_All = (-1),
	DataType_Default = (0),
	DataType_State = DataType_Default,
	DataType_RefState,
	DataType_Tracers,
	DataType_Auxiliary2D,
	DataType_Auxiliary3D,
	DataType_Jacobian,
	DataType_ElementArea,
	DataType_Topography,
	DataType_TopographyDeriv,
	DataType_Longitude,
	DataType_Latitude,
	DataType_Z,
	DataType_Pressure,
	DataType_SurfacePressure,
	DataType_KineticEnergy,
	DataType_Vorticity,
	DataType_Divergence,
	DataType_Temperature,
	DataType_RayleighStrength,
	DataType_Richardson,
	DataType_None
};

///////////////////////////////////////////////////////////////////////////////

#endif
