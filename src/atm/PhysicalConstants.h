///////////////////////////////////////////////////////////////////////////////
///
///	\file    PhysicalConstants.h
///	\author  Paul Ullrich
///	\version February 24, 2013
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

#ifndef _PHYSICALCONSTANTS_H_
#define _PHYSICALCONSTANTS_H_

#include "Exception.h"

#include <cmath>

#ifdef TEMPEST_NETCDF
#include <netcdfcpp.h>
#endif

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Physical constants which are relevant to the simulation.
///	</summary>
class PhysicalConstants {

private:
	///	<summary>
	///		Radius of the Earth.
	///	</summary>
	double m_dEarthRadius;

	///	<summary>
	///		Constant of surface gravity (m/s^2)
	///	</summary>
	double m_dG;

	///	<summary>
	///		Rotation rate of the Earth (1/s)
	///	</summary>
	double m_dOmega;

	///	<summary>
	///		Inclination of the grid
	///	</summary>
	double m_dAlpha;

	///	<summary>
	///		Ideal gas constant of dry air (J/kg/K)
	///	</summary>
	double m_dR;

	///	<summary>
	///		Specific heat capacity of dry air at constant pressure (J/kg/K)
	///	</summary>
	double m_dCp;

	///	<summary>
	///		Reference temperature (K)
	///	</summary>
	double m_dT0;

	///	<summary>
	///		Reference pressure (Pa)
	///	</summary>
	double m_dP0;

	///	<summary>
	///		Reference density of water (kg/m^3)
	///	</summary>
	double m_dRhoWater;

	///	<summary>
	///		Ideal gas constant for water vapor (J/kg/K)
	///	</summary>
	double m_dRvap;

	///	<summary>
	///		Ratio of molar mass of dry air to molar mass of water vapor
	///		(dimensionless)
	///	</summary>
	double m_dMvap;

	///	<summary>
	///		Latent heat of vaporization of water (J / kg)
	///	</summary>
	double m_dLvap;

private:
	///	<summary>
	///		Ratio of the ideal gas constant and specific heat capacity at
	///		constant pressure (dimensionless)
	///	</summary>
	double m_dKappa;

	///	<summary>
	///		Polytropic gas constant (dimensionless)
	///	</summary>
	double m_dGamma;

	///	<summary>
	///		Pressure scaling, used to calculate the pressure from the
	///		potential temperature density (Pa)
	///	</summary>
	double m_dPressureScaling;

public:
	///	<summary>
	///		Construct a new PhysicalConstants object with default values
	///		for all physical constants.
	///	</summary>
	PhysicalConstants() :
		m_dEarthRadius(6.37122e6),
		m_dG(9.80616),
		m_dOmega(7.29212e-5),
		m_dAlpha(0.0),
		m_dR(287.0),
		m_dCp(1004.5),
		m_dT0(300.0),
		m_dP0(100000.0),
		m_dRhoWater(1000.0),
		m_dRvap(461.5),
		m_dMvap(0.608),
		m_dLvap(2.5e6)
	{
		RecalculateKappa();
		RecalculateGamma();
		RecalculatePressureScaling();
	}

public:
	///	<summary>
	///		Return the radius of the Earth (m)
	///	</summary>
	inline double GetEarthRadius() const {
		return m_dEarthRadius;
	}

	///	<summary>
	///		Return the constant of surface gravity (m/s^2)
	///	</summary>
	inline double GetG() const {
		return m_dG;
	}

	///	<summary>
	///		Return the rotation rate of the Earth (1/s)
	///	</summary>
	inline double GetOmega() const {
		return m_dOmega;
	}

	///	<summary>
	///		Return the grid inclination (radians)
	///	</summary>
	inline double GetAlpha() const {
		return m_dAlpha;
	}

	///	<summary>
	///		Return the ideal gas constant for dry air (J/kg/K)
	///	</summary>
	inline double GetR() const {
		return m_dR;
	}

	///	<summary>
	///		Return the specific heat capacity at constant pressure for
	///		dry air (J/kg/K)
	///	</summary>
	inline double GetCp() const {
		return m_dCp;
	}

	///	<summary>
	///		Return the specific heat capacity at constant volume for
	///		dry air (J/kg/K)
	///	</summary>
	inline double GetCv() const {
		return (m_dCp - m_dR);
	}

	///	<summary>
	///		Return the reference temperature (K)
	///	</summary>
	inline double GetT0() const {
		return m_dT0;
	}

	///	<summary>
	///		Return the reference pressure (Pa)
	///	</summary>
	inline double GetP0() const {
		return m_dP0;
	}

	///	<summary>
	///		Return the density of water (kg/m^3)
	///	</summary>
	inline double GetRhoWater() const {
		return m_dRhoWater;
	}

	///	<summary>
	///		Return the ideal gas constant for water vapor (J/kg/K)
	///	</summary>
	inline double GetRvap() const {
		return m_dRvap;
	}

	///	<summary>
	///		Return the ratio of molar mass of dry air to molar mass of
	///		water vapor (dimensionless)
	///	</summary>
	inline double GetMvap() const {
		return m_dMvap;
	}

	///	<summary>
	///		Return the latent heat of vaporization of water vapor (J/kg)
	///	</summary>
	inline double GetLvap() const {
		return m_dLvap;
	}

public:
	///	<summary>
	///		Return the ratio of ideal gas constant to specific heat capacity
	///		at constant pressure for dry air
	///	</summary>
	inline double GetKappa() const {
		return m_dKappa;
	}

	///	<summary>
	///		Return the polytropic gas constant for dry air (dimensionless)
	///	</summary>
	inline double GetGamma() const {
		return m_dGamma;
	}

	///	<summary>
	///		Return the potential temperature pressure scaling coefficient (Pa)
	///	</summary>
	inline double GetPressureScaling() const {
		return m_dPressureScaling;
	}

public:
	///	<summary>
	///		Set the radius of the Earth (m)
	///	</summary>
	inline void SetEarthRadius(double dEarthRadius) {
		m_dEarthRadius = dEarthRadius;
	}

	///	<summary>
	///		Set the constant of surface gravity (m/s^2)
	///	</summary>
	inline void SetG(double dG) {
		m_dG = dG;
	}

	///	<summary>
	///		Set the rotation rate of the Earth (1/s)
	///	</summary>
	inline void SetOmega(double dOmega) {
		m_dOmega = dOmega;
	}

	///	<summary>
	///		Set the grid inclination (radians)
	///	</summary>
	inline void SetAlpha(double dAlpha) {
		m_dAlpha = dAlpha;
	}

	///	<summary>
	///		Set the ideal gas constant of dry air (J/kg/K).  Kappa, Gamma
	///		and the pressure scaling will be recalculated.
	///	</summary>
	inline void SetR(double dR) {
		m_dR = dR;

		RecalculateKappa();
		RecalculateGamma();
		RecalculatePressureScaling();
	}

	///	<summary>
	///		Set the specific heat capacity of dry air at constant pressure
	///		(J/kg/K).  Kappa, Gamma and the pressure scaling will be
	///		recalculated.
	///	</summary>
	inline void SetCp(double dCp) {
		m_dCp = dCp;

		RecalculateKappa();
		RecalculateGamma();
		RecalculatePressureScaling();
	}

	///	<summary>
	///		Set the reference temperature (K)
	///	</summary>
	inline void SetT0(double dT0) {
		m_dT0 = dT0;
	}

	///	<summary>
	///		Set the reference pressure (Pa)
	///	</summary>
	inline void SetP0(double dP0) {
		m_dP0 = dP0;

		RecalculatePressureScaling();
	}

	///	<summary>
	///		Set the density of water (kg/m^3)
	///	</summary>
	inline void SetRhoWater(double dRhoWater) {
		m_dRhoWater = dRhoWater;
	}

	///	<summary>
	///		Set the ideal gas constant of water vapor (J/kg/K)
	///	</summary>
	inline void SetRvap(double dRvap) {
		m_dRvap = dRvap;
	}

	///	<summary>
	///		Set the ratio of molar mass of dry air to molar mass of water
	///		vapor (dimensionless)
	///	</summary>
	inline void SetMvap(double dMvap) {
		m_dMvap = dMvap;
	}

	///	<summary>
	///		Set the latent heat of vaporization of water vapor (J/kg)
	///	</summary>
	inline void SetLvap(double dLvap) {
		m_dLvap = dLvap;
	}

private:
	///	<summary>
	///		Recalculate kappa due to changes in R or Cp.
	///	</summary>
	inline void RecalculateKappa() {
		m_dKappa = m_dR / m_dCp;
	}

	///	<summary>
	///		Recalculate gamma due to changes in R or Cp.
	///	</summary>
	inline void RecalculateGamma() {
		m_dGamma = m_dCp / (m_dCp - m_dR);
	}

	///	<summary>
	///		Recalculate pressure scaling due to changes in R, Cp or P0.
	///	</summary>
	inline void RecalculatePressureScaling() {
		m_dPressureScaling = m_dP0 * pow(m_dR / m_dP0, m_dGamma);
	}

public:
	///	<summary>
	///		Calculate the pressure from potential temperature density.
	///	</summary>
	inline double PressureFromRhoTheta(double dRhoTheta) const {
		return m_dPressureScaling * exp(log(dRhoTheta) * m_dGamma);
	}

	///	<summary>
	///		Calculate the potential temperature density from the pressure.
	///	</summary>
	inline double RhoThetaFromPressure(double dPressure) const {
		return exp(log(dPressure / m_dPressureScaling) / m_dGamma);
	}

public:
	///	<summary>
	///		Calculate the Exner pressure from potential temperature density.
	///	</summary>
	inline double ExnerPressureFromRhoTheta(double dRhoTheta) const {
		return m_dCp * exp(m_dR / (m_dCp - m_dR) * log(m_dR / m_dP0 * dRhoTheta));
	}

	///	<summary>
	///		Calculate the potential temperature density from the Exner
	///		pressure.
	///	</summary>
	inline double RhoThetaFromExnerPressure(double dPi) const {
		return m_dP0 / m_dR * exp((m_dCp - m_dR) / m_dR * log(dPi / m_dCp));
	}

public:
	///	<summary>
	///		Calculate the Exner pressure from pressure.
	///	</summary>
	inline double ExnerPressureFromPressure(double dP) const {
		return m_dCp * exp(m_dR / m_dCp * log(dP / m_dP0));
	}

	///	<summary>
	///		Calculate the pressure from Exner pressure.
	///	</summary>
	inline double PressureFromExnerPressure(double dPi) const {
		return m_dP0 * exp(m_dCp / m_dR * log(dPi / m_dCp));
	}

#ifdef TEMPEST_NETCDF
public:
	///	<summary>
	///		Output the set of physical constants to a NetCDF file.
	///	</summary>
	///	<param name="ncOut">
	///		NetCDF file to receive the output.
	///	</param>
	void NcOutputPhysicalConstants(NcFile * ncOut);

	///	<summary>
	///		Input the set of physical constants from a NetCDF file.
	///	</summary>
	///	<param name="ncIn">
	///		NetCDF file that contains the set of physical constants.
	///	</param>
	void NcInputPhysicalConstants(NcFile * ncIn);

#else
public:
	///	<summary>
	///		Stub function to avoid .cpp has no symbols error.
	///	</summary>
	void Stub();
#endif
};

///////////////////////////////////////////////////////////////////////////////

#endif

