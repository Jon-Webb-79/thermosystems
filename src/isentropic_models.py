# Import necessary packages here
from scipy.optimize import bisect
from typing import Dict, Tuple
from src.therm_prop import ThermProps

# ============================================================================
# ============================================================================
# Date:    September 17, 2020
# Purpose: This file contains classes describing the performance of
#          gas flow systems as classes

# Source Code Metadata
__author__ = "Jonathan A. Webb"
__copyright__ = "Copyright 2020, Jon Webb Inc."
__version__ = "1.0"

# REFERENCES
#
# 1. Hill, P. and Peterson C., "Mechanics and Thermodynamics of Propulsion,"
#    Addison Wesley Publishing Company, Reading, MA, 1992
# ============================================================================
# ============================================================================
# This class describes stagnation conditions based on static conditions


class Stagnation:
    """
    This class contains equations for generic adiabatic-isentropic
    relationships that describe the effect of work or heat addition
    on the properties of a gaseous flow field.  All functions are
    encoded as class methods and all relationships are derived from
    the following references

    1. Hill, P. and Peterson, C., "Mechanics and Thermodynamics of Propulsion,"
       Addison Wesley Publishing Co., Reading, MA, 1992
    """
    @classmethod
    def stagnation_pressure(cls, static_pressure: float, mach_number: float,
                            gamma: float) -> float:
        """

        :param static_pressure: The static pressure in units of Pascals
        :param mach_number: The Mach number
        :param gamma: Ratio of specific heats
        :return stag_pres: The stagnation pressure in units of Pascals

        This function determines the stagnation pressure
        using the relationship from pg. 72 or reference 1 shown
        below

        .. math::
           P_o = P\\left[1+\\frac{\gamma-1}{\gamma}M^2\\right]^{\gamma/\\left(\gamma-1\\right)}

        """
        stag_pres = static_pressure * (1.0 + ((gamma - 1.0) / 2.0) *
                                       mach_number ** 2.0) ** (gamma / (gamma - 1.0))
        return stag_pres
# ----------------------------------------------------------------------------

    @classmethod
    def stagnation_temperature(cls, static_temperature: float,
                               mach_number: float,
                               gamma: float) -> float:
        """

        :param static_temperature: The static temperature in units of Kelvins
        :param mach_number: The Mach number
        :param gamma: Ratio of specific heats
        :return stag_temp: The stagnation temperature

        This function determines the stagnation temperature
        using the relationship from pg. 72 or reference 1,
        shown below

        .. math::
           T_o=T\\left(1+\\frac{\gamma-1}{2}M^2\\right)
        """
        stag_temp = static_temperature * (1.0 + ((gamma - 1) / 2.0) *
                                          mach_number ** 2.0)
        return stag_temp
# ----------------------------------------------------------------------------

    @classmethod
    def static_pressure(cls, stagnation_pressure: float,
                        mach_number: float,
                        gamma: float) -> float:
        """

        :param stagnation_pressure: The stagnation pressure in units
                                       of Pascals
        :param mach_number: The Mach number
        :param gamma: Ratio of specific heats
        :return static_pres: The static pressure in units of Pascals

        This function determines the static pressure
        using the relationship from pg. 72 or reference 1, shown
        below

        .. math::
           P = \\frac{P_o}{\\left[1+\\frac{\gamma-1}{2}M^2\\right]^{\gamma/\\left(\gamma-1\\right)}}
        """
        denominator = (1.0 + ((gamma - 1.0) / 2.0) * mach_number ** 2.0)
        denominator = denominator ** (gamma / (gamma - 1.0))
        return stagnation_pressure / denominator
# ----------------------------------------------------------------------------

    @classmethod
    def static_temperature(cls, stagnation_temperature: float,
                           mach_number: float,
                           gamma: float) -> float:
        """

        :param stagnation_temperature: The stagnation temperature in units
                                       of Kelvins
        :param mach_number: The Mach number
        :param gamma: Ratio of specific heats
        :return static_temp: The static temperature in units of Kelvins

        This function determines the static temperature
        using the relationship from pg. 72 or reference 1,
        shown below

        .. math::
           T = \\frac{T_o}{\\left(1+\\frac{\gamma-1}{2}M^2\\right)}
        """
        denominator = 1.0 + ((gamma - 1.0) / 2.0) * mach_number ** 2.0
        return stagnation_temperature / denominator
# ----------------------------------------------------------------------------

    @classmethod
    def pressure_based_mach(cls, stag_pres: float, stat_pres: float,
                            gamma: float) -> float:
        """

        :param stag_pres: The stagnation pressure in units of Pascals
        :param stat_pres: The static pressure in units of Pascals
        :param gamma: The ratio of specific heats
        :return mach_number: The Mach number

        This function determines the mach number
        using the relationship from pg. 72 or reference 1,
        shown below.

        .. math::
           M = \sqrt[]{\\frac{\\left(P_o/P\\right)^{\\left(\gamma-1\\right)/\gamma} - 1}{\\frac{\gamma-1}{2}}}
        """
        one = ((stag_pres / stat_pres) ** ((gamma - 1.0) / gamma)) - 1.0
        return ((2 * one) / (gamma - 1.0)) ** 0.5
# ----------------------------------------------------------------------------

    @classmethod
    def temperature_based_mach(cls, stag_temp: float, stat_temp: float,
                               gamma: float) -> float:
        """

        :param stag_temp: The stagnation temperature in units of Kelvins
        :param stat_temp: The static temperature in units of Kelvins
        :param gamma: The ratio of specific heats
        :return mach_number: The Mach number

        This function determines the mach number
        using the relationship from pg. 72 or reference 1,
        shown below

        .. math::
           M = \sqrt[]{\\frac{\\left(T_o/T\\right) - 1}{\\frac{\gamma-1}{2}}}
        """
        one = (stag_temp / stat_temp) - 1.0
        two = one / ((gamma - 1.0) / 2.0)
        return two ** 0.5
# ----------------------------------------------------------------------------

    @classmethod
    def speed_of_sound(cls, gamma: float, static_temperature: float,
                       molar_mass: float) -> float:
        """

        :param gamma: The ratio of specific heat
        :param static_temperature: The static temperature in units of Kelvins
        :param molar_mass: Molar mass of the fluid in units of g/mol
        :return sos: The speed of sound in meters per second

        This function determines the speed of sound
        using the relationship from pg. 68 or reference 1,
        shown below

        .. math::
           a = \sqrt[]{\gamma R T}
        """
        universal_gas_constant = 8.314  # J/K-mol
        gas_constant = 1000.0 * (universal_gas_constant / molar_mass)
        return (gamma * gas_constant * static_temperature) ** 0.5
# ----------------------------------------------------------------------------

    @classmethod
    def mach_number(cls, gamma: float, static_temperature: float,
                    molar_mass: float, velocity: float) -> float:
        """

        :param gamma: The ratio of specific heats
        :param static_temperature: The static temperature in units of Kelvins
        :param molar_mass: The molar mass of the fluid in units of J/mol-K
        :param velocity: The magnitude of the velocity in units of
                         meters per second
        :return mach_number: The Mach number

        This function determines the mach number
        using the relationship from pg. 68 or reference 1, shown below

        .. math::
           M = \\frac{u}{\sqrt[]{\gamma R T}}
        """
        sos = Stagnation.speed_of_sound(gamma, static_temperature, molar_mass)
        return velocity / sos
# ----------------------------------------------------------------------------

    @classmethod
    def temp_exit_mach(cls, exit_mach: float, inlet_mach: float,
                       gamma: float, stag_inlet_temp: float,
                       stag_exit_temp: float) -> float:
        """

        :param exit_mach: Mach number at exit
        :param inlet_mach: Mach number at inlet
        :param gamma: Ratio of specific heats
        :param stag_inlet_temp: Stagnation temperature at inlet in units
                                of Kelvins
        :param stag_exit_temp: Stagnation temperature at exit in units
                               of Kelvins
        :return delta: The difference between the right hand and left
                       hand side of the equation

        This function assumes energy addition or loss between the
        entrance and exit of a component.  This function is to be
        iteratively solved for determine the inlet mach number with
        knowledge of the exit mach number and the stagnation temperatures
        at both ends of the component.  A component could be a combustion
        chamber, or any component with isentropic losses.  This equation
        was derived from page 75 or Ref. 1., shown below

        .. math::
           \\frac{T_{o1}}{T_{o2}} = \\left[\\frac{1 + \gamma M^2_2}{1 + \gamma M^2_1}
           \\left(\\frac{M_1}{M_2}\\right)\\right]^2\\left(\\frac{1 + \\frac{\gamma-1}{2}M^2_1}
           {1 + \\frac{\gamma-1}{2}M^2_2}\\right)
        """
        one = (1.0 + gamma * exit_mach ** 2.0) / (1.0 + gamma * inlet_mach ** 2.0)
        two = (one * (inlet_mach / exit_mach)) ** 2.0
        three = 1.0 + ((gamma - 1.0) / 2.0) * inlet_mach ** 2.0
        four = 1.0 + ((gamma - 1.0) / 2.0) * exit_mach ** 2.0
        five = three / four
        return two * five - (stag_inlet_temp / stag_exit_temp)
# ----------------------------------------------------------------------------

    @classmethod
    def press_exit_mach(cls, exit_mach: float, inlet_mach: float,
                        gamma: float, stag_inlet_pres: float,
                        stag_exit_pres: float) -> float:
        """

        :param exit_mach: Mach number at exit
        :param inlet_mach: Mach number at inlet
        :param gamma: Ratio of specific heats
        :param stag_inlet_pres: Stagnation pressure at inlet, in units
                                of Pascals
        :param stag_exit_pres: Stagnation pressure at exit, in units of
                               Pascals
        :return delta: The difference between the right hand and left
                       hand side of the equation

        This function assumes energy addition or loss between the
        entrance and exit of a component.  This function is to be
        iteratively solved for determine the inlet mach number with
        knowledge of the exit mach number and the stagnation pressures
        at both ends of the component.  A component could be a combustion
        chamber, or any component with isentropic losses.  This equation
        was derived from page 75 or Ref. 1., shown below

        .. math::
           \\frac{P_{o1}}{P_{o2}} = \\left[\\frac{1 + \gamma M^2_2}{1 + \gamma M^2_1}
           \\right]\\left(\\frac{1 + \\frac{\gamma-1}{2}M^2_1}
           {1 + \\frac{\gamma-1}{2}M^2_2}\\right)^{\gamma/\\left(\gamma-1\\right)}
        """
        one = (1.0 + gamma * exit_mach ** 2.0) / (1.0 + gamma * inlet_mach ** 2.0)
        two = 1.0 + ((gamma - 1.0) / 2.0) * inlet_mach ** 2.0
        three = 1.0 + ((gamma - 1.0) / 2.0) * exit_mach ** 2.0
        four = (two / three) ** (gamma / (gamma - 1.0))
        return (one * four) - (stag_inlet_pres / stag_exit_pres)
# ============================================================================
# ============================================================================
# This class describes the performance of a diffuser


class DiffuserNozzle(Stagnation):
    """
    This class contains functions that determine the exit conditions for
    a diffuser or nozzle assuming knowledge of the inlet properties.  These
    relationships are developed for sub-sonic, incompressible flows and
    assume a negligible change the ratio of specific heats throughout the
    process.  All relationships are developed from the following references.

    1. Hill, P. and Peterson, C., "Mechanics and Thermodynamics of Propulsion,"
       Addison Wesley Publishing Co., Reading, MA, 1992
    """
    def __init__(self, efficiency: float, inlet_area: float, exit_area: float):
        """

        :param efficiency: The isentropic efficiency of the diffuser or nozzle
        :param inlet_area: The diffuser or nozzle inlet area in units of
                           square meters
        :param exit_area: The diffuser or nozzle exit area in units of
                          square meters
        """
        self.efficiency = efficiency
        self.inlet_area = inlet_area
        self.exit_area = exit_area
# ----------------------------------------------------------------------------

    def exit_stagnation_pressure(self, gamma: float,
                                 inlet_static_pres: float,
                                 inlet_mach_number) -> float:
        """

        :param gamma: Ratio of specific heats
        :param inlet_static_pres: The stagnation tempressure at the entrance
                                  to the diffuser or nozzle, in units of Pa
        :param inlet_mach_number: The Mach number at the diffuser or
                                  nozzle inlet
        :return exit_stag_pres: The stagnation pressure at the exit of the
                                diffuser or nozzle, in units of Pa

        This function determines the diffuser or nozzle outlet stagnation
        temperature using the relationship from pg. 225 or reference 1,
        shown below

        .. math::
           P_{o}=\\left[\eta M^2_{inlet}\\frac{\gamma - 1}{2}\\right]^{\gamma / \\left(\gamma - 1\\right)}
        """
        gamma_factor = gamma / (gamma - 1.0)
        exit_stag_pres = ((self.efficiency * ((gamma - 1.0) / 2.0) *
                          inlet_mach_number ** 2.0) + 1.0) ** gamma_factor
        return exit_stag_pres * inlet_static_pres
# ----------------------------------------------------------------------------

    def exit_stagnation_temperature(self, gamma: float, inlet_static_temp: float,
                                    inlet_static_pres: float,
                                    inlet_mach_number: float) -> float:
        """

        :param gamma: The ratio of specific heats
        :param inlet_static_temp: The diffuser or nozzle inlet static
                                  temperature in units of Kelvins
        :param inlet_static_pres: The diffuser or nozzle inlet static
                                  pressure in units of Pascals
        :param inlet_mach_number: The diffuser or nozzle inlet Mach number
        :return exit_stag_temp: The diffuser or nozzle exit stagnation
                                temperature in units of Kelvins

        This function determines the diffuser or nozzle  outlet stagnation
        temperature using the relationship from pg. 225 or reference 1,
        shown below

        .. math::
           T_{o2} = T_{o1}\\left(\\frac{P_{o2}}{P_{o1}}\\right)^{\\frac{\gamma-1}{\gamma}}
        """
        exit_stag_pres = self.exit_stagnation_pressure(gamma, inlet_static_pres,
                                                       inlet_mach_number)
        ratio = (exit_stag_pres / inlet_static_pres) ** ((gamma - 1.0) / gamma)
        return inlet_static_temp * ratio
# ----------------------------------------------------------------------------

    def exit_velocity(self, inlet_velocity: float) -> float:
        """

        :param inlet_velocity: The diffuser inlet velocity in units of meters
                               per second
        :return exit_velocity: The diffuser exit velocity in units of meters

        This function determines the fluid velocity leaving the diffuser
        via equation a re-arrangement of the law of continuity described
        on page 28 of Ref. 1, shown below

        .. math::
           u_2 = \\frac{u_1A_1}{A_2}
        """
        return inlet_velocity * self.inlet_area / self.exit_area
# ----------------------------------------------------------------------------

    def exit_static_temperature(self, specific_heat: float, gamma: float,
                                inlet_static_temp: float,
                                inlet_static_pres: float,
                                inlet_mach_number: float,
                                inlet_velocity: float) -> float:
        """

        :param specific_heat: The specific heat at constant pressure in units
                              of Joules/kg-K
        :param gamma: Ratio of specific heats
        :param inlet_static_temp: The static temperature at the inlet
                                  of the diffuser or nozzle in units of
                                  Kelvins
        :param inlet_static_pres: The static pressure at the inlet of the
                                  diffuser or nozzle in units of Pascals
        :param inlet_mach_number: The Mach number at the inlet of the
                                  diffuser or nozzle
        :param inlet_velocity: The velocity at the inlet to the diffuser or
                               nozzle
        :return static_temp: The static temperature at the diffuser or nozzle
                             exit, in units of Kelvins

        This function determines the fluid static temperature leaving the
        diffuser or nozzle.  The equation is derived from Eq. 3.9 on page
        70 of Ref. 1, shown below.

        .. math::
           T_2 = T_{o2}-\\frac{u^2}{2c_p}
        """
        stag_temp = self.exit_stagnation_temperature(gamma, inlet_static_temp,
                                                     inlet_static_pres,
                                                     inlet_mach_number)
        velocity = self.exit_velocity(inlet_velocity)
        return stag_temp - ((velocity ** 2) / (2 * specific_heat))
# ----------------------------------------------------------------------------

    def exit_mach_number(self, specific_heat: float, gamma: float,
                         inlet_static_temp: float, inlet_static_pres: float,
                         inlet_mach_number: float, inlet_velocity: float,
                         molar_mass: float) -> float:
        """

        :param specific_heat: The specific heat at constant pressure in units
                              of Joules/kg-K
        :param gamma: Ratio of specific heats
        :param inlet_static_temp: The static temperature at the inlet to the
                                  diffuser or nozzle, in units of Kelvins
        :param inlet_static_pres: The static pressure at the inlet to the
                                  diffuser or nozzle, in units of Pascals
        :param inlet_mach_number: The Mach number at the inlet to the diffuser
                                  or nozzle
        :param inlet_velocity: The velocity at the inlet to the diffuser or
                               nozzle, in units of m/s
        :param molar_mass: The molar mass of the fluid in units of J/mol-K
        :return mach_number: The Mach number at the exit to the diffuser or
                             nozzle

        This function calculates the Mach number at the diffuser or nozzle
        exit.  This function assumes sub-sonic, incompressible flow and uses
        Eq. 3.7 and 3.8 from page 69 of Ref. 1, shown below

        .. math::
           M_2 = \\frac{u_2}{\sqrt[]{\gamma R T_2}}
        """
        stat_temp = self.exit_static_temperature(specific_heat, gamma,
                                                 inlet_static_temp,
                                                 inlet_static_pres,
                                                 inlet_mach_number,
                                                 inlet_velocity)
        velocity = self.exit_velocity(inlet_velocity)
        mach_number = self.mach_number(gamma, stat_temp, molar_mass, velocity)
        return mach_number
# ----------------------------------------------------------------------------

    def exit_static_pressure(self, gamma: float, inlet_static_pres: float,
                             inlet_mach_number, specific_heat: float,
                             inlet_static_temp: float, inlet_velocity: float,
                             molar_mass: float) -> float:
        """

        :param gamma: Ratio of specific heats
        :param inlet_static_pres: The static temperature at the inlet to
                                  the diffuser or nozzle, in units of Kelvins
        :param inlet_mach_number: The Mach number at the diffuser/nozzle
                                  inlet
        :param specific_heat: The specific heat at constant pressure in units
                              of Joules/kg-K
        :param inlet_static_temp: The static temperature at the inlet
                                  of the diffuser or nozzle, in units of
                                  Kelvins
        :param inlet_velocity: The velocity at the inlet to the diffuser or
                               nozzle in units of m/s
        :param molar_mass: The molar mass of the fluid in units of J/mol-K
        :return stat_pres: The static pressure at the exit to the diffuser
                           or nozzle in units of Pascals

        This function calculates the static pressure at the exit of the
        diffuser or nozzle.  This function assumes sub-sonic, incompressible
        flow, assuming the following relationship

        .. math::
           P = \\frac{P_o}{\\left[1+\\frac{\gamma-1}{2}M^2\\right]^{\gamma/\\left(\gamma-1\\right)}}
        """
        stag_pres = self.exit_stagnation_pressure(gamma, inlet_static_pres,
                                                  inlet_mach_number)
        mach_num = self.exit_mach_number(specific_heat, gamma,
                                         inlet_static_temp, inlet_static_pres,
                                         inlet_mach_number, inlet_velocity,
                                         molar_mass)
        stat_pres = self.static_pressure(stag_pres, mach_num, gamma)
        return stat_pres
# ----------------------------------------------------------------------------

    def exit_density(self, inlet_velocity: float, mass_flow_rate: float) -> float:
        """

        :param inlet_velocity: The velocity at the inlet to the diffuser or
                               nozzle, in units of m/s
        :param mass_flow_rate: The fluid mass flow rate in units of
                               kg/s
        :return exit_density: The fluid density at the diffuser/nozzle exit
                              im units of kg/s

        This function calculates the fluid density at the ext to the diffuser
        or nozzle with the following relationship.

        .. math::
           \\rho = \\frac{\dot{m}}{uA}
        """
        velocity = self.exit_velocity(inlet_velocity)
        return mass_flow_rate / (velocity * self.exit_area)
# ============================================================================
# ============================================================================
# This class describes the performance of a compressor


class Compressor(Stagnation):
    """
    This class contains functions that determine the exit conditions for
    a gas compressor assuming knowledge of the inlet properties.  These
    relationships are developed for sub-sonic, incompressible flows and
    assume a negligible change the ratio of specific heats throughout the
    process.  All relationships are developed from the following references.

    1. Hill, P. and Peterson, C., "Mechanics and Thermodynamics of Propulsion,"
       Addison Wesley Publishing Co., Reading, MA, 1992
    """
    def __init__(self, compression_ratio: float, efficiency: float):
        """

        :param compression_ratio: The ratio of outlet stagnation pressure
                                  to inlet stagnation pressure
        :param efficiency: The compressor isentropic efficiency
        """
        self.compression_ratio = compression_ratio
        self.efficiency = efficiency
# ----------------------------------------------------------------------------

    def exit_stagnation_pressure(self, inlet_stagnation_pressure: float) -> float:
        """

        :param inlet_stagnation_pressure: The inlet stagnation pressure in units
                               of Pascals
        :return outlet_stagnation_pressure: The stagnation pressure
                                            leaving the compressure in
                                            units of Pascals

        This function calculates the stagnation pressure leaving the compressor
        using Equation 5.30 from page 161 of Ref. 1, shown below

        .. math::
           p_{cr}=\\frac{P_{oexit}}{P_{oinlet}}

        """
        return self.compression_ratio * inlet_stagnation_pressure
# ----------------------------------------------------------------------------

    def exit_stagnation_temperature(self, inlet_stagnation_temperature: float,
                                    gamma: float) -> float:
        """

        :param inlet_stagnation_temperature: The stagnation temperature at the
                                             compressor inlet in units of Kelvins
        :param gamma: Ratio of specific heats
        :return exit_stag_temp: The stagnation temperature at the exit of the
                                compressor in units of Kelvins

        This function calculates the stagnation temperature leaving the
        compressor using Equation 5.40 from page 171 of Ref. 1 shown
        below

        .. math::
           T_{oexit} = T_{oinlet}\\left[1 + \\frac{1}{\eta_c}
           \\left(p_{cr}^{\\frac{\gamma-1}{\gamma}}\\right) - 1\\right]
        """
        stag_temp = inlet_stagnation_temperature * ((1.0 + (1.0 / self.efficiency) *
                    (self.compression_ratio ** ((gamma - 1.0) / gamma) - 1.0)))
        return stag_temp
# ----------------------------------------------------------------------------

    def compressor_work(self, mass_flow_rate: float, specific_heat: float,
                        inlet_stagnation_temperature: float, gamma:float) -> float:
        """

        :param mass_flow_rate: The mass flow rate in units of kg/s
        :param specific_heat: The specific heat at constant pressure,
                              in units of J/kg-K
        :param inlet_stagnation_temperature: The stagnation temperature at
                                             the compressor inlet
        :param gamma: Ratio of specific heats
        :return work: The work done on the fluid in units of Watts

        This function calculates the work done on the fluid by a compressor
        using Equation 5.43 on page 172 of Ref 1 shown below.

        .. math::
           W = \\dot{m}c_p\\left(T_{o2} - T_{o1}\\right)
        """
        exit_stag_temp = self.exit_stagnation_temperature(inlet_stagnation_temperature,
                                                          gamma)
        work = mass_flow_rate * specific_heat * (exit_stag_temp -
                                                 inlet_stagnation_temperature)
        return work / self.efficiency
# ----------------------------------------------------------------------------

    def exit_mach_number(self, inlet_mach_number: float,
                         inlet_stagnation_temperature: float,
                         gamma: float, upper: float = 0.001,
                         lower: float = 0.999) -> float:
        """

        :param inlet_mach_number: The mach number at the compressor inlet
        :param inlet_stagnation_temperature: The stagnation temperature at the
                                             compressor inlet in units of
                                             Kelvins
        :param gamma: Ratio of specific heats
        :param upper: The upper bounds of the possible solution for Mach number,
                      defaulted to 0.999
        :param lower: The lower bounds of the possible solution for Mach number,
                      defaulted to 0.001
        :return mach: The mach number at the compressor exit

        This function determines the Mach number at the compressor exit using
        Eq. 3.19 from page 75 of Ref. 1 as shown below.  The function cannot directly
        solve for the exit Mach number and instead uses the bisect method to
        solve for the root value.

        .. math::
           \\frac{T_{o2}}{T_{o1}} = \\left[\\frac{1+\gamma M^2_2}{1+\gamma M^2_1}
           \\left(\\frac{M_2}{M_1} \\right)\\right]^2 \\left(\\frac{1+\\frac{\gamma-1}{2}M^2_2}
           {1+\\frac{\gamma-1}{2}M^2_1} \\right)
        """
        exit_stag_temp = self.exit_stagnation_temperature(inlet_stagnation_temperature,
                                                          gamma)
        mach = bisect(self.temp_exit_mach, lower, upper,
                      args=(inlet_mach_number, gamma, inlet_stagnation_temperature,
                            exit_stag_temp), full_output=True)
        return mach[0]
# ----------------------------------------------------------------------------

    def exit_static_pressure(self, inlet_stagnation_pressure: float,
                             inlet_stagnation_temperature: float,
                             gamma: float, inlet_mach_number: float) -> float:
        """

        :param inlet_stagnation_pressure: The stagnation pressure at the
                                          compressor entrance in units of Pascals
        :param inlet_stagnation_temperature: The stagnation temperature at the
                                             compressor entrance in unites of Kelvins
        :param gamma: The ratio of specific heats
        :param inlet_mach_number: The Mach number at the compressor entrance
        :return stat_press: The static pressure at the exit in units of Pascals

        This function solves for the compressor exit static pressure using
        the relationship from pg. 72 or reference 1, shown
        below

        .. math::
           P = \\frac{P_o}{\\left[1+\\frac{\gamma-1}{2}M^2\\right]^{\gamma/\\left(\gamma-1\\right)}}
        """
        mach_number = self.exit_mach_number(inlet_mach_number,
                                            inlet_stagnation_temperature,
                                            gamma)
        exit_stag_pres = self.exit_stagnation_pressure(inlet_stagnation_pressure)
        stat_press = self.static_pressure(exit_stag_pres,
                                          mach_number, gamma)
        return stat_press
# ----------------------------------------------------------------------------

    def exit_static_temperature(self, inlet_stagnation_temperature: float,
                                gamma: float, inlet_mach_number: float) -> float:
        """

        :param inlet_stagnation_temperature: The stagnation temperature at the
                                             compressor inlet in units of Kelvins
        :param gamma: The ratio of specific heats
        :param inlet_mach_number: The Mach number at the inlet to the compressor
        :return stat_temp: The static temperature at the outlet of the compressor

        This function determines the static temperature
        using the relationship from pg. 72 or reference 1,
        shown below

        .. math::
           T = \\frac{T_o}{\\left(1+\\frac{\gamma-1}{2}M^2\\right)}
        """
        mach_number = self.exit_mach_number(inlet_mach_number,
                                            inlet_stagnation_temperature,
                                            gamma)
        exit_stag_temp = self.exit_stagnation_temperature(inlet_stagnation_temperature,
                                                          gamma)
        stat_temp = self.static_temperature(exit_stag_temp, mach_number, gamma)
        return stat_temp
# ----------------------------------------------------------------------------

    def exit_velocity(self, inlet_stagnation_temperature: float, gamma: float,
                      inlet_mach_number: float, molar_mass: float):
        """

        :param inlet_stagnation_temperature: The stagnation temperature at the
                                             compressor inlet in units of Kelvins
        :param gamma: The ratio of specific heats
        :param inlet_mach_number: The Mach number at the compressor inlet
        :param molar_mass: The fluid molar mass in units of J/mol-K
        :return velocity: The fluid velocity in units of meters per second

        This function determines the fluid velocity via the relationship shown
        below.

        .. math::
           u = M\\sqrt[]{\gamma R T}
        """
        stat_temp = self.exit_static_temperature(inlet_stagnation_temperature,
                                                 gamma, inlet_mach_number)
        sos = self.speed_of_sound(gamma, stat_temp, molar_mass)
        mach_number = self.exit_mach_number(inlet_mach_number, inlet_stagnation_temperature, gamma)
        return mach_number * sos
# ============================================================================
# ============================================================================
# - This class describes the performance of a heat source such as would
#   be seen with a combustion chamber or heat exchanger


class HeatAddition(Stagnation):
    """
    This class contains functions that determine the exit conditions for
    a heat exchanger or combustion chamber with knowledge of the inlet properties.
    These relationships are developed for sub-sonic, incompressible flows and
    assume a negligible change the ratio of specific heats throughout the
    process.  All relationships are developed from the following references.

    1. Hill, P. and Peterson, C., "Mechanics and Thermodynamics of Propulsion,"
       Addison Wesley Publishing Co., Reading, MA, 1992
    """
    def __init__(self, efficiency: float):
        """

        :param efficiency: The isentropic efficiency of the heat exchanger
        """
        self.efficiency = efficiency
# ----------------------------------------------------------------------------

    def input_power(self, heat_addition: float) -> float:
        """

        :param heat_addition: The heat that must be added to the fluid in
                              units of Watts
        :return power: The power that must be supplied to the heat exchanger
                       in units of Watts

        This function determines the power that must be supplied to the heat
        exchanger in order to supply the correct heat addition to the fluid.
        This function determines the input power via the equation below.

        .. math::
           \dot{Q}_t=\\frac{\dot{Q}_{he}}{\eta_{he}}
        """
        return heat_addition / self.efficiency
# ----------------------------------------------------------------------------

    @classmethod
    def exit_stagnation_temperature(cls, inlet_stagnation_temperature: float,
                                    heat_addition: float, specific_heat: float,
                                    mass_flow_rate: float) -> float:
        """

        :param inlet_stagnation_temperature: The stagnation conditions at the
                                             inlet to the heat addition source
                                             in units of Kelvins
        :param heat_addition: The heat added to the fluid in units of Watts
        :param specific_heat: The specific heat of the fluid in units of J/kg-K
        :param mass_flow_rate: The fluid mass flow rate in units of kg/s
        :return stag_temp: The stagnation temperature at the exit to the
                           heat addition in units of Kelvins

        This function determines the stagnation temperature at the heat
        addition component.  This function solves teh following equation

        .. math::
           T_{o2} = T_{o1} + \\frac{\dot{Q}_{he}}{\dot{m}c_p}
        """
        stag_temp = inlet_stagnation_temperature + heat_addition / \
                    (mass_flow_rate * specific_heat)
        return stag_temp
# ----------------------------------------------------------------------------

    def exit_mach_number(self, inlet_stagnation_temperature: float,
                         heat_addition: float, specific_heat: float,
                         mass_flow_rate: float, inlet_mach_number: float,
                         gamma: float, lower: float = 0.001,
                         upper: float = 0.999) -> float:
        """

        :param inlet_stagnation_temperature: The stagnation temperature at the
                                             heat addition inlet in units of
                                             Kelvins
        :param heat_addition: The heat addition in units of Watts
        :param specific_heat: The fluid specific heat in units of J/kg-K
        :param mass_flow_rate: The mass flow rate in units of kg/s
        :param inlet_mach_number: The Mach number at the inlet to the
                                  heat addition
        :param gamma: The ratio of specific heats
        :param lower: The lower bound for the Mach number solution
        :param upper: The upper bound for the Mach number solution
        :return exit_mach: The Mach number at the exit to the heat
                           addition unit.

        This function solves for the Mach number at the exit of the
        heat addition unit.  This equation was derived from page
        75 or Ref. 1., shown below

        .. math::
           \\frac{T_{o1}}{T_{o2}} = \\left[\\frac{1 + \gamma M^2_2}{1 + \gamma M^2_1}
           \\left(\\frac{M_1}{M_2}\\right)\\right]^2\\left(\\frac{1 + \\frac{\gamma-1}{2}M^2_1}
           {1 + \\frac{\gamma-1}{2}M^2_2}\\right)
        """
        exit_stag_temp = self.exit_stagnation_temperature(inlet_stagnation_temperature,
                                                          heat_addition, specific_heat,
                                                          mass_flow_rate)
        exit_mach = bisect(self.temp_exit_mach, lower, upper,
                           args=(inlet_mach_number, gamma, inlet_stagnation_temperature,
                                 exit_stag_temp), full_output=True)
        return exit_mach[0]
# ----------------------------------------------------------------------------

    def exit_stagnation_pressure(self, inlet_stagnation_pressure: float,
                                 inlet_stagnation_temperature: float,
                                 heat_addition: float, specific_heat: float,
                                 mass_flow_rate: float, inlet_mach_number: float,
                                 gamma: float) -> float:
        """

        :param inlet_stagnation_pressure: The stagnation pressure at the inlet
                                          to the heat addition component in
                                          units of Pascals
        :param inlet_stagnation_temperature: The stagnation temperature at the
                                             heat addition inlet in units of
                                             Kelvins
        :param heat_addition: The heat addition in units of Watts
        :param specific_heat: The fluid specific heat in units of J/kg-K
        :param mass_flow_rate: The mass flow rate in units of kg/s
        :param inlet_mach_number: The Mach number at the inlet to the
                                  heat addition
        :param gamma: The ratio of specific heats
        :return exit_stag_pressure: The stagnation pressure at the exit
                                    to the heat addition component in
                                    units of Pascals

        This function determines the stagnation at the exit of the heat
        addition component.

        This equation
        was derived from page 75 of Ref. 1., shown below

        .. math::
           P_{o2} = P_{o1}\\left[\\frac{1 + \gamma M^2_1}{1 + \gamma M^2_2}
           \\right]\\left(\\frac{1 + \\frac{\gamma-1}{2}M^2_2}
           {1 + \\frac{\gamma-1}{2}M^2_1}\\right)^{\gamma/\\left(\gamma-1\\right)}
        """
        exit_mach = self.exit_mach_number(inlet_stagnation_temperature,
                                          heat_addition, specific_heat,
                                          mass_flow_rate, inlet_mach_number,
                                          gamma)
        one = (1.0 + gamma * exit_mach ** 2.0) / (1.0 + gamma * inlet_mach_number ** 2.0)
        two = 1.0 + ((gamma - 1.0) / 2.0) * inlet_mach_number ** 2.0
        three = 1.0 + ((gamma - 1.0) / 2.0) * exit_mach ** 2.0
        four = (two / three) ** (gamma / (gamma - 1.0))
        return inlet_stagnation_pressure / (one * four)
# ----------------------------------------------------------------------------

    def exit_static_temperature(self, inlet_stagnation_temperature: float,
                                heat_addition: float, specific_heat: float,
                                mass_flow_rate: float, gamma: float,
                                inlet_mach_number: float) -> float:
        """

        :param inlet_stagnation_temperature: The stagnation conditions at the
                                             heat addition component inlet in
                                             units of Kelvins
        :param heat_addition: The heat addition in units of Watts
        :param specific_heat: The fluid specific heat in units of J/kg-K
        :param mass_flow_rate: The mass flow rate in units of kg/s
        :param gamma: The ratio of specific heats
        :param inlet_mach_number: The Mach number at the heat addition
                                  component inlet
        :return exit_static_temp: The static temperature at the heat
                                  addition component in units of Kelvins

        This function determines the static temperature at the heat
        addition component exit.  This function determines the static
        temperature using the relationship from pg. 72 or reference 1,
        shown below

        .. math::
           T = \\frac{T_o}{\\left(1+\\frac{\gamma-1}{2}M^2\\right)}
        """
        exit_stag_temp = self.exit_stagnation_temperature(inlet_stagnation_temperature,
                                                          heat_addition, specific_heat,
                                                          mass_flow_rate)
        exit_mach_number = self.exit_mach_number(inlet_stagnation_temperature,
                                                 heat_addition,
                                                 specific_heat,
                                                 mass_flow_rate, inlet_mach_number,
                                                 gamma)
        exit_static_temp = self.static_temperature(exit_stag_temp,
                                                   exit_mach_number, gamma)
        return exit_static_temp
# ----------------------------------------------------------------------------

    def exit_static_pressure(self, inlet_stagnation_pressure: float,
                             inlet_stagnation_temperature: float,
                             heat_addition: float, specific_heat: float,
                             mass_flow_rate: float, gamma: float,
                             inlet_mach_number: float) -> float:
        """

        :param inlet_stagnation_pressure: The stagnation pressure at the inlet
                                          to the heat addition component in
                                          units of Pascals
        :param inlet_stagnation_temperature: The stagnation conditions at the
                                             heat addition component inlet in
                                             units of Kelvins
        :param heat_addition: The heat addition in units of Watts
        :param specific_heat: The fluid specific heat in units of J/kg-K
        :param mass_flow_rate: The mass flow rate in units of kg/s
        :param gamma: The ratio of specific heats
        :param inlet_mach_number: The Mach number at the heat addition
                                  component inlet
        :return static_pressure: The static pressure at the heat addition
                                 component inlet in units of Pascals

        This function determines the static pressure at the exit to the
        heat addition component. This function determines the static pressure
        using the relationship from pg. 72 or reference 1, shown
        below

        .. math::
           P = \\frac{P_o}{\\left[1+\\frac{\gamma-1}{2}M^2\\right]^{\gamma/\\left(\gamma-1\\right)}}
        """
        exit_mach_number = self.exit_mach_number(inlet_stagnation_temperature,
                                                 heat_addition, specific_heat,
                                                 mass_flow_rate, inlet_mach_number,
                                                 gamma)
        static_pressure = self.static_pressure(inlet_stagnation_pressure,
                                               exit_mach_number, gamma)
        return static_pressure
# ----------------------------------------------------------------------------

    def exit_velocity(self, inlet_stagnation_temperature: float,
                      heat_addition: float, specific_heat: float,
                      mass_flow_rate: float, inlet_mach_number: float,
                      gamma: float, molar_mass: float) -> float:
        """

        :param inlet_stagnation_temperature: The stagnation conditions at the
                                             heat addition component inlet in
                                             units of Kelvins
        :param heat_addition: The heat addition in units of Watts
        :param specific_heat: The fluid specific heat in units of J/kg-K
        :param mass_flow_rate: The mass flow rate in units of kg/s
        :param inlet_mach_number: The Mach number at the heat addition
                                  component inlet
        :param gamma: The ratio of specific heats
        :param molar_mass: Molar mass of the fluid in units of J/mol-K
        :return exit_velocity: The fluid velocity at the exit to the
                               heat addition component in units of
                               meters per second

        This function determines the fluid velocity leaving the heat
        addition component by solving the equation below

        .. math::
           u = M\sqrt[]{\gamma R T}
        """
        exit_mach = self.exit_mach_number(inlet_stagnation_temperature,
                                          heat_addition, specific_heat,
                                          mass_flow_rate, inlet_mach_number,
                                          gamma)
        static_temp = self.exit_static_temperature(inlet_stagnation_temperature,
                                                   heat_addition, specific_heat,
                                                   mass_flow_rate, gamma,
                                                   inlet_mach_number)
        sos = self.speed_of_sound(gamma, static_temp, molar_mass)
        return exit_mach * sos
# ============================================================================
# ============================================================================
# This class describes the performance of a compressor


class Turbine(Stagnation):
    """
    This class contains functions that determine the exit conditions for
    a gas turbine assuming knowledge of the inlet properties.  These
    relationships are developed for sub-sonic, incompressible flows and
    assume a negligible change the ratio of specific heats throughout the
    process.  All relationships are developed from the following references.

    1. Hill, P. and Peterson, C., "Mechanics and Thermodynamics of Propulsion,"
       Addison Wesley Publishing Co., Reading, MA, 1992

    This class can also be used to model a heat sink.
    """
    def __init__(self, efficiency: float):
        """

        :param efficiency: The isentropic efficiency of the turbine
        """
        self.efficiency = efficiency
# ----------------------------------------------------------------------------

    def work_extraction(self, turbine_work: float) -> float:
        """

        :param turbine_work: The usable work extracted by the turbine in
                             units of Watts
        :return total_work: The total work extracted by a turbine in
                            units of

        This function solves for the total work to be extracted from the
        system with knowledge of the work required for an alternator
        or a compressor.  This function solves the equation below.

        .. math::
           \dot{W}_{total}=\\frac{\dot{W}_t}{\eta_t}
        """
        return turbine_work / self.efficiency
# ----------------------------------------------------------------------------

    def exit_stagnation_temperature(self, turbine_work: float,
                                    specific_heat: float,
                                    mass_flow_rate: float,
                                    inlet_stagnation_temperature: float) -> float:
        """

        :param turbine_work: The total work to be extracted by the turbine
                             in units of Watts
        :param specific_heat: The fluid specific heat at constant pressure
                              in units of J/kg-K
        :param mass_flow_rate: The fluid mass flow rate in units of kg/s
        :param inlet_stagnation_temperature: The stagnation temperature at
                                             the turbine inlet in units
                                             of Kelvins
        :return stag_temp: The stagnation temperature at the turbine exit
                           in units of Kelvins

        This function calculates the stagnation temperature at the exit to
        the compressor with the equation below

        .. math::
           T_{o2} = T_{o1} - \\frac{\dot{W}_{total}}{\dot{m}c_p}
        """
        total_work = self.work_extraction(turbine_work)
        stag_temp = inlet_stagnation_temperature - (total_work /
                                                    (mass_flow_rate * specific_heat))
        return stag_temp
# ---------------------------------------------------------------------------

    def exit_stagnation_pressure(self, inlet_stagnation_pressure: float,
                                 inlet_stagnation_temperature: float,
                                 specific_heat: float, mass_flow_rate: float,
                                 turbine_work: float, gamma: float) -> float:
        """

        :param inlet_stagnation_pressure: The stagnation pressure at the turbine
                                          inlet in units of Pascals
        :param inlet_stagnation_temperature: The stagnation temperature at the
                                             turbine inlet in units of Kelvins
        :param specific_heat: The fluid specific heat in units of J/kg-K
        :param mass_flow_rate: The fluid mass flow rate in units of kg/s
        :param turbine_work: The usable work extracted by the turbine in
                             units of Watts
        :param gamma: The ratio of specific heats
        :return stag_pres: The stagnation pressure at the exit to the
                           diffuser in units of Pascals

        This function solves for the stagnation pressure at the turbine exit
        using the equation below.

        .. math::
           P_{o2} = P_{o1}\\left[1 - \\frac{1}{\eta_t}\\left(1 -
           \\frac{T_{o2}}{T_{o1}} \\right) \\right]^{\\frac{\gamma}{\gamma - 1}}
        """
        stag_temp = self.exit_stagnation_temperature(turbine_work, specific_heat,
                                                     mass_flow_rate, inlet_stagnation_temperature)
        one = 1 - (1/self.efficiency) * (1 - (stag_temp / inlet_stagnation_temperature))
        stag_pres = inlet_stagnation_pressure * one ** (gamma / (gamma - 1))
        return stag_pres
# ---------------------------------------------------------------------------

    def exit_mach_number(self, inlet_stagnation_temperature: float,
                         inlet_mach_number: float, gamma: float,
                         turbine_work: float, specific_heat: float,
                         mass_flow_rate: float, lower: float = 0.001,
                         upper: float = 0.999) -> float:
        """

        :param inlet_stagnation_temperature: The stagnation temperature at the
                                             turbine inlet in units of Kelvins
        :param inlet_mach_number: The Mach number at the inlet to the
                                  turbine
        :param gamma: The ratio of specific heats
        :param turbine_work: The usable work extracted by the turbine
        :param specific_heat: The fluid specific heat in units of J/kg-K
        :param mass_flow_rate: The fluid mass flow rate in units of kg/s
        :param lower: The lower bound for Mach number preset to 0.001
        :param upper: The upper bound for Mach number preset to 0.999
        :return exit_mach: The Mach number at the exit to the turbine

        This function determines the Mach number at the exit to the
        turbine using the equation shown below.


        .. math::
           \\frac{T_{o1}}{T_{o2}} = \\left[\\frac{1 + \gamma M^2_2}{1 + \gamma M^2_1}
           \\left(\\frac{M_1}{M_2}\\right)\\right]^2\\left(\\frac{1 + \\frac{\gamma-1}{2}M^2_1}
           {1 + \\frac{\gamma-1}{2}M^2_2}\\right)
        """
        exit_stag_temp = self.exit_stagnation_temperature(turbine_work, specific_heat,
                                                          mass_flow_rate,
                                                          inlet_stagnation_temperature)
        exit_mach = bisect(self.temp_exit_mach, lower, upper,
                           args=(inlet_mach_number, gamma, inlet_stagnation_temperature,
                                 exit_stag_temp), full_output=True)
        return exit_mach[0]
# ---------------------------------------------------------------------------

    def exit_static_temperature(self, inlet_stagnation_temperature: float,
                                inlet_mach_number: float, gamma: float,
                                turbine_work: float, specific_heat: float,
                                mass_flow_rate: float) -> float:
        """

        :param inlet_stagnation_temperature: The stagnation temperature at
                                             the inlet to the turbine in
                                             units of Kelvins
        :param inlet_mach_number: The Mach number at the inlet to the
                                  turbine
        :param gamma: The ratio of specific heats
        :param turbine_work: The usable work extracted by the turbine
        :param specific_heat: The fluid specific heat in units of J/kg-K
        :param mass_flow_rate: The mass flow rate in units of kg/s
        :return stat_temp: The static temperature at the exit to the
                           turbine in units of Kelvins.

        This function determines the static temperature at the exit to the
        turbine by solving the equation below.


        .. math::
           T = \\frac{T_o}{\\left(1+\\frac{\gamma-1}{2}M^2\\right)}
        """
        exit_mach = self.exit_mach_number(inlet_stagnation_temperature,
                                          inlet_mach_number, gamma, turbine_work,
                                          specific_heat, mass_flow_rate)
        stag_temp = self.exit_stagnation_temperature(turbine_work, specific_heat,
                                                     mass_flow_rate,
                                                     inlet_stagnation_temperature)
        stat_temp = self.static_temperature(stag_temp, exit_mach, gamma)
        return stat_temp
# ---------------------------------------------------------------------------

    def exit_static_pressure(self, inlet_stagnation_temperature: float,
                             inlet_mach_number: float, gamma: float,
                             turbine_work: float, specific_heat: float,
                             mass_flow_rate: float,
                             inlet_stagnation_pressure: float) -> float:
        """

        :param inlet_stagnation_temperature: The stagnation temperature at the
                                             inlet to the turbine in units of
                                             Kelvins
        :param inlet_mach_number: The Mach number at the inlet to the turbine
        :param gamma: The ratio of specific heats
        :param turbine_work: The usable work extracted by the turbine
        :param specific_heat: The fluid specific heat in units of J/kg-K
        :param mass_flow_rate: The fluid mass flow rate in units of kg/s
        :param inlet_stagnation_pressure: The stagnation pressure at the inlet
                                          to the turbine in units of Pascals
        :return stat_pres: The static pressure at the exit to the turbine
                           in units of Pascals.

        This function determines the static pressure at the turbine exit
        using the equation below.

        .. math::
           P = \\frac{P_o}{\\left[1+\\frac{\gamma-1}{2}M^2\\right]^{\gamma/\\left(\gamma-1\\right)}}
        """
        exit_mach = self.exit_mach_number(inlet_stagnation_temperature,
                                          inlet_mach_number, gamma, turbine_work,
                                          specific_heat, mass_flow_rate)
        stag_pres = self.exit_stagnation_pressure(inlet_stagnation_pressure,
                                                  inlet_stagnation_temperature,
                                                  specific_heat, mass_flow_rate,
                                                  turbine_work, gamma)
        stat_pres = self.static_pressure(stag_pres, exit_mach,
                                         gamma)
        return stat_pres
# ----------------------------------------------------------------------------

    def exit_velocity(self, inlet_stagnation_temperature: float,
                      inlet_mach_number: float, gamma: float,
                      turbine_work: float, specific_heat: float,
                      mass_flow_rate: float, molar_mass: float) -> float:
        """

        :param inlet_stagnation_temperature: The stagnation temperature at
                                             the turbine inlet in units of
                                             Kelvins
        :param inlet_mach_number: The Mach number at the inlet to the turbine
        :param gamma: The ratio of specific heats
        :param turbine_work: The usable work extracted by the turbine
        :param specific_heat: The fluid specific heat in units of J/kg-K
        :param mass_flow_rate: The fluid mass flow rate in units of kg/s
        :param molar_mass: The molar mass in units of J/mol-K
        :return velocity: The fluid velocity in meters per second

        This function determines the fluid velocity using the equation
        below

        .. math::
           u = M \sqrt[]{\gamma R T}
        """
        stat_temp = self.exit_static_temperature(inlet_stagnation_temperature,
                                                 inlet_mach_number, gamma, turbine_work,
                                                 specific_heat, mass_flow_rate)
        exit_mach = self.exit_mach_number(inlet_stagnation_temperature,
                                          inlet_mach_number, gamma, turbine_work,
                                          specific_heat, mass_flow_rate)
        sos = self.speed_of_sound(gamma, stat_temp, molar_mass)
        return sos * exit_mach
# ============================================================================
# ============================================================================


class Propeller(Stagnation):
    """
    This class contains functions that determine the exit conditions for
    a propeller assuming knowledge of the inlet properties.  These
    relationships are developed for sub-sonic, incompressible flows and
    assume a negligible change the ratio of specific heats throughout the
    process.  All relationships are developed from the following references.

    1. Hill, P. and Peterson, C., "Mechanics and Thermodynamics of Propulsion,"
       Addison Wesley Publishing Co., Reading, MA, 1992
    """
    def __init__(self, efficiency: float):
        """

        :param efficiency: The isentropic efficiency of the propeller
        """
        self.efficiency = efficiency
# ----------------------------------------------------------------------------

    def propeller_work(self, extracted_work: float) -> float:
        """

        :param extracted_work: The usable work extracted by the
                               propeller in units of Watts
        :return total_work: The total amount of work extracted by the
                            propeller in units of Watts.

        This function determines the total amount of work extracted
        by the propeller using the equation below

        .. math::
           \dot{W}_p = \\frac{\dot{W}}{\eta_p}
        """
        return extracted_work / self.efficiency
# ----------------------------------------------------------------------------

    def exit_stagnation_temperature(self, inlet_stagnation_temperature: float,
                                    mass_flow_rate: float, specific_heat: float,
                                    extracted_work: float) -> float:
        """

        :param inlet_stagnation_temperature: The stagnation temperature at
                                             the propeller inlet in units
                                             of Kelvins
        :param mass_flow_rate: The mass flow rate in units of kg/s
        :param specific_heat: The fluid specific heat at constant pressure
                              in units of J/kg-K
        :param extracted_work: The usable work extracted by the propeller
                               in units of Watts
        :return stag_temp: The stagnation temperature at the propeller
                           exit in units of Kelvins

        This function determines the stagnation temperature leaving
        the propeller using the equation below

        .. math::
           T_{o2} = T_{o1} + \\frac{\dot{W}_p}{\dot{m}c_p}
        """
        work = self.propeller_work(extracted_work)
        temp = inlet_stagnation_temperature + (work / (mass_flow_rate * specific_heat))
        return temp
# ----------------------------------------------------------------------------

    def exit_stagnation_pressure(self, inlet_stagnation_pressure: float,
                                 inlet_stagnation_temperature: float,
                                 gamma: float, mass_flow_rate: float,
                                 specific_heat: float, extracted_work: float) -> float:
        """

        :param inlet_stagnation_pressure: The stagnation pressure at the
                                          propeller inlet in units of Pascals
        :param inlet_stagnation_temperature: The stagnation temperature
                                             at the propeller inlet in
                                             units of Kelvins
        :param gamma: The ratio of specific heats
        :param mass_flow_rate: The fluid mass flow rate in units of kg/s
        :param specific_heat: The fluid specific heat at constant pressure
                              in units of J/kg-K
        :param extracted_work: The usable work extracted by the propeller
                               in units of Watts
        :return stag_pres: The stagnation pressure at the propeller
                           exit in units of Pascals

        This function determines the stagnation pressure at the propeller
        exit using the equation shown below;

        .. math::
           P_{o2} = P_{o1} \\left[1 + \eta_d\\left(\\frac{T_{o2}}{T_{o1}} - 1 \\right)
           \\right]^{\\frac{\gamma}{\gamma - 1}}
        """
        stag_temp = self.exit_stagnation_temperature(inlet_stagnation_temperature,
                                                     mass_flow_rate, specific_heat,
                                                     extracted_work)
        one = 1 + self.efficiency * ((inlet_stagnation_temperature / stag_temp) - 1.0)
        return inlet_stagnation_pressure * one ** (gamma / (gamma - 1.0))
# ----------------------------------------------------------------------------

    def exit_mach_number(self, inlet_stagnation_temperature: float,
                         inlet_mach_number, mass_flow_rate: float,
                         specific_heat: float, extracted_work: float,
                         gamma: float, lower: float = 0.001,
                         upper: float = 0.999) -> float:
        """

        :param inlet_stagnation_temperature: The stagnation temperature at the
                                             inlet to the propeller in units
                                             of Kelvins
        :param inlet_mach_number: The Mach number at the propeller inlet
        :param mass_flow_rate: The fluid mass flow rate in units of kg/s
        :param specific_heat: The fluid specific heat at constant pressure in
                              in units of J/kg-K
        :param extracted_work: The usable work extracted by the propeller in
                               units of Watts
        :param gamma: Ratio of specific heats
        :param lower: The lower bounds for Mach number
        :param upper: The upper bounds for Mach number
        :return mach_number: The Mach number at the propeller exit

        This equation determines the Mach number at the propeller exit.
        This equation
        was derived from page 75 or Ref. 1., shown below

        .. math::
           \\frac{T_{o1}}{T_{o2}} = \\left[\\frac{1 + \gamma M^2_2}{1 + \gamma M^2_1}
           \\left(\\frac{M_1}{M_2}\\right)\\right]^2\\left(\\frac{1 + \\frac{\gamma-1}{2}M^2_1}
           {1 + \\frac{\gamma-1}{2}M^2_2}\\right)
        """
        exit_temp = self.exit_stagnation_temperature(inlet_stagnation_temperature,
                                                     mass_flow_rate, specific_heat,
                                                     extracted_work)
        exit_mach = bisect(self.temp_exit_mach, lower, upper,
                           args=(inlet_mach_number, gamma, inlet_stagnation_temperature,
                                 exit_temp), full_output=True)
        return exit_mach[0]
# ----------------------------------------------------------------------------

    def exit_static_temperature(self, inlet_stagnation_temperature: float,
                                inlet_mach_number, mass_flow_rate: float,
                                specific_heat: float, extracted_work: float,
                                gamma: float) -> float:
        """

        :param inlet_stagnation_temperature: The stagnation temperature at the
                                             inlet to the propeller in units
                                             of Kelvins
        :param inlet_mach_number: The Mach number at the entrance to the
                                  propeller
        :param mass_flow_rate: The mass flow rate in units of kg/s
        :param specific_heat: The fluid specific heat in units of J/kg-K
        :param extracted_work: The usable work extracted by the propeller
                               in units of Watts
        :param gamma: The ratio of specific heats
        :return stat_temp: The static temperature at the propeller exit
                           in units of Kelvins

        This function determines the static temperature at the propeller
        exit using the relationship from pg. 72 or reference 1,
        shown below

        .. math::
           T = \\frac{T_o}{\\left(1+\\frac{\gamma-1}{2}M^2\\right)}
        """
        mach = self.exit_mach_number(inlet_stagnation_temperature, inlet_mach_number,
                                     mass_flow_rate, specific_heat, extracted_work,
                                     gamma)
        stag_temp = self.exit_stagnation_temperature(inlet_stagnation_temperature,
                                                     mass_flow_rate, specific_heat,
                                                     extracted_work)
        stat_temp = self.static_temperature(stag_temp, mach, gamma)
        return stat_temp
# ----------------------------------------------------------------------------

    def exit_static_pressure(self, inlet_stagnation_temperature: float,
                             inlet_mach_number, mass_flow_rate: float,
                             specific_heat: float, extracted_work: float,
                             gamma: float, inlet_stagnation_pressure: float):
        """

        :param inlet_stagnation_temperature: The stagnation temperature at the
                                             propeller inlet in units of Kelvins
        :param inlet_mach_number: The Mach number at the inlet to the propeller
        :param mass_flow_rate: The fluid mass flow rate in units of kg/s
        :param specific_heat: The fluid specific heat in units if J/kg-K
        :param extracted_work: The usable work extracted by the propeller in
                               units of Watts
        :param gamma: The ratio of specific heat
        :param inlet_stagnation_pressure: The stagnation pressure at the
                                          propeller inlet in units of Pascals
        :return stat_pressure: The stagnation pressure at the propeller
                               exit in units of Pascals

        This function determines the static pressure at the propeller exit
        using the relationship from pg. 72 or reference 1, shown
        below

        .. math::
           P = \\frac{P_o}{\\left[1+\\frac{\gamma-1}{2}M^2\\right]^{\gamma/\\left(\gamma-1\\right)}}
        """
        mach = self.exit_mach_number(inlet_stagnation_temperature, inlet_mach_number,
                                     mass_flow_rate, specific_heat, extracted_work,
                                     gamma)
        stag_pres = self.exit_stagnation_pressure(inlet_stagnation_pressure,
                                                  inlet_stagnation_temperature,
                                                  gamma, mass_flow_rate,
                                                  specific_heat, extracted_work)
        stat_pres = self.static_pressure(stag_pres, mach, gamma)
        return stat_pres
# ----------------------------------------------------------------------------

    def exit_velocity(self, inlet_stagnation_temperature: float,
                      inlet_mach_number, mass_flow_rate: float,
                      specific_heat: float, extracted_work: float,
                      gamma: float, molar_mass: float) -> float:
        """

        :param inlet_stagnation_temperature: The stagnation temperature at
                                             the propeller inlet in units
                                             of Kelvins
        :param inlet_mach_number: The Mach number at the propeller inlet
        :param mass_flow_rate: The mass flow rate in units of kg/s
        :param specific_heat: The fluid specific heat in units of J/kg-K
        :param extracted_work: The usable work extracted by the propeller
                               in units of Watts
        :param gamma: The ratio of specific heats
        :param molar_mass: The molar mass in units of J/mol-K
        :return velocity: The fluid velocity leaving the propeller
                          in units of meters per second

        This function determines the fluid velocity leaving the
        propeller using the following equation

        .. math::
           u = M\sqrt[]{\gamma R T}
        """
        stat_temp = self.exit_static_temperature(inlet_stagnation_temperature,
                                                 inlet_mach_number, mass_flow_rate,
                                                 specific_heat, extracted_work,
                                                 gamma)
        sos = self.speed_of_sound(gamma, stat_temp, molar_mass)
        mach = self.exit_mach_number(inlet_stagnation_temperature,
                                     inlet_mach_number, mass_flow_rate,
                                     specific_heat, extracted_work,
                                     gamma)
        return sos * mach
# ============================================================================
# ============================================================================


class DiffuserNozzleComponent(DiffuserNozzle):
    """

    This class combined all functionality of the DiffuserNozzle class
    into a single function for ease of user access
    """
    def __init__(self, efficiency: float, inlet_area: float, exit_area: float):
        """

        :param efficiency: The isentropic efficiency of the diffuser or nozzle
        :param inlet_area: The diffuser or nozzle inlet area in units of
                           square meters
        :param exit_area: The diffuser or nozzle exit area in units of
                          square meters
        """
        DiffuserNozzle.__init__(self, efficiency, inlet_area, exit_area)
# ----------------------------------------------------------------------------

    def outlet_conditions(self, gamma: float, specific_heat: float,
                          molar_mass: float,
                          inlet_static_temperature: float,
                          inlet_static_pressure: float,
                          inlet_mach_number: float, mass_flow_rate: float,
                          inlet_velocity: float) -> Dict[str, float]:
        """

        :param gamma: The ratio of specific heats
        :param specific_heat: The fluid specific heat at constant pressure
                              in units of J/kg-K
        :param molar_mass: The fluid molar mass in units of J/mol-K
        :param inlet_static_temperature: The inlet static temperature in
                                         units of Kelvins
        :param inlet_static_pressure: The inlet static pressure in units
                                      of Pascals
        :param inlet_mach_number: The inlet Mach number
        :param mass_flow_rate: The fluid mass flow in units of kg/s
        :param inlet_velocity: The inlet fluid velocity in units of
                               m/s
        :return dict: A dictionary containing all fluid exit
                      properties

        The dictionary for the diffuser contains the following key words

        * **static_temperature** in units of Kelvins
        * **static_pressure** in units of Pascals
        * **stagnation_temperature** in units of Kelvins
        * **stagnation_pressure** in units of Pascals
        * **velocity** in units of meters per second
        * **mach_number**
        * **density** in units of g/cc
        """
        exit_stag_pres = self.exit_stagnation_pressure(gamma, inlet_static_pressure,
                                                       inlet_mach_number)
        exit_stag_temp = self.exit_stagnation_temperature(gamma, inlet_static_temperature,
                                                          inlet_static_pressure,
                                                          inlet_mach_number)
        exit_velocity = self.exit_velocity(inlet_velocity)
        exit_stat_temp = self.exit_static_temperature(specific_heat, gamma,
                                                      inlet_static_temperature,
                                                      inlet_static_pressure,
                                                      inlet_mach_number,
                                                      inlet_velocity)
        exit_mach_number = self.exit_mach_number(specific_heat, gamma,
                                                 inlet_static_temperature,
                                                 inlet_static_pressure,
                                                 inlet_mach_number,
                                                 inlet_velocity, molar_mass)
        exit_stat_pres = self.exit_static_pressure(gamma, inlet_static_pressure,
                                                   inlet_mach_number, specific_heat,
                                                   inlet_static_temperature, inlet_velocity,
                                                   molar_mass)
        exit_density = self.exit_density(inlet_velocity, mass_flow_rate)
        dictionary = {'static_pressure': exit_stat_pres, 'static_temperature': exit_stat_temp,
                      'stagnation_pressure': exit_stag_pres, 'stagnation_temperature': exit_stag_temp,
                      'velocity': exit_velocity, 'mach_number': exit_mach_number,
                      'density': exit_density}
        return dictionary
# ============================================================================
# ============================================================================


class CompressorComponent(Compressor):
    """

    This class combined all functionality of the Compressor class
    into a single function for ease of user access
    """
    def __init__(self, compression_ratio: float, efficiency: float):
        """

        :param compression_ratio: The ratio of outlet stagnation pressure
                                  to the inlet stagnation pressure
        :param efficiency: The isentropic efficiency of the compressor
        """
        Compressor.__init__(self, compression_ratio, efficiency)
# ----------------------------------------------------------------------------

    def outlet_conditions(self, gamma: float, specific_heat: float,
                          molar_mass: float,
                          inlet_mach_number: float, mass_flow_rate: float,
                          inlet_stagnation_pressure: float,
                          inlet_stagnation_temperature: float) -> Dict[str, float]:
        """

        :param gamma: The ratio of specific heats
        :param specific_heat: The fluid specific heat at constant pressure
                              in units of J/kg-K
        :param molar_mass: The fluid molar mass in units of J/mol-K
        :param inlet_mach_number: The inlet Mach number
        :param mass_flow_rate: The fluid mass flow in units of kg/s
        :param inlet_stagnation_pressure: The stagnation pressure at the
                                          inlet to the compressor in units
                                          of Pascals
        :param inlet_stagnation_temperature: The stagnation temperature at the
                                             inlet to the compressor in units of
                                             Kelvins
        :return dict: A dictionary containing all fluid exit
                      properties

        The dictionary representing compressor contains
        the following key words

        * **static_temperature** in units of Kelvins
        * **static_pressure** in units of Pascals
        * **stagnation_temperature** in units of Kelvins
        * **stagnation_pressure** in units of Pascals
        * **velocity** in units of meters per second
        * **mach_number**
        * **work** in units of Watts
        """
        exit_stag_pres = self.exit_stagnation_pressure(inlet_stagnation_pressure)
        exit_stag_temp = self.exit_stagnation_temperature(inlet_stagnation_temperature,
                                                          gamma)
        comp_work = self.compressor_work(mass_flow_rate, specific_heat,
                                         inlet_stagnation_temperature, gamma)
        exit_mach_number = self.exit_mach_number(inlet_mach_number,
                                                 inlet_stagnation_temperature,
                                                 gamma)
        exit_stat_pres = self.exit_static_pressure(inlet_stagnation_pressure,
                                                   inlet_stagnation_temperature,
                                                   gamma, inlet_mach_number)
        exit_stat_temp = self.exit_static_temperature(inlet_stagnation_temperature,
                                                      gamma, inlet_mach_number)
        exit_velocity = self.exit_velocity(inlet_stagnation_temperature, gamma,
                                           inlet_mach_number, molar_mass)
        dictionary = {'static_pressure': exit_stat_pres, 'static_temperature': exit_stat_temp,
                      'stagnation_pressure': exit_stag_pres, 'stagnation_temperature': exit_stag_temp,
                      'velocity': exit_velocity, 'mach_number': exit_mach_number,
                      'work': comp_work}
        return dictionary
# ============================================================================
# ============================================================================


class HeatAdditionComponent(HeatAddition):
    """

    This class combined all functionality of the Compressor class
    into a single function for ease of user access
    """
    def __init__(self, efficiency: float):
        """

        :param efficiency: The isentropic efficiency of the heat
                           addition component
        """
        HeatAddition.__init__(self, efficiency)
# ----------------------------------------------------------------------------

    def outlet_conditions(self, gamma: float, specific_heat: float,
                          molar_mass: float,
                          inlet_mach_number: float, mass_flow_rate: float,
                          inlet_stagnation_pressure: float,
                          inlet_stagnation_temperature: float,
                          heat_addition: float) -> Dict[str, float]:
        """

        :param gamma: The ratio of specific heats
        :param specific_heat: The fluid specific heat at constant pressure
                              in units of J/kg-K
        :param molar_mass: The fluid molar mass in units of J/mol-K
        :param inlet_mach_number: The inlet Mach number
        :param mass_flow_rate: The fluid mass flow in units of kg/s
        :param inlet_stagnation_pressure: The stagnation pressure at the
                                          inlet to the compressor in units
                                          of Pascals
        :param inlet_stagnation_temperature: The stagnation temperature at the
                                             inlet to the compressor in units of
                                             Kelvins
        :param heat_addition: The heat added to the fluid in units of
                              Watts
        :return dict: A dictionary containing all fluid exit
                      properties

        The dictionary containing the heat exchanger parameters
        contains the following key words

        * **static_temperature** in units of Kelvins
        * **static_pressure** in units of Pascals
        * **stagnation_temperature** in units of Kelvins
        * **stagnation_pressure** in units of Pascals
        * **velocity** in units of meters per second
        * **mach_number**
        * **power** in units of Watts
        """
        input_power = self.input_power(heat_addition)
        exit_stag_temp = self.exit_stagnation_temperature(inlet_stagnation_temperature,
                                                          heat_addition, specific_heat,
                                                          mass_flow_rate)
        exit_mach_number = self.exit_mach_number(inlet_stagnation_temperature,
                                                 heat_addition, specific_heat,
                                                 mass_flow_rate, inlet_mach_number,
                                                 gamma)
        exit_stag_pres = self.exit_stagnation_pressure(inlet_stagnation_pressure,
                                                       inlet_stagnation_temperature,
                                                       heat_addition, specific_heat,
                                                       mass_flow_rate, inlet_mach_number,
                                                       gamma)
        exit_stat_temp = self.exit_static_temperature(inlet_stagnation_temperature,
                                                      heat_addition, specific_heat,
                                                      mass_flow_rate, gamma,
                                                      inlet_mach_number)
        exit_stat_pres = self.exit_static_pressure(inlet_stagnation_pressure,
                                                   inlet_stagnation_temperature,
                                                   heat_addition, specific_heat,
                                                   mass_flow_rate, gamma,
                                                   inlet_mach_number)
        exit_velocity = self.exit_velocity(inlet_stagnation_temperature,
                                           heat_addition, specific_heat,
                                           mass_flow_rate, inlet_mach_number,
                                           gamma, molar_mass)
        dictionary = {'static_pressure': exit_stat_pres, 'static_temperature': exit_stat_temp,
                      'stagnation_pressure': exit_stag_pres, 'stagnation_temperature': exit_stag_temp,
                      'velocity': exit_velocity, 'mach_number': exit_mach_number,
                      'power': input_power}
        return dictionary
# ============================================================================
# ============================================================================


class TurbineComponent(Turbine):
    """

        This class combined all functionality of the Turbine class
        into a single function for ease of user access.  This class can
        also be used to model a heat sink.
        """

    def __init__(self, efficiency: float):
        """

        :param efficiency: The isentropic efficiency of the turbine
        """
        Turbine.__init__(self, efficiency)
# ----------------------------------------------------------------------------

    def outlet_conditions(self, gamma: float, specific_heat: float,
                          molar_mass: float,
                          inlet_mach_number: float, mass_flow_rate: float,
                          inlet_stagnation_pressure: float,
                          inlet_stagnation_temperature: float,
                          turbine_work: float) -> Dict[str, float]:
        """

        :param gamma: The ratio of specific specific heats
        :param specific_heat: The fluid specific heat at constant pressure
                              in units of J/kg-K
        :param molar_mass: The fluid molar mass in units of J/mol-K
        :param inlet_mach_number: The Mach number at the turbine inlet
        :param mass_flow_rate: The fluid mass flow rate in units of kg/s
        :param inlet_stagnation_pressure: The stagnation pressure at the turbine
                                          inlet in units of Pascals
        :param inlet_stagnation_temperature: The stagnation temperature at the
                                             turbine inlet in units of Pascals
        :param turbine_work: The usable work extracted by the turbine in
                             units of Watts
        :return dict: A dictionary containing all fluid exit
                      properties

        The dictionary representing the turbine parameters contains the
        following key words

        * **static_temperature** in units of Kelvins
        * **static_pressure** in units of Pascals
        * **stagnation_temperature** in units of Kelvins
        * **stagnation_pressure** in units of Pascals
        * **velocity** in units of meters per second
        * **mach_number**
        * **extracted_work** in units of Watts
        """
        work_extracted = self.work_extraction(turbine_work)
        exit_stag_temp = self.exit_stagnation_temperature(turbine_work, specific_heat,
                                                          mass_flow_rate,
                                                          inlet_stagnation_temperature)
        exit_stag_pres = self.exit_stagnation_pressure(inlet_stagnation_pressure,
                                                       inlet_stagnation_temperature,
                                                       specific_heat, mass_flow_rate,
                                                       turbine_work, gamma)
        exit_mach_number = self.exit_mach_number(inlet_stagnation_temperature,
                                                 inlet_mach_number, gamma,
                                                 turbine_work, specific_heat,
                                                 mass_flow_rate)
        exit_stat_temp = self.exit_static_temperature(inlet_stagnation_temperature,
                                                      inlet_mach_number, gamma,
                                                      turbine_work, specific_heat,
                                                      mass_flow_rate)
        exit_stat_pres = self.exit_static_pressure(inlet_stagnation_temperature,
                                                   inlet_mach_number, gamma,
                                                   turbine_work, specific_heat,
                                                   mass_flow_rate,
                                                   inlet_stagnation_pressure)
        exit_velocity = self.exit_velocity(inlet_stagnation_temperature,
                                           inlet_mach_number, gamma,
                                           turbine_work, specific_heat,
                                           mass_flow_rate, molar_mass)
        dictionary = {'static_pressure': exit_stat_pres, 'static_temperature': exit_stat_temp,
                      'stagnation_pressure': exit_stag_pres, 'stagnation_temperature': exit_stag_temp,
                      'velocity': exit_velocity, 'mach_number': exit_mach_number,
                      'extracted_work': work_extracted}
        return dictionary
# ============================================================================
# ============================================================================


class PropellerComponent(Propeller):
    """

        This class combined all functionality of the Propeller class
        into a single function for ease of user access
        """

    def __init__(self, efficiency: float):
        """

        :param efficiency: The isentropic efficiency of the propeller
        """
        Propeller.__init__(self, efficiency)
# ----------------------------------------------------------------------------

    def outlet_conditions(self, extracted_work: float,
                          inlet_stagnation_temperature: float,
                          mass_flow_rate: float, specific_heat: float,
                          inlet_stagnation_pressure: float, gamma: float,
                          inlet_mach_number: float, molar_mass: float) -> Dict[str, float]:
        """

        :param extracted_work: The total usable work done by the propeller
        :param inlet_stagnation_temperature: The stagnation temperature at
                                             the propeller inlet in units
                                             of Kelvins
        :param mass_flow_rate: The fluid mass flow rate in units of
                               kg/s
        :param specific_heat: The fluid specific heat in units of J/kg-K
        :param inlet_stagnation_pressure: The stagnation pressure at the
                                          propeller inlet in units of Pascals
        :param gamma: The ratio of specific heats
        :param inlet_mach_number: The Mach number at the propeller inlet
        :param molar_mass: The fluid molar mass
        :return dict: A dictionary containing all fluid exit properties

        The dictionary representing the propeller exit conditions
        contains the following key words

        * **static_temperature** in units of Kelvins
        * **static_pressure** in units of Pascals
        * **stagnation_temperature** in units of Kelvins
        * **stagnation_pressure** in units of Pascals
        * **velocity** in units of meters per second
        * **mach_number**
        * **total_work** The total work done by the propeller in units of Watts
        """
        total_work = self.propeller_work(extracted_work)
        exit_stag_temp = self.exit_stagnation_temperature(inlet_stagnation_temperature,
                                                          mass_flow_rate, specific_heat,
                                                          extracted_work)
        exit_stag_pres = self.exit_stagnation_pressure(inlet_stagnation_pressure,
                                                       inlet_stagnation_temperature,
                                                       gamma, mass_flow_rate,
                                                       specific_heat, extracted_work)
        exit_mach_number = self.exit_mach_number(inlet_stagnation_temperature, inlet_mach_number,
                                                 mass_flow_rate, specific_heat, extracted_work,
                                                 gamma)
        exit_stat_temp = self.exit_static_temperature(inlet_stagnation_temperature,
                                                      inlet_mach_number, mass_flow_rate,
                                                      specific_heat, extracted_work,
                                                      gamma)
        exit_stat_pres = self.exit_static_pressure(inlet_stagnation_temperature,
                                                   inlet_mach_number, mass_flow_rate,
                                                   specific_heat, extracted_work,
                                                   gamma, inlet_stagnation_pressure)
        exit_velocity = self.exit_velocity(inlet_stagnation_temperature,
                                           inlet_mach_number, mass_flow_rate,
                                           specific_heat, extracted_work,
                                           gamma, molar_mass)
        dictionary = {'static_pressure': exit_stat_pres, 'static_temperature': exit_stat_temp,
                      'stagnation_pressure': exit_stag_pres, 'stagnation_temperature': exit_stag_temp,
                      'velocity': exit_velocity, 'mach_number': exit_mach_number,
                      'total_work': total_work}
        return dictionary
# ============================================================================
# ============================================================================


class RamJet:
    def __init__(self, dif_eff: float, dif_inlet_area: float, dif_exit_area: float,
                 heat_eff: float, noz_eff: float, noz_inlet_area: float,
                 noz_exit_area: float, species: str):
        """

        :param dif_eff: The diffuser isentropic efficiency
        :param dif_inlet_area: The diffuser inlet area in square meters
        :param dif_exit_area: The diffuser exit area in square meters
        :param heat_eff: The head addition isentropic efficiency
        :param noz_eff: The nozzle isentropic efficiency
        :param noz_inlet_area: The nozzle inlet area in square meters
        :param noz_exit_area: The nozzle exit area in square meters
        :param species: The gas species
        """
        self.diffuser = DiffuserNozzleComponent(dif_eff, dif_inlet_area, dif_exit_area)
        self.heat = HeatAdditionComponent(heat_eff)
        self.nozzle = DiffuserNozzleComponent(noz_eff, noz_inlet_area, noz_exit_area)
        self.gas = ThermProps(species)
# ----------------------------------------------------------------------------

    def performance(self, inlet_static_temp: float, inlet_static_pressure: float,
                    inlet_mach_number: float, inlet_velocity: float,
                    mass_flow_rate: float, heat_addition: float) -> Tuple[Dict[str, float],
                                                                          Dict[str, float],
                                                                          Dict[str, float]]:
        """

        :param inlet_static_temp: The inlet static temperature in units of Kelvins
        :param inlet_static_pressure: The inlet static pressure in units of Pascals
        :param inlet_mach_number: The inlet mach number
        :param inlet_velocity: The inlet flow velocity
        :param mass_flow_rate: The inlet mass flow rate
        :param heat_addition: The work added to the heat addition chamber in units
                              of Watts
        :return dictionary: A Tuple of three dictionaries containing the
                            gas properties at the exit to each component of
                            the ram jet

        The first dictionary for the diffuser contains the following key words

        * **static_temperature** in units of Kelvins
        * **static_pressure** in units of Pascals
        * **stagnation_temperature** in units of Kelvins
        * **stagnation_pressure** in units of Pascals
        * **velocity** in units of meters per second
        * **mach_number**
        * **density** in units of g/cc

        The second dictionary represents heat exchanger parameters and contains
        the following key words

        * **static_temperature** in units of Kelvins
        * **static_pressure** in units of Pascals
        * **stagnation_temperature** in units of Kelvins
        * **stagnation_pressure** in units of Pascals
        * **velocity** in units of meters per second
        * **mach_number**
        * **power** in units of Watts

        The third dictionary represents nozzle parameters and contains the
        followingf key words

        * **static_temperature** in units of Kelvins
        * **static_pressure** in units of Pascals
        * **stagnation_temperature** in units of Kelvins
        * **stagnation_pressure** in units of Pascals
        * **velocity** in units of meters per second
        * **mach_number**
        * **density** in units of g/cc
        """
        cp = self.gas.spec_heat_const_pressure(inlet_static_temp, inlet_static_pressure)
        gamma = self.gas.ratio_specific_heats(inlet_static_temp, inlet_static_pressure)
        mw = self.gas.molar_mass(inlet_static_temp, inlet_static_pressure)

        # Calculate diffuser properties
        dif_prop = self.diffuser.outlet_conditions(gamma, cp, mw,
                                                   inlet_static_temp,
                                                   inlet_static_pressure,
                                                   inlet_mach_number, mass_flow_rate,
                                                   inlet_velocity)
        cp = self.gas.spec_heat_const_pressure(dif_prop['static_temperature'],
                                               dif_prop['static_pressure'])
        gamma = self.gas.ratio_specific_heats(dif_prop['static_temperature'],
                                              dif_prop['static_pressure'])
        mw = self.gas.molar_mass(dif_prop['static_temperature'],
                                 dif_prop['static_pressure'])
        # Calculate heat addition properties
        he_props = self.heat.outlet_conditions(gamma, cp, mw, dif_prop['mach_number'],
                                               mass_flow_rate, dif_prop['stagnation_pressure'],
                                               dif_prop['stagnation_temperature'],
                                               heat_addition)
        cp = self.gas.spec_heat_const_pressure(he_props['static_temperature'],
                                               he_props['static_pressure'])
        gamma = self.gas.ratio_specific_heats(he_props['static_temperature'],
                                              he_props['static_pressure'])
        mw = self.gas.molar_mass(he_props['static_temperature'],
                                 he_props['static_pressure'])
        # Calculate nozzle properties
        noz_props = self.nozzle.outlet_conditions(gamma, cp, mw, he_props['static_temperature'],
                                                  he_props['static_pressure'],
                                                  he_props['mach_number'],
                                                  mass_flow_rate, he_props['velocity'])
        return dif_prop, he_props, noz_props
# ============================================================================
# ============================================================================


class TurboJet:
    def __init__(self, dif_eff: float, dif_inlet_area: float, dif_exit_area: float,
                 comp_ratio: float, comp_eff: float, heat_eff: float, turb_eff: float,
                 noz_eff: float, noz_inlet_area: float, noz_exit_area: float,
                 species: str):
        """

        :param dif_eff: The diffuser isentropic efficiency
        :param dif_inlet_area: The diffuser inlet area in units of square meters
        :param dif_exit_area: The diffuser exit area in units of square meters
        :param comp_ratio: The compressor ratio
        :param comp_eff: The compressor isentropic efficiency
        :param heat_eff: The heat addition component isentropic efficiency
        :param turb_eff: The turbine isentropic efficiency
        :param noz_eff: The nozzle isentropic efficiency
        :param noz_inlet_area: The nozzle inlet area in units of square meters
        :param noz_exit_area: The nozzle exit area in units of square meters
        :param species: The gas species
        """
        self.diffuser = DiffuserNozzleComponent(dif_eff, dif_inlet_area, dif_exit_area)
        self.compressor = CompressorComponent(comp_ratio, comp_eff)
        self.heat = HeatAdditionComponent(heat_eff)
        self.turbine = TurbineComponent(turb_eff)
        self.nozzle = DiffuserNozzleComponent(noz_eff, noz_inlet_area, noz_exit_area)
        self.gas = ThermProps(species)
# ----------------------------------------------------------------------------

    def performance(self, inlet_static_temp: float, inlet_static_pressure: float,
                    inlet_mach_number: float, mass_flow_rate: float,
                    inlet_velocity: float, heat_addition: float,
                    turbine_work: float) -> Tuple[Dict[str, float],
                                                  Dict[str, float],
                                                  Dict[str, float],
                                                  Dict[str, float],
                                                  Dict[str, float]]:
        """
        
        :param inlet_static_temp: The inlet static temperature in units of Kelvins
        :param inlet_static_pressure: The inlet static temperature in units of Pascals
        :param inlet_mach_number: The inlet Mach number
        :param mass_flow_rate: The mass flow rate in units of kg/s
        :param inlet_velocity: The inlet velocity in units of meters per second
        :param heat_addition: The work added to the heat addition chamber in units
                              of Watts
        :return dict: A tuple containing five dictionaries

        The first dictionary for the diffuser contains the following key words

        * **static_temperature** in units of Kelvins
        * **static_pressure** in units of Pascals
        * **stagnation_temperature** in units of Kelvins
        * **stagnation_pressure** in units of Pascals
        * **velocity** in units of meters per second
        * **mach_number**
        * **density** in units of g/cc

        The second dictionary represents compressor parameters and contains
        the following key words

        * **static_temperature** in units of Kelvins
        * **static_pressure** in units of Pascals
        * **stagnation_temperature** in units of Kelvins
        * **stagnation_pressure** in units of Pascals
        * **velocity** in units of meters per second
        * **mach_number**
        * **work** in units of Watts

        The third dictionary contains heat exchanger parameters and
        contains the following key words

        * **static_temperature** in units of Kelvins
        * **static_pressure** in units of Pascals
        * **stagnation_temperature** in units of Kelvins
        * **stagnation_pressure** in units of Pascals
        * **velocity** in units of meters per second
        * **mach_number**
        * **power** in units of Watts

        The fourth dictionary represents turbine parameters and contains the
        following key words

        * **static_temperature** in units of Kelvins
        * **static_pressure** in units of Pascals
        * **stagnation_temperature** in units of Kelvins
        * **stagnation_pressure** in units of Pascals
        * **velocity** in units of meters per second
        * **mach_number**
        * **extracted_work** in units of Watts

        The final dictionary for the nozzle contains the following
        key words

        * **static_temperature** in units of Kelvins
        * **static_pressure** in units of Pascals
        * **stagnation_temperature** in units of Kelvins
        * **stagnation_pressure** in units of Pascals
        * **velocity** in units of meters per second
        * **mach_number**
        * **density** in units of g/cc
        """

        cp = self.gas.spec_heat_const_pressure(inlet_static_temp, inlet_static_pressure)
        gamma = self.gas.ratio_specific_heats(inlet_static_temp, inlet_static_pressure)
        mw = self.gas.molar_mass(inlet_static_temp, inlet_static_pressure)

        # Calculate diffuser properties
        dif_prop = self.diffuser.outlet_conditions(gamma, cp, mw,
                                                   inlet_static_temp,
                                                   inlet_static_pressure,
                                                   inlet_mach_number, mass_flow_rate,
                                                   inlet_velocity)
        # Calculate compressor properties
        cp = self.gas.spec_heat_const_pressure(dif_prop['static_temperature'],
                                               dif_prop['static_pressure'])
        gamma = self.gas.ratio_specific_heats(dif_prop['static_temperature'],
                                              dif_prop['static_pressure'])
        mw = self.gas.molar_mass(dif_prop['static_temperature'],
                                 dif_prop['static_pressure'])
        comp_prop = self.compressor.outlet_conditions(gamma, cp, mw,
                                                      dif_prop['mach_number'],
                                                      mass_flow_rate,
                                                      dif_prop['stagnation_pressure'],
                                                      dif_prop['stagnation_temperature'])
        # Calculate heat addition properties
        cp = self.gas.spec_heat_const_pressure(comp_prop['static_temperature'],
                                               comp_prop['static_pressure'])
        gamma = self.gas.ratio_specific_heats(comp_prop['static_temperature'],
                                              comp_prop['static_pressure'])
        mw = self.gas.molar_mass(comp_prop['static_temperature'],
                                 comp_prop['static_pressure'])

        he_prop = self.heat.outlet_conditions(gamma, cp, mw, comp_prop['mach_number'],
                                              mass_flow_rate, comp_prop['stagnation_pressure'],
                                              comp_prop['stagnation_temperature'],
                                              heat_addition)
        # Calculate turbine properties
        cp = self.gas.spec_heat_const_pressure(he_prop['static_temperature'],
                                               he_prop['static_pressure'])
        gamma = self.gas.ratio_specific_heats(he_prop['static_temperature'],
                                              he_prop['static_pressure'])
        mw = self.gas.molar_mass(he_prop['static_temperature'],
                                 he_prop['static_pressure'])
        turb_prop = self.turbine.outlet_conditions(gamma, cp, mw,
                                                   he_prop['mach_number'],
                                                   mass_flow_rate,
                                                   he_prop['stagnation_pressure'],
                                                   he_prop['stagnation_temperature'],
                                                   turbine_work)
        # Calculate nozzle properties
        cp = self.gas.spec_heat_const_pressure(turb_prop['static_temperature'],
                                               turb_prop['static_pressure'])
        gamma = self.gas.ratio_specific_heats(turb_prop['static_temperature'],
                                              turb_prop['static_pressure'])
        mw = self.gas.molar_mass(turb_prop['static_temperature'],
                                 turb_prop['static_pressure'])
        noz_prop = self.nozzle.outlet_conditions(gamma, cp, mw,
                                                 turb_prop['static_temperature'],
                                                 turb_prop['static_pressure'],
                                                 turb_prop['mach_number'],
                                                 mass_flow_rate,
                                                 turb_prop['velocity'])
        return dif_prop, comp_prop, he_prop, turb_prop, noz_prop
# ============================================================================
# ============================================================================


class TurboProp:
    def __init__(self, dif_eff: float, dif_inlet_area: float, dif_exit_area: float,
                 comp_ratio: float, comp_eff: float, heat_eff: float, turb_eff: float,
                 noz_eff: float, noz_inlet_area: float, noz_exit_area: float,
                 prop_eff: float, species: str):
        """

        :param dif_eff: The diffuser isentropic efficiency
        :param dif_inlet_area: The diffuser inlet area in units of square meters
        :param dif_exit_area: The diffuser exit area in units of square meters
        :param comp_ratio: The compressor ratio
        :param comp_eff: The compressor isentropic efficiency
        :param heat_eff: The heat addition component isentropic efficiency
        :param turb_eff: The turbine isentropic efficiency
        :param noz_eff: The nozzle isentropic efficiency
        :param noz_inlet_area: The nozzle inlet area in units of square meters
        :param noz_exit_area: The nozzle exit area in units of square meters
        :param species: The gas species
        """
        self.propeller = PropellerComponent(prop_eff)
        self.diffuser = DiffuserNozzleComponent(dif_eff, dif_inlet_area, dif_exit_area)
        self.compressor = CompressorComponent(comp_ratio, comp_eff)
        self.heat = HeatAdditionComponent(heat_eff)
        self.turbine = TurbineComponent(turb_eff)
        self.nozzle = DiffuserNozzleComponent(noz_eff, noz_inlet_area, noz_exit_area)
        self.gas = ThermProps(species)
# ----------------------------------------------------------------------------

    def performance(self, work: float, inlet_stagnation_temperature: float,
                    mass_flow_rate: float, inlet_stagnation_pressure: float,
                    inlet_mach_number: float, inlet_static_temp: float,
                    inlet_static_pressure: float, heat_addition: float,
                    turbine_work: float) -> Tuple[Dict[str, float],
                                                  Dict[str, float],
                                                  Dict[str, float],
                                                  Dict[str, float],
                                                  Dict[str, float],
                                                  Dict[str, float]]:
        """

        :param work: The work done by the propeller in units of Watts
        :param inlet_stagnation_temperature: The stagnation temperature
                                             at the propeller inlet in
                                             units of Kelvins
        :param mass_flow_rate: The mass flow rate in units of kg/s
        :param inlet_stagnation_pressure: The stagnation pressure in units
                                          of Pascals
        :param inlet_mach_number: The inlet mach number
        :param inlet_static_temp: The static temperature in units of Kelvins
        :param inlet_static_pressure: The inlet static pressure in units of
                                      Pascals
        :param heat_addition: The heat added in the head addition chamber in
                              units of Watts
        turb_work: The work extracted by the turbine in units of Watts
        :return: A tuple containing six dictionaries

        The first dictionary representing the propeller exit conditions
        contains the following key words

        * **static_temperature** in units of Kelvins
        * **static_pressure** in units of Pascals
        * **stagnation_temperature** in units of Kelvins
        * **stagnation_pressure** in units of Pascals
        * **velocity** in units of meters per second
        * **mach_number**
        * **total_work** The total work done by the propeller in units of Watts

        The second dictionary for the diffuser contains the following key words

        * **static_temperature** in units of Kelvins
        * **static_pressure** in units of Pascals
        * **stagnation_temperature** in units of Kelvins
        * **stagnation_pressure** in units of Pascals
        * **velocity** in units of meters per second
        * **mach_number**
        * **density** in units of g/cc

        The third dictionary represents compressor parameters and contains
        the following key words

        * **static_temperature** in units of Kelvins
        * **static_pressure** in units of Pascals
        * **stagnation_temperature** in units of Kelvins
        * **stagnation_pressure** in units of Pascals
        * **velocity** in units of meters per second
        * **mach_number**
        * **work** in units of Watts

        The fourth dictionary contains heat exchanger parameters and
        contains the following key words

        * **static_temperature** in units of Kelvins
        * **static_pressure** in units of Pascals
        * **stagnation_temperature** in units of Kelvins
        * **stagnation_pressure** in units of Pascals
        * **velocity** in units of meters per second
        * **mach_number**
        * **power** in units of Watts

        The fifth dictionary represents turbine parameters and contains the
        following key words

        * **static_temperature** in units of Kelvins
        * **static_pressure** in units of Pascals
        * **stagnation_temperature** in units of Kelvins
        * **stagnation_pressure** in units of Pascals
        * **velocity** in units of meters per second
        * **mach_number**
        * **extracted_work** in units of Watts

        The final dictionary for the nozzle contains the following
        key words

        * **static_temperature** in units of Kelvins
        * **static_pressure** in units of Pascals
        * **stagnation_temperature** in units of Kelvins
        * **stagnation_pressure** in units of Pascals
        * **velocity** in units of meters per second
        * **mach_number**
        * **density** in units of g/cc
        """
        # Calculate propeller outlet conditions
        cp = self.gas.spec_heat_const_pressure(inlet_static_temp, inlet_static_pressure)
        gamma = self.gas.ratio_specific_heats(inlet_static_temp, inlet_static_pressure)
        mw = self.gas.molar_mass(inlet_static_temp, inlet_static_pressure)
        prop_prop = self.propeller.outlet_conditions(work, inlet_stagnation_temperature,
                                                     mass_flow_rate, cp,
                                                     inlet_stagnation_pressure, gamma,
                                                     inlet_mach_number, mw)
        # Calculate Diffuser properties
        cp = self.gas.spec_heat_const_pressure(prop_prop['static_temperature'],
                                               prop_prop['static_pressure'])
        gamma = self.gas.ratio_specific_heats(prop_prop['static_temperature'],
                                              prop_prop['static_pressure'])
        mw = self.gas.molar_mass(prop_prop['static_temperature'],
                                 prop_prop['static_pressure'])
        diff_prop = self.diffuser.outlet_conditions(gamma, cp, mw,
                                                    prop_prop['static_temperature'],
                                                    prop_prop['static_pressure'],
                                                    prop_prop['mach_number'],
                                                    mass_flow_rate,
                                                    prop_prop['velocity'])
        # Calculate Compressor properties
        cp = self.gas.spec_heat_const_pressure(diff_prop['static_temperature'],
                                               diff_prop['static_pressure'])
        gamma = self.gas.ratio_specific_heats(diff_prop['static_temperature'],
                                              diff_prop['static_pressure'])
        mw = self.gas.molar_mass(diff_prop['static_temperature'],
                                 diff_prop['static_pressure'])
        comp_prop = self.compressor.outlet_conditions(gamma, cp, mw,
                                                      diff_prop['mach_number'],
                                                      mass_flow_rate,
                                                      diff_prop['stagnation_pressure'],
                                                      diff_prop['stagnation_temperature'])
        # Calculate heat addition properties
        cp = self.gas.spec_heat_const_pressure(comp_prop['static_temperature'],
                                               comp_prop['static_pressure'])
        gamma = self.gas.ratio_specific_heats(comp_prop['static_temperature'],
                                              comp_prop['static_pressure'])
        mw = self.gas.molar_mass(comp_prop['static_temperature'],
                                 comp_prop['static_pressure'])

        he_prop = self.heat.outlet_conditions(gamma, cp, mw, comp_prop['mach_number'],
                                              mass_flow_rate, comp_prop['stagnation_pressure'],
                                              comp_prop['stagnation_temperature'],
                                              heat_addition)
        # Calculate turbine properties
        cp = self.gas.spec_heat_const_pressure(he_prop['static_temperature'],
                                               he_prop['static_pressure'])
        gamma = self.gas.ratio_specific_heats(he_prop['static_temperature'],
                                              he_prop['static_pressure'])
        mw = self.gas.molar_mass(he_prop['static_temperature'],
                                 he_prop['static_pressure'])
        turb_prop = self.turbine.outlet_conditions(gamma, cp, mw,
                                                   he_prop['mach_number'],
                                                   mass_flow_rate,
                                                   he_prop['stagnation_pressure'],
                                                   he_prop['stagnation_temperature'],
                                                   turbine_work)
        # Calculate nozzle properties
        cp = self.gas.spec_heat_const_pressure(turb_prop['static_temperature'],
                                               turb_prop['static_pressure'])
        gamma = self.gas.ratio_specific_heats(turb_prop['static_temperature'],
                                              turb_prop['static_pressure'])
        mw = self.gas.molar_mass(turb_prop['static_temperature'],
                                 turb_prop['static_pressure'])
        noz_prop = self.nozzle.outlet_conditions(gamma, cp, mw,
                                                 turb_prop['static_temperature'],
                                                 turb_prop['static_pressure'],
                                                 turb_prop['mach_number'],
                                                 mass_flow_rate,
                                                 turb_prop['velocity'])
        return prop_prop, diff_prop, comp_prop, he_prop, turb_prop, noz_prop
# ============================================================================
# ============================================================================


class ReactorLoop:
    def __init__(self, heat_eff: float, turb_eff: float, species: str,
                 sink_eff: float = 1):
        """

        :param heat_eff: The heat addition efficiency
        :param turb_eff: The turbine efficiency
        :param sink_eff: The heat sink efficiency
        :param species: The gas species
        """
        self.heat = HeatAdditionComponent(heat_eff)
        self.turbine = TurbineComponent(turb_eff)
        self.sink = Turbine(sink_eff)
        self.gas = ThermProps(species)
        self.turb_eff = turb_eff
        self.heat_eff = heat_eff
# ----------------------------------------------------------------------------

    def performance(self, inlet_stag_temp: float, inlet_stag_pres: float,
                   inlet_mach_number: float, inlet_stat_temp: float,
                   inlet_stat_pres: float, elec_power: float,
                   comp_ratio_guess: float, comp_eff: float,
                   conv_eff: float, mass_flow_rate) -> Tuple[Dict[str, float],
                                                             Dict[str, float],
                                                             Dict[str, float],
                                                             Dict[str, float]]:

        # Calculate compressor fluid properties
        cp = self.gas.spec_heat_const_pressure(inlet_stat_temp, inlet_stat_pres)
        gamma = self.gas.ratio_specific_heats(inlet_stat_temp, inlet_stat_pres)
        mw = self.gas.molar_mass(inlet_stat_temp, inlet_stat_pres)
        compressor = CompressorComponent(comp_ratio_guess, comp_eff)
        comp_prop = compressor.outlet_conditions(gamma, cp, mw, inlet_mach_number,
                                                 mass_flow_rate, inlet_stag_pres,
                                                 inlet_stag_temp)

        extract_power = (elec_power + comp_prop['work']) / self.heat_eff
        total_power = extract_power / conv_eff
        # Calculate reactor fluid properties
        cp = self.gas.spec_heat_const_pressure(comp_prop['static_temperature'],
                                               comp_prop['static_pressure'])
        gamma = self.gas.ratio_specific_heats(comp_prop['static_temperature'],
                                              comp_prop['static_pressure'])
        mw = self.gas.molar_mass(comp_prop['static_temperature'],
                                 comp_prop['static_pressure'])
        he_prop = self.heat.outlet_conditions(gamma, cp, mw, comp_prop['mach_number'],
                                              mass_flow_rate, comp_prop['stagnation_pressure'],
                                              comp_prop['stagnation_temperature'],
                                              total_power)

        # Calculate turbine fluid properties
        cp = self.gas.spec_heat_const_pressure(he_prop['static_temperature'],
                                               he_prop['static_pressure'])
        gamma = self.gas.ratio_specific_heats(he_prop['static_temperature'],
                                              he_prop['static_pressure'])
        mw = self.gas.molar_mass(he_prop['static_temperature'],
                                 he_prop['static_pressure'])
        turb_prop = self.turbine.outlet_conditions(gamma, cp, mw,
                                                   he_prop['mach_number'],
                                                   mass_flow_rate,
                                                   he_prop['stagnation_pressure'],
                                                   he_prop['stagnation_temperature'],
                                                   extract_power)
        # Calculate turbine fluid properties
        cp = self.gas.spec_heat_const_pressure(turb_prop['static_temperature'],
                                               turb_prop['static_pressure'])
        gamma = self.gas.ratio_specific_heats(turb_prop['static_temperature'],
                                              turb_prop['static_pressure'])
        mw = self.gas.molar_mass(turb_prop['static_temperature'],
                                 turb_prop['static_pressure'])

        # calculate rejected power
        reject = he_prop['power'] - turb_prop['extracted_work']
        reject_prop = self.turbine.outlet_conditions(gamma, cp, mw,
                                                     turb_prop['mach_number'],
                                                     mass_flow_rate,
                                                     turb_prop['stagnation_pressure'],
                                                     turb_prop['stagnation_temperature'],
                                                     reject)
        return comp_prop, he_prop, turb_prop, reject_prop
# ============================================================================
# ============================================================================
# eof
