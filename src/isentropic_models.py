# Import necessary packages here

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
        :param molar_mass: Molar mass of the fluid in units of J/mol-K
        :return sos: The speed of sound in meters per second

        This function determines the speed of sound
        using the relationship from pg. 68 or reference 1,
        shown below

        .. math::
           a = \sqrt[]{\gamma R T}
        """
        universal_gas_constant = 8.314  # J/K-mol
        gas_constant = universal_gas_constant / molar_mass
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
        return stag_temp - ((velocity ** 2) / 2 * specific_heat)
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
# ============================================================================
# ============================================================================
# eof
