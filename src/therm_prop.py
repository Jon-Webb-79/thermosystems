# Import necessary packages here
from thermo.chemical import Chemical
# ============================================================================
# ============================================================================
# Date:    October 3, 2020
# Purpose: This file contains classes describing the thermodynamic
#          properties of gases

# Source Code Metadata
__author__ = "Jonathan A. Webb"
__copyright__ = "Copyright 2020, Jon Webb Inc."
__version__ = "1.0"

# ============================================================================
# ============================================================================git


class ThermProps:
    def __init__(self, species: str):
        """

        :param species: The species of gas to be modeled
        """
        self.chem = Chemical(species)
# ----------------------------------------------------------------------------

    def spec_heat_const_pressure(self, temperature: float, pressure: float) -> float:
        """

        :param temperature: The static temperature in units of Kelvins
        :param pressure: The static pressure in units of Pascals
        :return cp: The specific heat at constant pressure in units
                    of J/kg-k
        """
        self.chem.calculate(T=temperature, P=pressure)
        return self.chem.Cp
# ----------------------------------------------------------------------------

    def spec_heat_const_volume(self, temperature: float, pressure: float) -> float:
        """

        :param temperature: The static temperature in units of Kelvins
        :param pressure: The static pressure in units of Pascals
        :return cv: The specific heat at constant volume in units
                    of J/kg-k
        """
        self.chem.calculate(T=temperature, P=pressure)
        return self.chem.Cvg
# ----------------------------------------------------------------------------

    def ratio_specific_heats(self, temperature: float, pressure: float) -> float:
        """

        :param temperature: The static temperature in units of Kelvins
        :param pressure: The static pressure in units of Pascals
        :return gamma: The ratio of specific heats
        """
        cp = self.spec_heat_const_pressure(temperature, pressure)
        cv = self.spec_heat_const_volume(temperature, pressure)
        return cp / cv
# ----------------------------------------------------------------------------

    def molar_mass(self, temperature: float, pressure: float) -> float:
        """

        :param temperature: The static temperature in units of Kelvins
        :param pressure: The static pressure in units of Pascals
        :return mw: THe molecular weight in units of g/mol
        """
        self.chem.calculate(T=temperature, P=pressure)
        return self.chem.MW
# ============================================================================
# ============================================================================
# eof
