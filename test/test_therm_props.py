# Import modules here
from src.therm_prop import ThermProps

from math import isclose
# ==============================================================================
# ==============================================================================
# Date:    October 13, 2020
# Purpose: This code contains functions that test the functions in the
#          ThermProps class

# Source Code Metadata
__author__ = "Jonathan A. Webb"
__copyright__ = "Copyright 2020, Jon Webb Inc."
__version__ = "1.0"
# ==============================================================================
# ==============================================================================
# Test the ThermProps with nitrogen


species = 'nitrogen'
gas = ThermProps(species)
temp = 340.0
pres = 8000000.0


def test_specific_heat_const_pres():
    """

    This function tests the spec_heat_const_pressure() function
    """
    cp = gas.spec_heat_const_pressure(temp, pres)
    assert isclose(1040.0, cp, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_specific_heat_const_volume():
    """

    This function tests the spec_heat_const_volume() function
    """
    cv = gas.spec_heat_const_volume(temp, pres)
    assert isclose(743.9, cv, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_ratio_of_specific_heats():
    """

    This function tests the ratio_of_specific_heats() function
    """
    gamma = gas.ratio_specific_heats(temp, pres)
    assert isclose(1.3989, gamma, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_molar_mass():
    """

    This function test the molar_mass() function
    """
    mw = gas.molar_mass(temp, pres)
    assert isclose(28.0134, mw, rel_tol=1.0e-3)
# ==============================================================================
# ==============================================================================
# eof
