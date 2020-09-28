# Import modules here
from src.isentropic_models import DiffuserNozzle, Stagnation, Compressor

from math import isclose
# ==============================================================================
# ==============================================================================
# Date:    September 17, 2020
# Purpose: This code contains functions that test the functions and classes
#          in the isentropic_models.py file

# Source Code Metadata
__author__ = "Jonathan A. Webb"
__copyright__ = "Copyright 2020, Jon Webb Inc."
__version__ = "1.0"
# ==============================================================================
# ==============================================================================
# Test the Stagnation class
gamma = 1.4
mach_number = 0.15


def test_stagnation_pressure():
    """

    This function tests the stagnation_pressure() function
    """
    static_pressure = 6000000.0
    stag_pres = Stagnation.stagnation_pressure(static_pressure, mach_number, gamma)
    assert isclose(stag_pres, 6095032.0, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_stagnation_temperature():
    """

    This function tests the stagnation_temperature() function
    """
    static_temperature = 800.0
    stag_temp = Stagnation.stagnation_temperature(static_temperature,
                                                  mach_number, gamma)
    assert isclose(stag_temp, 803.59, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_static_pressure():
    """

    This function tests the static_pressure() function
    """
    stag_pres = 6095032.0
    stat_pres = Stagnation.static_pressure(stag_pres, mach_number, gamma)
    assert isclose(stat_pres, 6000000.0, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_static_temperature():
    """

    This function tests the static_temperature() function
    """
    stag_temp = 803.59
    stat_temp = Stagnation.static_temperature(stag_temp, mach_number, gamma)
    assert isclose(stat_temp, 800.0, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_pressure_based_mach():
    """

    This function tests the pressure_based_mach() function
    """
    stag_pressure = 6095032.0
    stat_pressure = 6000000.0
    mach = Stagnation.pressure_based_mach(stag_pressure, stat_pressure, gamma)
    assert isclose(mach, 0.15, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_temperature_based_mach():
    """


    This function tests the temperature_based_mach() function
    """
    stag_temp = 803.59
    stat_temp = 800.0
    mach = Stagnation.temperature_based_mach(stag_temp, stat_temp, gamma)
    assert isclose(mach, 0.15, rel_tol=1.0e-2)
# ------------------------------------------------------------------------------


def test_speed_of_sound():
    """

    This function tests the speed_of_sound() function
    """
    temp = 273.0
    sos = Stagnation.speed_of_sound(gamma, temp, 0.0289645)
    assert isclose(331.22, sos, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_mach_number():
    """

    This function tests the mach_number() function
    """
    temp = 273.0
    velocity = 662.44
    molar_mass = 0.0289645
    mach = Stagnation.mach_number(gamma, temp, molar_mass, velocity)
    assert isclose(2.0, mach, rel_tol=1.0e-2)
# ------------------------------------------------------------------------------


def test_temp_exit_mach():
    """

    This function tests the temp_exit_mach() function
    """
    entrance_mach = 0.19
    exit_mach = 0.2151
    inlet_temp = 800.0
    exit_temp = 1000.0
    remainder = Stagnation.temp_exit_mach(exit_mach, entrance_mach,
                                          gamma, inlet_temp, exit_temp)
    assert isclose(-9.00597e-5, remainder, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_press_exit_mach():
    """

    This function tests the press_exit_mach() function
    """
    entrance_mach = 0.25
    exit_mach = 0.05
    inlet_press = 6000000.0
    exit_press = 6236000.0
    remainder = Stagnation.press_exit_mach(exit_mach, entrance_mach, gamma,
                                           inlet_press, exit_press)
    assert isclose(remainder, -7.58682e-5, rel_tol=1.0e-3)
# ==============================================================================
# ==============================================================================
# Test the Diffuser class


efficiency = 0.94
inlet_static_pressure = 6000000.0
inlet_static_temperature = 800.0
inlet_area = 1.8241
exit_area = 7.2965
inlet_velocity = 3.4
diff = DiffuserNozzle(efficiency, inlet_area, exit_area)


def test_diff_exit_stag_pressure():
    """

    Test the exit_stagnation_pressure() function for the Diffuser class
    """
    exit_stag_pres = diff.exit_stagnation_pressure(gamma, inlet_static_pressure,
                                                   mach_number)
    assert isclose(6089300.0, exit_stag_pres, rel_tol=1.0e-2)
# ------------------------------------------------------------------------------


def test_diff_exit_stag_temperature():
    """

    Test the exit_stagnation_temperature() function for the Diffuser class
    """
    exit_stag_temp = diff.exit_stagnation_temperature(gamma, inlet_static_temperature,
                                                      inlet_static_pressure, mach_number)
    assert isclose(exit_stag_temp, 803.384, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_diff_exit_velocity():
    """

    Test the exit_velocity() function for the Diffuser class
    """
    exit_velocity = diff.exit_velocity(inlet_velocity)
    assert isclose(0.8499, exit_velocity, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_diff_exit_stat_temp():
    """

    This function tests the exit_static_temperature() function
    """
    cp = 1.0
    stat_temp = diff.exit_static_temperature(cp, gamma, inlet_static_temperature,
                                             inlet_static_pressure, mach_number,
                                             inlet_velocity)
    assert isclose(stat_temp, 803.022, rel_tol=1.0e-4)
# ------------------------------------------------------------------------------


def test_diff_exit_mach_number():
    """

    This function tests the exit_mach_number() function
    """
    cp = 1.0
    molar_mass = 0.0289645
    mach_num = diff.exit_mach_number(cp, gamma, inlet_static_temperature,
                                     inlet_static_pressure, mach_number,
                                     inlet_velocity, molar_mass)
    assert isclose(mach_num, 0.001496, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_diff_exit_static_pressure():
    """

    This function tests the exit_static_pressure() function
    """
    cp = 1.0
    molar_mass = 0.0289645
    stat_pres = diff.exit_static_pressure(gamma, inlet_static_pressure,
                                          mach_number, cp, inlet_static_temperature,
                                          inlet_velocity, molar_mass)
    assert isclose(stat_pres, 6089291.0, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_diff_exit_density():
    """

    This function tests the exit_density() function
    """
    mass_flow_rate = 3.0
    density = diff.exit_density(inlet_velocity, mass_flow_rate)
    assert isclose(0.4837, density, rel_tol=1.0e-3)
# ==============================================================================
# ==============================================================================
# Test the Compressor class


compressor_ratio = 3.2
inlet_stag_pressure = 6000000.0
inlet_stag_temperature = 800.0
efficiency = 0.9
inlet_mach_number = 0.1
exit_area = 1.0
comp = Compressor(compressor_ratio, efficiency, exit_area)


def test_comp_exit_stag_pressure():
    """

    This function tests the exit_stagnation_pressure() function
    """
    exit_stag_pres = comp.exit_stagnation_pressure(inlet_stag_pressure)
    assert isclose(19200000.0, exit_stag_pres, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_comp_exit_stag_temperature():
    """

    This function tests the exit_stagnation_temperature() function
    """
    exit_stag_temp = comp.exit_stagnation_temperature(inlet_stag_temperature, gamma)
    assert isclose(exit_stag_temp, 1150.4, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_comp_work():
    """

    This function tests the compressor_work() function
    """
    mass_flow_rate = 1.10
    specific_heat = 1.0
    comp_work = comp.compressor_work(mass_flow_rate, specific_heat,
                                     inlet_stag_temperature, gamma)
    assert isclose(428.278, comp_work, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_comp_mach_number():
    """

    This function tests the compressor_mach_number() function
    """
    mach = comp.exit_mach_number(inlet_mach_number, inlet_stag_temperature,
                                 gamma, 0.001, 0.999)
    assert isclose(mach[0], 0.120615, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_comp_static_pressure():
    """

    This function tests the compressor exit_static_pressure() function
    """
    press = comp.exit_static_pressure(inlet_stag_pressure, inlet_stag_temperature,
                                      gamma, inlet_mach_number)
    assert isclose(press, 19005747.0, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_comp_static_temperature():
    """

    This function tests the compressor exit_static_temperature() function
    """
    temp = comp.exit_static_temperature(inlet_stag_temperature, gamma, inlet_mach_number)
    assert isclose(temp, 1147.07, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_comp_fluid_velocity():
    """

    This function tests the compressor exit_velocity() function
    """
    vel = comp.exit_velocity(inlet_stag_temperature, gamma, inlet_mach_number,
                             1.0)
    assert isclose(vel, 13.936, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_comp_fluid_density():
    """

    This function tests the compressor exit_density() function
    """
    den = comp.exit_density(1.0, inlet_stag_temperature, inlet_mach_number, 1.0,
                            gamma)
    assert isclose(den, 0.07175, rel_tol=1.0e-3)
# ==============================================================================
# ==============================================================================
# eof
