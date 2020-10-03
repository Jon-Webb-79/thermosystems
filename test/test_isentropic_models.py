# Import modules here
from src.isentropic_models import DiffuserNozzle, Stagnation, Compressor
from src.isentropic_models import HeatAddition, Turbine, DiffuserNozzleComponent
from src.isentropic_models import CompressorComponent, HeatAdditionComponent
from src.isentropic_models import TurbineComponent, Propeller, PropellerComponent

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
    assert isclose(mach, 0.120615, rel_tol=1.0e-3)
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
    den = comp.exit_density(1.10, inlet_stag_temperature, inlet_mach_number, 1.0,
                            gamma)
    assert isclose(den, 0.07892, rel_tol=1.0e-3)
# ==============================================================================
# ==============================================================================
# Test the HeatAddition class


heat = 1000000.0
efficiency = 0.9
he = HeatAddition(efficiency, 10.0)


def test_he_input_power():
    """

    This function tests the input_power() function
    """
    power = he.input_power(heat)
    assert isclose(power, 1111111.11, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_he_exit_stagnation_temperature():
    """

    This function tests the exit_stagnation_temperature() function
    """
    stag_temp = he.exit_stagnation_temperature(inlet_stag_temperature, heat, 1.0, 1000.0)
    assert isclose(1800.0, stag_temp, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_he_exit_mach_number():
    """

    This function tests the exit_mach_number() function
    """
    mach = he.exit_mach_number(inlet_stag_temperature, heat, 1.0, 1000.0,
                               inlet_mach_number, gamma)
    assert isclose(mach, 0.1525, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_he_exit_stagnation_pressure():
    """

    This function tests the exit_stagnation_pressure() function
    """
    stag_pres = he.exit_stagnation_pressure(inlet_stag_pressure,
                                            inlet_stag_temperature,
                                            heat, 1.0, 1000.0,
                                            inlet_mach_number, gamma)
    assert isclose(stag_pres, 5946850.0, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_he_exit_static_temperature():
    """

    This function tests the exit_static_temperature() function
    """
    stat_temp = he.exit_static_temperature(inlet_stag_temperature, heat, 1.0,
                                           1000.0, gamma, inlet_mach_number)
    assert isclose(stat_temp, 1791.6, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_he_exit_static_pressure():
    """

    This function tests the exit_static_pressure() function
    """
    stat_pres = he.exit_static_pressure(inlet_stag_pressure, inlet_stag_temperature,
                                        heat, 1.0, 1000.0, gamma, inlet_mach_number)
    assert isclose(5903279.0, stat_pres, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_he_exit_velocity():
    """

    This function tests the exit_velocity() function
    """
    velocity = he.exit_velocity(inlet_stag_temperature, heat, 1.0,
                                1000.0, inlet_mach_number, gamma, 1.0)
    assert isclose(velocity, 22.029, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_he_exit_density():
    """

    This function tests the exit_density() function
    """
    density = he.exit_density(inlet_stag_temperature, heat, 1.0,
                              1000.0, inlet_mach_number, gamma, 1.0)
    assert isclose(density, 4.539, rel_tol=1.0e-3)
# ==============================================================================
# ==============================================================================
# Test Turbine class


work = 1000.0
turb = Turbine(0.9, 10.0)


def test_turbine_work():
    """

    :This function tests the work_extraction() function
    """
    total_work = turb.work_extraction(work)
    assert isclose(1111.11, total_work, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_turbine_stagnation_temperature():
    """

    This function tests the exit_stagnation_temperature() function
    """
    temp = turb.exit_stagnation_temperature(work, 1.0, 1000.0, inlet_stag_temperature)
    assert isclose(798.88, temp, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_turbine_stagnation_pressure():
    """

    This function tests the exit_stagnation_pressure() function
    """
    pres = turb.exit_stagnation_pressure(inlet_stag_pressure, inlet_stag_temperature,
                                         1.0, 1000.0, work, gamma)
    assert isclose(pres, 5967655.0, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_turbine_exit_mach_number():
    """

    This function tests the exit_mach_number() function
    """
    mach = turb.exit_mach_number(inlet_stag_temperature, inlet_mach_number, gamma,
                                 work, 1.0, 1000.0)
    assert isclose(mach, 0.09928707, rel_tol=1.0e-2)
# ------------------------------------------------------------------------------


def test_turbine_exit_static_temperature():
    """

    This function tests the exit_static_temperature() function
    """
    temp = turb.exit_static_temperature(inlet_stag_temperature, inlet_mach_number,
                                        gamma, work, 1.0, 1000.0)
    assert isclose(797.296, temp, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_turbine_exit_static_pressure():
    """

    This function tests the exit_static_pressure() function
    """
    pres = turb.exit_static_pressure(inlet_stag_temperature, inlet_mach_number,
                                     gamma, work, 1.0, 1000.0, inlet_stag_pressure)
    assert isclose(5958247.0, pres, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_turbine_exit_velocity():
    """

    This function tests the exit_velocity() function
    """
    vel = turb.exit_velocity(inlet_stag_temperature, inlet_mach_number, gamma,
                             work, 1.0, 1000.0, 1.0)
    assert isclose(vel, 9.626, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_turbine_exit_density():
    """

    This function tests the exit_density() function
    """
    den = turb.exit_density(inlet_stag_temperature, inlet_mach_number, gamma,
                             work, 1.0, 1000.0, 1.0)
    assert isclose(10.387, den, rel_tol=1.0e-3)
# ==============================================================================
# ==============================================================================
# Test the Propeller class


efficiency = 0.9
prop = Propeller(efficiency, 10.0)
inlet_stag_temp = 345.0
inlet_stag_pres = 101325.0
usable_work = 1000.0
mass_flow_rate = 100.0
specific_heat = 1.0

def test_propeller_work():
    """

    This function tests the propeller_work() function
    """
    work = prop.propeller_work(usable_work)
    assert isclose(work, 1111.11, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_propeller_exit_stag_temp():
    """

    This function tests the exit_stagnation_temperature() function
    """
    stag_temp = prop.exit_stagnation_temperature(inlet_stag_temp, mass_flow_rate,
                                                 specific_heat, work)
    assert isclose(stag_temp, 356.11, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_propeller_exit_stag_pres():
    """

    This function tests the exit_stagnation_pressure() function
    """
    stag_pres = prop.exit_stagnation_pressure(inlet_stag_pres, inlet_stag_temp,
                                              gamma, mass_flow_rate, specific_heat,
                                              work)
    assert isclose(91711.0, stag_pres, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_propeller_exit_mach_number():
    """

    This function tests the exit_mach_number() function
    """
    mach = prop.exit_mach_number(inlet_stag_temp, inlet_mach_number, mass_flow_rate,
                                 specific_heat, work, gamma)
    assert isclose(mach, 0.1016405, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_propeller_exit_static_temp():
    """

    This function tests the exit_static_temperature() function
    """
    stat_temp = prop.exit_static_temperature(inlet_stag_temp, inlet_mach_number,
                                             mass_flow_rate, specific_heat,
                                             work, gamma)
    assert isclose(stat_temp, 355.37, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_propeller_exit_static_pres():
    """

    This function tests the exit_static_pressure() function
    """
    stat_pres = prop.exit_static_pressure(inlet_stag_temp, inlet_mach_number,
                                          mass_flow_rate, specific_heat,
                                          work, gamma,
                                          inlet_stag_pres)
    assert isclose(stat_pres, 91050.91, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_propeller_exit_velocity():
    """

    This function tests the exit_velocity() function
    """
    velocity = prop.exit_velocity(inlet_stag_temp, inlet_mach_number,
                                  mass_flow_rate, specific_heat,
                                  work, gamma, 1.0)
    assert isclose(velocity, 6.537, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_propeller_exit_density():
    """

    This function tests the exit_density() function
    """
    density = prop.exit_density(inlet_stag_temp, inlet_mach_number,
                                mass_flow_rate, specific_heat,
                                work, gamma, 1.0)
    assert isclose(density, 1.5297, rel_tol=1.0e-3)
# ==============================================================================
# ==============================================================================
# Test the DiffuserNozzleComponent class


def test_diffusernozzle_component():
    """

    This function tests the DiffuserNozzleComponent class
    """
    mass_flow_rate = 3.0
    gamma = 1.4
    mach_number = 0.15
    efficiency = 0.94
    inlet_static_pressure = 6000000.0
    inlet_static_temperature = 800.0
    inlet_area = 1.8241
    exit_area = 7.2965
    inlet_velocity = 3.4
    molar_mass = 0.0289645
    diff_comp = DiffuserNozzleComponent(efficiency, inlet_area, exit_area)
    exit_cond = diff_comp.outlet_conditions(gamma, 1.0, molar_mass, inlet_static_temperature,
                                            inlet_static_pressure, mach_number,
                                            mass_flow_rate, inlet_velocity)
    assert isclose(exit_cond['static_pressure'], 6088971.0, rel_tol=1.0e-3)
    assert isclose(exit_cond['static_temperature'], 803.022, rel_tol=1.0e-3)
    assert isclose(exit_cond['stagnation_pressure'], 6089300.0, rel_tol=1.0e-3)
    assert isclose(exit_cond['stagnation_temperature'], 803.384, rel_tol=1.0e-3)
    assert isclose(exit_cond['velocity'], 0.8499, rel_tol=1.0e-3)
    assert isclose(exit_cond['mach_number'], 0.001496, rel_tol=1.0e-3)
    assert isclose(exit_cond['density'], 0.4837, rel_tol=1.0e-3)
# ==============================================================================
# ==============================================================================
# Test the DiffuserNozzleComponent class


def test_compressor_component():
    """

    This function tests the CompressorComponent class
    """
    compressor_ratio = 3.2
    efficiency = 0.9
    inlet_mach_number = 0.1
    exit_area = 1.0
    gamma = 1.4
    mass_flow_rate = 1.10
    specific_heat = 1.0
    molar_mass = 1.0

    comp = CompressorComponent(compressor_ratio, efficiency, exit_area)
    exit_cond = comp.outlet_conditions(gamma, specific_heat, molar_mass,
                                       inlet_mach_number, mass_flow_rate,
                                       inlet_stag_pressure,
                                       inlet_stag_temperature)

    assert isclose(exit_cond['static_pressure'], 19005747.0, rel_tol=1.0e-3)
    assert isclose(exit_cond['static_temperature'], 1147.07, rel_tol=1.0e-3)
    assert isclose(exit_cond['stagnation_pressure'], 19200000.0, rel_tol=1.0e-3)
    assert isclose(exit_cond['stagnation_temperature'], 1150.4, rel_tol=1.0e-3)
    assert isclose(exit_cond['velocity'], 13.936, rel_tol=1.0e-3)
    assert isclose(exit_cond['mach_number'], 0.120615, rel_tol=1.0e-3)
    assert isclose(exit_cond['density'], 0.07892, rel_tol=1.0e-3)
    assert isclose(exit_cond['work'], 428.278, rel_tol=1.0e-3)
# ==============================================================================
# ==============================================================================
# Test the HeatAdditionComponent class


def test_heataddition_component():
    """

    This function tests the HeatAdditionComponent class
    """
    heat = 1000000.0
    efficiency = 0.9
    inlet_mach_number = 0.1
    exit_area = 10.0
    gamma = 1.4
    mass_flow_rate = 1000.0
    specific_heat = 1.0
    molar_mass = 1.0

    he = HeatAdditionComponent(efficiency, exit_area)
    exit_cond = he.outlet_conditions(gamma, specific_heat, molar_mass,
                                     inlet_mach_number, mass_flow_rate,
                                     inlet_stag_pressure,
                                     inlet_stag_temperature, heat)

    assert isclose(exit_cond['power'], 1111111.11, rel_tol=1.0e-3)
    assert isclose(exit_cond['static_temperature'], 1791.6, rel_tol=1.0e-3)
    assert isclose(exit_cond['static_pressure'], 5903279.0, rel_tol=1.0e-3)
    assert isclose(exit_cond['stagnation_pressure'], 5946850.0, rel_tol=1.0e-3)
    assert isclose(exit_cond['stagnation_temperature'], 1800.0, rel_tol=1.0e-3)
    assert isclose(exit_cond['velocity'], 22.029, rel_tol=1.0e-3)
    assert isclose(exit_cond['mach_number'], 0.1525, rel_tol=1.0e-3)
    assert isclose(exit_cond['density'], 4.539, rel_tol=1.0e-3)
# ==============================================================================
# ==============================================================================


def test_turbine_component():
    """

    This function tests the TurbineComponent() class
    """
    turbine_work = 1000.0
    efficiency = 0.9
    inlet_mach_number = 0.1
    exit_area = 10.0
    gamma = 1.4
    mass_flow_rate = 1000.0
    specific_heat = 1.0
    molar_mass = 1.0

    turb = TurbineComponent(efficiency, exit_area)
    exit_cond = turb.outlet_conditions(gamma, specific_heat, molar_mass,
                                       inlet_mach_number, mass_flow_rate,
                                       inlet_stag_pressure,
                                       inlet_stag_temperature,
                                       turbine_work)
    assert isclose(exit_cond['extracted_work'], 1111.11, rel_tol=1.0e-3)
    assert isclose(exit_cond['static_temperature'], 797.296, rel_tol=1.0e-3)
    assert isclose(exit_cond['static_pressure'], 5958247.0, rel_tol=1.0e-3)
    assert isclose(exit_cond['stagnation_pressure'], 5967655.0, rel_tol=1.0e-3)
    assert isclose(exit_cond['stagnation_temperature'], 798.88, rel_tol=1.0e-3)
    assert isclose(exit_cond['velocity'], 9.626, rel_tol=1.0e-3)
    assert isclose(exit_cond['mach_number'], 0.099928707, rel_tol=1.0e-3)
    assert isclose(exit_cond['density'], 10.387, rel_tol=1.0e-3)
# ==============================================================================
# ==============================================================================


def test_propeller_component():
    """

    This function tests the PropellerComponent() class
    """
    efficiency = 0.9
    inlet_mach_number = 0.1
    exit_area = 10.0
    gamma = 1.4
    mass_flow_rate = 100.0
    specific_heat = 1.0
    molar_mass = 1.0
    work = 1000.0
    inlet_stag_temperature = 345.0
    inlet_stag_pressure = 101325.0


    propc = PropellerComponent(efficiency, exit_area)
    exit_cond = propc.outlet_conditions(work, inlet_stag_temperature,
                                        mass_flow_rate, specific_heat,
                                        inlet_stag_pressure, gamma,
                                        inlet_mach_number, molar_mass)
    assert isclose(exit_cond['total_work'], 1111.11, rel_tol=1.0e-3)
    assert isclose(exit_cond['static_temperature'], 355.37, rel_tol=1.0e-3)
    assert isclose(exit_cond['static_pressure'], 91050.91, rel_tol=1.0e-3)
    assert isclose(exit_cond['stagnation_pressure'], 91711.0, rel_tol=1.0e-3)
    assert isclose(exit_cond['stagnation_temperature'], 356.11, rel_tol=1.0e-3)
    assert isclose(exit_cond['velocity'], 6.537, rel_tol=1.0e-3)
    assert isclose(exit_cond['mach_number'], 0.1016405, rel_tol=1.0e-3)
    assert isclose(exit_cond['density'], 1.5297, rel_tol=1.0e-3)
# ==============================================================================
# ==============================================================================
# eof
