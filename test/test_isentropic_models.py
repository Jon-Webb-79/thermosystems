# Import modules here
from src.isentropic_models import DiffuserNozzle, Stagnation, Compressor
from src.isentropic_models import HeatAddition, Turbine, DiffuserNozzleComponent
from src.isentropic_models import CompressorComponent, HeatAdditionComponent
from src.isentropic_models import TurbineComponent, Propeller, PropellerComponent
from src.isentropic_models import RamJet, TurboJet

from thermo.chemical import Chemical
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
    temp = 800.0
    sos = Stagnation.speed_of_sound(gamma, temp, 28.014)
    assert isclose(576.53, sos, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_mach_number():
    """

    This function tests the mach_number() function
    """
    temp = 800.0
    velocity = 600.0
    molar_mass = 28.014
    mach = Stagnation.mach_number(gamma, temp, molar_mass, velocity)
    assert isclose(1.04, mach, rel_tol=1.0e-2)
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
inlet_velocity = 86.5
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
    assert isclose(21.624, exit_velocity, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_diff_exit_stat_temp():
    """

    This function tests the exit_static_temperature() function
    """
    cp = 1030
    stat_temp = diff.exit_static_temperature(cp, gamma, inlet_static_temperature,
                                             inlet_static_pressure, mach_number,
                                             inlet_velocity)
    assert isclose(stat_temp, 803.157, rel_tol=1.0e-4)
# ------------------------------------------------------------------------------


def test_diff_exit_mach_number():
    """

    This function tests the exit_mach_number() function
    """
    cp = 1030.0
    molar_mass = 28.014
    mach_num = diff.exit_mach_number(cp, gamma, inlet_static_temperature,
                                     inlet_static_pressure, mach_number,
                                     inlet_velocity, molar_mass)
    assert isclose(mach_num, 0.03743, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_diff_exit_static_pressure():
    """

    This function tests the exit_static_pressure() function
    """
    cp = 1030.0
    molar_mass = 28.014
    stat_pres = diff.exit_static_pressure(gamma, inlet_static_pressure,
                                          mach_number, cp, inlet_static_temperature,
                                          inlet_velocity, molar_mass)
    assert isclose(stat_pres, 6083331.69, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_diff_exit_density():
    """

    This function tests the exit_density() function
    """
    mass_flow_rate = 3.0
    density = diff.exit_density(inlet_velocity, mass_flow_rate)
    assert isclose(0.01901, density, rel_tol=1.0e-3)
# ==============================================================================
# ==============================================================================
# Test the Compressor class


compressor_ratio = 3.2
inlet_stag_pressure = 6000000.0
inlet_stag_temperature = 800.0
efficiency = 0.9
inlet_mach_number = 0.1
comp = Compressor(compressor_ratio, efficiency)


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
    assert isclose(vel, 440.72, rel_tol=1.0e-3)
# ==============================================================================
# ==============================================================================
# Test the HeatAddition class


heat = 1000000.0
efficiency = 0.9
he = HeatAddition(efficiency)


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
    stag_temp = he.exit_stagnation_temperature(803.384, heat, 1030.0, 3.0)
    assert isclose(1127.0, stag_temp, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_he_exit_mach_number():
    """

    This function tests the exit_mach_number() function
    """
    mach = he.exit_mach_number(803.384, heat, 1030.0, 3.0,
                               0.03743, gamma)
    assert isclose(mach, 0.044365, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_he_exit_stagnation_pressure():
    """

    This function tests the exit_stagnation_pressure() function
    """
    stag_pres = he.exit_stagnation_pressure(6089300.0,
                                            803.384,
                                            heat, 1030.0, 3.0,
                                            0.03743, gamma)
    assert isclose(stag_pres, 6086893.10, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_he_exit_static_temperature():
    """

    This function tests the exit_static_temperature() function
    """
    stat_temp = he.exit_static_temperature(803.384, heat, 1030.0,
                                           3.0, gamma, 0.03743)
    assert isclose(stat_temp, 1126.565, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_he_exit_static_pressure():
    """

    This function tests the exit_static_pressure() function
    """
    stat_pres = he.exit_static_pressure(6089300.0, 803.384,
                                        heat, 1030.0, 3.0, gamma, 0.03743)
    assert isclose(6080917.69, stat_pres, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_he_exit_velocity():
    """

    This function tests the exit_velocity() function
    """
    velocity = he.exit_velocity(803.384, heat, 1030.0,
                                3.0, 0.03743, gamma, 28.014)
    assert isclose(velocity, 30.352, rel_tol=1.0e-3)
# ==============================================================================
# ==============================================================================
# Test Turbine class


work = 1000.0
turb = Turbine(0.9)


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
    assert isclose(5926127.77, pres, rel_tol=1.0e-3)
# ------------------------------------------------------------------------------


def test_turbine_exit_velocity():
    """

    This function tests the exit_velocity() function
    """
    vel = turb.exit_velocity(inlet_stag_temperature, inlet_mach_number, gamma,
                             work, 1.0, 1000.0, 1.0)
    assert isclose(vel, 304.417, rel_tol=1.0e-3)
# ==============================================================================
# ==============================================================================
# Test the Propeller class


efficiency = 0.9
prop = Propeller(efficiency)
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
    assert isclose(velocity, 206.719, rel_tol=1.0e-3)
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
    inlet_velocity = 86.5
    molar_mass = 28.014
    diff_comp = DiffuserNozzleComponent(efficiency, inlet_area, exit_area)
    exit_cond = diff_comp.outlet_conditions(gamma, 1030.0, molar_mass, inlet_static_temperature,
                                            inlet_static_pressure, mach_number,
                                            mass_flow_rate, inlet_velocity)
    assert isclose(exit_cond['static_pressure'], 6083331.69, rel_tol=1.0e-3)
    assert isclose(exit_cond['static_temperature'], 803.157, rel_tol=1.0e-3)
    assert isclose(exit_cond['stagnation_pressure'], 6089300.0, rel_tol=1.0e-3)
    assert isclose(exit_cond['stagnation_temperature'], 803.384, rel_tol=1.0e-3)
    assert isclose(exit_cond['velocity'], 21.624, rel_tol=1.0e-3)
    assert isclose(exit_cond['mach_number'], 0.03743, rel_tol=1.0e-3)
    assert isclose(exit_cond['density'], 0.01901, rel_tol=1.0e-3)
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

    comp = CompressorComponent(compressor_ratio, efficiency)
    exit_cond = comp.outlet_conditions(gamma, specific_heat, molar_mass,
                                       inlet_mach_number, mass_flow_rate,
                                       inlet_stag_pressure,
                                       inlet_stag_temperature)
    assert isclose(exit_cond['static_pressure'], 19005747.0, rel_tol=1.0e-3)
    assert isclose(exit_cond['static_temperature'], 1147.07, rel_tol=1.0e-3)
    assert isclose(exit_cond['stagnation_pressure'], 19200000.0, rel_tol=1.0e-3)
    assert isclose(exit_cond['stagnation_temperature'], 1150.4, rel_tol=1.0e-3)
    assert isclose(exit_cond['velocity'], 440.72, rel_tol=1.0e-3)
    assert isclose(exit_cond['mach_number'], 0.120615, rel_tol=1.0e-3)
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
    inlet_mach_number = 0.03743
    gamma = 1.4
    mass_flow_rate = 3.0
    specific_heat = 1030.0
    molar_mass = 28.014

    he = HeatAdditionComponent(efficiency)
    exit_cond = he.outlet_conditions(gamma, specific_heat, molar_mass,
                                     inlet_mach_number, mass_flow_rate,
                                     6089300.0, 803.384, heat)
    assert isclose(exit_cond['power'], 1111111.11, rel_tol=1.0e-3)
    assert isclose(exit_cond['static_temperature'], 1126.565, rel_tol=1.0e-3)
    assert isclose(exit_cond['static_pressure'], 6080917.69, rel_tol=1.0e-3)
    assert isclose(exit_cond['stagnation_pressure'], 6086893.14, rel_tol=1.0e-3)
    assert isclose(exit_cond['stagnation_temperature'], 1127.008, rel_tol=1.0e-3)
    assert isclose(exit_cond['velocity'], 30.352, rel_tol=1.0e-3)
    assert isclose(exit_cond['mach_number'], 0.044365, rel_tol=1.0e-3)
# ==============================================================================
# ==============================================================================


def test_turbine_component():
    """

    This function tests the TurbineComponent() class
    """
    turbine_work = 1000.0
    efficiency = 0.9
    inlet_mach_number = 0.1
    gamma = 1.4
    mass_flow_rate = 1000.0
    specific_heat = 1.0
    molar_mass = 1.0

    turb = TurbineComponent(efficiency)
    exit_cond = turb.outlet_conditions(gamma, specific_heat, molar_mass,
                                       inlet_mach_number, mass_flow_rate,
                                       inlet_stag_pressure,
                                       inlet_stag_temperature,
                                       turbine_work)
    assert isclose(exit_cond['extracted_work'], 1111.11, rel_tol=1.0e-3)
    assert isclose(exit_cond['static_temperature'], 797.296, rel_tol=1.0e-3)
    assert isclose(exit_cond['static_pressure'], 5926127.77, rel_tol=1.0e-3)
    assert isclose(exit_cond['stagnation_pressure'], 5967655.0, rel_tol=1.0e-3)
    assert isclose(exit_cond['stagnation_temperature'], 798.88, rel_tol=1.0e-3)
    assert isclose(exit_cond['velocity'], 304.417, rel_tol=1.0e-3)
    assert isclose(exit_cond['mach_number'], 0.099928707, rel_tol=1.0e-3)
# ==============================================================================
# ==============================================================================


def test_propeller_component():
    """

    This function tests the PropellerComponent() class
    """
    efficiency = 0.9
    inlet_mach_number = 0.1
    gamma = 1.4
    mass_flow_rate = 100.0
    specific_heat = 1.0
    molar_mass = 1.0
    work = 1000.0
    inlet_stag_temperature = 345.0
    inlet_stag_pressure = 101325.0

    propc = PropellerComponent(efficiency)
    exit_cond = propc.outlet_conditions(work, inlet_stag_temperature,
                                        mass_flow_rate, specific_heat,
                                        inlet_stag_pressure, gamma,
                                        inlet_mach_number, molar_mass)
    assert isclose(exit_cond['total_work'], 1111.11, rel_tol=1.0e-3)
    assert isclose(exit_cond['static_temperature'], 355.37, rel_tol=1.0e-3)
    assert isclose(exit_cond['static_pressure'], 91050.91, rel_tol=1.0e-3)
    assert isclose(exit_cond['stagnation_pressure'], 91711.0, rel_tol=1.0e-3)
    assert isclose(exit_cond['stagnation_temperature'], 356.11, rel_tol=1.0e-3)
    assert isclose(exit_cond['velocity'], 206.719, rel_tol=1.0e-3)
    assert isclose(exit_cond['mach_number'], 0.1016405, rel_tol=1.0e-3)
# ==============================================================================
# ==============================================================================
# Test RamJet


species = 'nitrogen'
gas = Chemical(species)

dif_eff = 0.95
diff_inlet_area = 0.0729
diff_outlet_area = 0.462
heat_eff = 0.98
nozzle_eff = 0.95
nozzle_inlet_area = 0.462
nozzle_outlet_area = 0.074
temp = 81.05  # Kelvins
pres = 49300.0  # Pascals
gas.calculate(T=temp, P=pres)
density = gas.rho
gamma1 = gas.Cp / gas.Cvg
mw = gas.MW
sos = (gamma1 * ((1000.0 * 8.314) / mw) * temp) ** 0.5
velocity = 400.0  # m/s
mach = velocity / sos
mdot = density * velocity * diff_inlet_area  # kg/s
jet = RamJet(dif_eff, diff_inlet_area, diff_outlet_area, heat_eff,
             nozzle_eff, nozzle_inlet_area, nozzle_outlet_area, species)


def test_ram_jet():
    """

    This function tests the RamJet class
    """
    dif, heat, noz = jet.performance(temp, pres, mach, velocity, mdot, 100000.0)
    assert isclose(dif['static_temperature'], 152.297, rel_tol=1.0e-3)
    assert isclose(dif['static_pressure'], 448365.0, rel_tol=1.0e-3)
    assert isclose(dif['stagnation_pressure'], 468436.0, rel_tol=1.0e-3)
    assert isclose(dif['stagnation_temperature'], 154.21, rel_tol=1.0e-3)
    assert isclose(dif['velocity'], 63.116, rel_tol=1.0e-3)
    assert isclose(dif['mach_number'], 0.250906235, rel_tol=1.0e-3)
    assert isclose(dif['density'], 2.0858, rel_tol=1.0e-3)

    assert isclose(heat['static_temperature'], 153.836, rel_tol=1.0e-3)
    assert isclose(heat['static_pressure'], 448365.38, rel_tol=1.0e-3)
    assert isclose(heat['stagnation_pressure'], 468436.95, rel_tol=1.0e-3)
    assert isclose(heat['stagnation_temperature'], 155.79, rel_tol=1.0e-3)
    assert isclose(heat['velocity'], 63.817, rel_tol=1.0e-3)
    assert isclose(heat['mach_number'], 0.2524, rel_tol=1.0e-3)
    assert isclose(heat['power'], 102040.81, rel_tol=1.0e-3)

    assert isclose(noz['static_temperature'], 79.29, rel_tol=1.0e-3)
    assert isclose(noz['static_pressure'], 44051.23, rel_tol=1.0e-3)
    assert isclose(noz['stagnation_pressure'], 467406.27, rel_tol=1.0e-3)
    assert isclose(noz['stagnation_temperature'], 155.699, rel_tol=1.0e-3)
    assert isclose(noz['velocity'], 398.42, rel_tol=1.0e-3)
    assert isclose(noz['mach_number'], 2.1958, rel_tol=1.0e-3)
    assert isclose(noz['density'], 2.0629, rel_tol=1.0e-3)
# ==============================================================================
# ==============================================================================
# Test TurboJet


def test_turbo_jet():
    species = 'nitrogen'
    gas = Chemical(species)

    dif_eff = 0.95
    diff_inlet_area = 0.0729
    diff_outlet_area = 0.462

    comp_ratio = 1.2
    comp_eff = 0.9

    heat_eff = 0.98
    nozzle_eff = 0.95

    turb_eff = 0.95
    turb_work = 89286.33

    nozzle_inlet_area = 0.462
    nozzle_outlet_area = 0.12

    temp = 81.05  # Kelvins
    pres = 49300.0  # Pascals
    gas.calculate(T=temp, P=pres)
    density = gas.rho
    gamma1 = gas.Cp / gas.Cvg
    mw = gas.MW
    sos = (gamma1 * ((1000.0 * 8.314) / mw) * temp) ** 0.5
    velocity = 50.0  # m/s
    mach = velocity / sos
    mdot = density * velocity * diff_inlet_area  # kg/s
    jet = TurboJet(dif_eff, diff_inlet_area, diff_outlet_area, comp_ratio, comp_eff,
                   heat_eff, turb_eff, nozzle_eff, nozzle_inlet_area, nozzle_outlet_area,
                   species)
    dif, comp, heat, turb, noz = jet.performance(temp, pres, mach, mdot, velocity, 10000.0, turb_work)
    assert isclose(dif['static_temperature'], 82.163, rel_tol=1.0e-3)
    assert isclose(dif['stagnation_pressure'], 51777.00, rel_tol=1.0e-3)
    assert isclose(dif['static_pressure'], 51710.97, rel_tol=1.0e-3)
    assert isclose(dif['stagnation_temperature'], 82.193, rel_tol=1.0e-3)
    assert isclose(dif['velocity'], 7.889, rel_tol=1.0e-3)
    assert isclose(dif['mach_number'], 0.0427, rel_tol=1.0e-3)
    assert isclose(dif['density'], 2.0858, rel_tol=1.0e-3)

    assert isclose(comp['static_temperature'], 87.042, rel_tol=1.0e-3)
    assert isclose(comp['stagnation_pressure'], 62132.40, rel_tol=1.0e-3)
    assert isclose(comp['static_pressure'], 62048.22, rel_tol=1.0e-3)
    assert isclose(comp['stagnation_temperature'], 87.07, rel_tol=1.0e-3)
    assert isclose(comp['velocity'], 8.359, rel_tol=1.0e-3)
    assert isclose(comp['mach_number'], 0.043956, rel_tol=1.0e-3)
    assert isclose(comp['work'], 42854.99, rel_tol=1.0e-3)

    assert isclose(heat['static_temperature'], 88.308, rel_tol=1.0e-3)
    assert isclose(heat['stagnation_pressure'], 62131.186, rel_tol=1.0e-3)
    assert isclose(heat['static_pressure'], 62047.22, rel_tol=1.0e-3)
    assert isclose(heat['stagnation_temperature'], 88.342, rel_tol=1.0e-3)
    assert isclose(heat['velocity'], 8.481, rel_tol=1.0e-3)
    assert isclose(heat['mach_number'], 0.044267, rel_tol=1.0e-3)
    assert isclose(heat['power'], 10204.081, rel_tol=1.0e-3)

    assert isclose(turb['static_temperature'], 76.416, rel_tol=1.0e-3)
    assert isclose(turb['stagnation_pressure'], 36381.49, rel_tol=1.0e-3)
    assert isclose(turb['static_pressure'], 36338.84, rel_tol=1.0e-3)
    assert isclose(turb['stagnation_temperature'], 76.442, rel_tol=1.0e-3)
    assert isclose(turb['velocity'], 7.336, rel_tol=1.0e-3)
    assert isclose(turb['mach_number'], 0.04117, rel_tol=1.0e-3)
    assert isclose(turb['extracted_work'], 93985.61, rel_tol=1.0e-3)

    assert isclose(noz['static_temperature'], 76.057, rel_tol=1.0e-3)
    assert isclose(noz['stagnation_pressure'], 36379.82, rel_tol=1.0e-3)
    assert isclose(noz['static_pressure'], 35744.16, rel_tol=1.0e-3)
    assert isclose(noz['stagnation_temperature'], 76.44, rel_tol=1.0e-3)
    assert isclose(noz['velocity'], 28.24, rel_tol=1.0e-3)
    assert isclose(noz['mach_number'], 0.1588, rel_tol=1.0e-3)
    assert isclose(noz['density'], 2.243, rel_tol=1.0e-3)
# ==============================================================================
# ==============================================================================
# eof
