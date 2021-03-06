ThermSystems
============

The ThermSystems code repository was created to provide a simple
framework to model the performance of gas-flow thermodynamic systems
such as turbojet-engines, turbofan-engines, turbo-prop systems, and gas
flow loops such as those seen in coal-fired and nuclear power-plants.

The code includes pre-built configurations to model the various thermodynamic
systems.  In addition, a developr or user can also access the lower
level models to predict adiabatic-isentropic process for sub-sonic and
incompressible regimes, as well as models for a diffuser, nozzle,
turbine, heat exhcanger/combustion chamber and piping.

The code contains the following classes

* **Stagnation:** which contains classmethods to be inherited in other classes
* **DiffuserNozzle**: which contains functions that determine the fluid exit conditions leaving a diffuser or nozzle
* **Compressor**: which contains functions that determine the fluid exit conditions leaving a compressor
* **HeatAddition**: which contains functions that determine the fluid exit conditions leaving a combustion chamber
* **Turbine**: which contains functions that determine the fluid exit conditions leaving a turbine
* **Propeller**: which contains functions that determine the fluid exit conditions leaving a propeller