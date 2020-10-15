Diffuser/NozzleComponent
========================
This section documents the ``DiffuserNozzleComponent`` class encoded in the
``isentropic_models.py`` file.  This class contains models for
all generic adiabatic-isentropic process that may be relevant
to the determination of exit properties for a fluid leaving
a diffuser or a nozzle.  This class assumes sub-sonic, incompressible
gas flow.  This class reduces the complexity of using the DiffuserNozzle
class by allowing a user to obtain all outlet properties with a single
function call.

.. autoclass:: src.isentropic_models.DiffuserNozzleComponent
   :members: