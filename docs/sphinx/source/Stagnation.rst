Stagnation
==========
This section documents the ``Stagnation`` class encoded in the
``isentropic_models.py`` file.  This class contains models for
all generic adiabatic-isentropic process that may be relevant
to multiple components.  The functions are encoded as classmethods
owing to the fact that they share no data or functionality with
other member functions, but need to be in a class to allow for
proper inheretance in other classes.

.. autoclass:: src.isentropic_models.Stagnation
   :members: