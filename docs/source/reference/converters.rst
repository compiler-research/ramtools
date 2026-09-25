Converters
==========

The writers. Each calls ``RAMNTupleRecord::InitializeRefs()``, fills records
through the setters, checks the coordinate order and the longest span as it
goes, and writes ``METADATA`` and the ``headers`` key when the input is
exhausted. The conversion options are described in :doc:`../user/converting`.

samtoramntuple
--------------

.. doxygenfunction:: samtoramntuple

samtoramntuple_split_by_chromosome
----------------------------------

.. doxygenfunction:: samtoramntuple_split_by_chromosome

bamtoramntuple
--------------

.. doxygenfunction:: bamtoramntuple
