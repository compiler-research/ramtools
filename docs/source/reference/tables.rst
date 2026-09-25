Name tables and order check
===========================

Both are owned by ``RAMNTupleRecord`` as static members and are serialised
into the ``METADATA`` ntuple. After ``OpenRAMFile()`` you reach the tables
through ``GetRnameRefs()`` and ``GetRnextRefs()``, and the order through
``IsCoordinateSorted()``.

RAMNTupleRefs
-------------

.. doxygenclass:: RAMNTupleRefs

.. doxygenfunction:: RAMNTupleRefs::GetRefId

.. doxygenfunction:: RAMNTupleRefs::FindRefId

.. doxygenfunction:: RAMNTupleRefs::GetRefName

.. doxygenfunction:: RAMNTupleRefs::Size

.. doxygenfunction:: RAMNTupleRefs::GetRefs

.. doxygenfunction:: RAMNTupleRefs::AddRef

.. doxygenfunction:: RAMNTupleRefs::SetRefs

.. doxygenfunction:: RAMNTupleRefs::Clear

RAMCoordinateOrder
------------------

.. doxygenstruct:: RAMCoordinateOrder
   :members:
