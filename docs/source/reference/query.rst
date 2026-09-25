Region scan
===========

``ramntuplescan`` is the one implementation of "which records overlap this
region". ``ramntupleview`` and ``ramdump`` are both built on it, and code
that needs the records of a region should call it rather than search the
columns by hand. The region syntax and the overlap rule are described in
:doc:`../user/querying`.

ramntuplescan
-------------

.. doxygenfunction:: ramntuplescan

ramntupleview
-------------

.. doxygenfunction:: ramntupleview

.. doxygenstruct:: RAMNTupleViewOpts
   :members:
   :undoc-members:
