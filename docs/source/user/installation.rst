Installation
============

Requirements
------------

- `ROOT <https://root.cern>`_ 6.38 or newer, built with RNTuple (every binary
  release is). ``find_package(ROOT)`` needs ``thisroot.sh`` sourced or
  ``ROOT_DIR`` set.
- `htslib <https://github.com/samtools/htslib>`_, found through
  ``pkg-config``. On Debian and Ubuntu that is ``libhts-dev``.
- CMake 3.16 or newer and a C++17 compiler.
- ``samtools`` is not required, but it is the reference the tests and this
  documentation compare against.

Build
-----

.. code-block:: bash

   git clone https://github.com/compiler-research/ramtools
   cd ramtools
   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
   cmake --build build -j

The tools end up in ``build/tools``. CMake options:

.. list-table::
   :header-rows: 1
   :widths: 30 10 60

   * - Option
     - Default
     - Effect
   * - ``RAMTOOLS_BUILD_TOOLS``
     - ON
     - the four command-line tools
   * - ``RAMTOOLS_BUILD_TESTS``
     - ON
     - googletest binaries and ``ctest`` entries
   * - ``RAMTOOLS_BUILD_BENCHMARKS``
     - ON
     - Google Benchmark binaries (see :doc:`benchmarks`)
   * - ``RAMTOOLS_BUILD_DOCS``
     - OFF
     - a ``docs`` target that builds this site
   * - ``ENABLE_COVERAGE``
     - OFF
     - compile with coverage instrumentation

googletest and Google Benchmark are fetched at configure time.

Run the tests
-------------

.. code-block:: bash

   cd build && ctest --output-on-failure

The suite converts generated SAM files, queries them, dumps them back, and
compares counts and text with what samtools would produce. A run takes about
a minute.

Install
-------

``cmake --install build`` puts the tools under ``bin`` of the install prefix.
