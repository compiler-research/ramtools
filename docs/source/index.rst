RAMTools
========

RAMTools stores genomic alignments in RAM, a columnar file format built on
ROOT's `RNTuple <https://root.cern/doc/master/group__NTuple.html>`_. A RAM
file holds the same records as a SAM or BAM file, but each field lives in its
own compressed column, so a region query reads only the columns it needs and
a whole-file scan streams at disk speed.

The tools do four things:

- convert SAM or BAM to RAM (``samtoramntuple``, ``bamtoramntuple``)
- count the records overlapping a region (``ramntupleview``)
- write a RAM file, or a region of it, back out as SAM, as ``samtools view``
  would (``ramdump``)
- split an input into one RAM file per chromosome (``samtoramntuple -split``)

Five commands
-------------

.. code-block:: bash

   # convert on four threads
   samtoramntuple reads.sam reads.ram -threads 4

   # count the records overlapping a region (1-based, inclusive)
   ramntupleview reads.ram chr1:10000-20000

   # dump a region as SAM, header included
   ramdump -h reads.ram chr1:10000-20000 > region.sam

   # the dump of the whole file reproduces the input
   ramdump -h reads.ram | diff reads.sam -

   # and its counts agree with samtools
   ramdump -c reads.ram chr1:10000-20000
   samtools view -c reads.bam chr1:10000-20000

Where to go next
----------------

- :doc:`user/installation` builds the tools and runs the tests.
- :doc:`user/converting` covers every conversion flag, compression codes, and
  what the input has to look like for region queries to be fast.
- :doc:`user/querying` explains region syntax, what "overlapping" means, and
  how to check any result against samtools.
- :doc:`user/format` describes what is inside a RAM file. Read it if you want
  to understand why queries are fast or to read the files from your own code.
- :doc:`user/benchmarks` runs the benchmark binaries and gives the HG00154
  numbers against BAM and CRAM.
- :doc:`reference/index` is the C++ API, generated from the headers, for
  using ``ramcore`` from your own code.
- :doc:`dev/architecture` explains how the code is put together: the
  layout, how a record gets in, how a query gets out, and the shared state.
- :doc:`dev/contributing` is about getting a change in: tests, style, CI,
  pull requests.

.. toctree::
   :maxdepth: 2
   :hidden:

   user/installation
   user/converting
   user/querying
   user/format
   user/benchmarks
   reference/index
   dev/architecture
   dev/contributing
