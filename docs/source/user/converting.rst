Converting to RAM
=================

Two tools write RAM files: ``samtoramntuple`` reads SAM text and
``bamtoramntuple`` reads BAM through htslib. They produce the same file
layout and take the same options, except that only ``samtoramntuple`` has
``-split``.

.. code-block:: bash

   samtoramntuple <input.sam> [output] [options]
   bamtoramntuple <input.bam> [output] [options]

If ``output`` is omitted the input name is used with its extension removed.
Unless the name already ends in ``.ram`` or ``.root``, ``.ram`` is appended.
A RAM file is an ordinary ROOT file; either extension works.

Options
-------

``-threads N``
   Convert on N threads; the default is 1. ``samtoramntuple`` reads the
   input in 64 MB blocks of whole lines, and every thread parses, encodes
   and compresses blocks of its own through ROOT's
   ``RNTupleParallelWriter``. The blocks are committed to the file in input
   order, so the records come out in the order they went in whatever the
   thread count; only the cluster layout differs. Memory use is about
   2 × N × 64 MB on top of ROOT's own. ``bamtoramntuple``, and
   ``samtoramntuple`` with ``-split``, parse on one thread and use the
   others to compress pages.

``-compression N``
   The ROOT compression code, ``algorithm * 100 + level``: 1 is ZLIB, 2
   LZMA, 4 LZ4, 5 ZSTD, with levels 1 to 9, and 0 means no compression. The
   default is 505, ZSTD level 5. Anything else is rejected with a message.
   Measured on HG00154 (196 million records, 72 GB of SAM; see
   :doc:`benchmarks`):

   ========  =========  =========  ============================================
   Code      Algorithm  Size (GB)  When to use it
   ========  =========  =========  ============================================
   ``505``   ZSTD 5     11.43      the default
   ``503``   ZSTD 3     11.74      half the CPU of the default, 3% larger
   ``501``   ZSTD 1     12.27      fastest to write; with ``-threads 4``
                                   quicker than ``samtools`` writes a BAM
   ``509``   ZSTD 9     10.38      smallest, about ten times the CPU of 505
   ``404``   LZ4 4      16.54      larger than BAM (15.22 GB)
   ``101``   ZLIB 1     14.14      for readers built without ZSTD
   ``0``     none       83.71      for measuring the cost of compression
   ========  =========  =========  ============================================

   The codec makes little difference to query time.

``-illumina``
   Store quality scores with Illumina's 8-level binning instead of the full
   Phred range. Every score maps to one of 0, 1, 6, 15, 22, 27, 33, 37 or 40,
   which compresses far better and is what most variant callers are
   evaluated against. The mapping is lossy.

``-dropqual``
   Store no quality scores at all. Reads come back with ``*`` in the QUAL
   column.

``-split``
   ``samtoramntuple`` only. Write one file per reference sequence instead of
   one file, named ``<output>_<rname>.root``. Records with no reference
   (``*``) are not written. See :ref:`split`.

What the input has to look like
-------------------------------

**Sort it for fast queries.** A region query on a sorted file finds the
region by binary search and stops at the first record past it. Both steps
assume the records are in coordinate order, the order ``samtools sort``
produces. The conversion checks the order as it goes and records the answer
in the file. Unsorted input is stored as it arrives and marked unsorted; the
conversion says so on standard error. Region queries on such a file are
still exact, but read every record. For fast queries, sort first:

.. code-block:: bash

   samtools sort -o sorted.bam reads.bam
   bamtoramntuple sorted.bam reads.ram

**Malformed records are skipped, not stored.** The SAM parser validates every
mandatory column: integer fields must be integers in range, the CIGAR must be
``*`` or a run of ``<length><op>`` pairs with the operators ``MIDNSHP=X``,
and no mandatory column may be empty. A record that fails is reported with
its line number on standard error and left out. The counts printed at the
end tell you how many records were written.

**Everything else round-trips.** Read names, FLAG, positions, MAPQ, CIGAR,
mate fields, TLEN, sequence, quality and every optional tag are stored so
that ``ramdump -h`` reproduces the input line for line. Lower-case bases are
stored as their upper-case code; a base outside ``ACGTN`` and the IUPAC
codes is stored as ``N``.

.. _split:

Splitting by chromosome
-----------------------

``samtoramntuple reads.sam out -split`` writes ``out_chr1.root``,
``out_chr2.root`` and so on, one complete RAM file per reference sequence,
each with its own metadata and a copy of the header. Records stream to
their file as they are parsed, in the order they arrive, so the split needs
no more memory than a plain conversion. Each file records whether its own
records are in coordinate order. Query and dump them like any other RAM
file. Records with no reference sequence are not written to any of them.

What gets printed
-----------------

Both tools print the name of the file written, the number of records in it
and the reference table. ``samtoramntuple`` also prints how many header
lines and records it read, so the difference from the number written is the
number of skipped records.
