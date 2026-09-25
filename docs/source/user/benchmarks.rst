Benchmarks
==========

Four Google Benchmark binaries are built with ``RAMTOOLS_BUILD_BENCHMARKS``
and live in ``build/benchmark``:

======================================  ======================================
``sam_to_ram_benchmark``                conversion throughput on generated SAM
``conversion_time_benchmark``           conversion time by input size
``chromosome_split_benchmark``          ``-split`` against ``samtools`` splitting
``region_query_benchmark``              region queries on a file you provide
======================================  ======================================

The first three generate their own input. The region query benchmark
measures a real file, so it takes it from the environment:

.. code-block:: bash

   export RAMTOOLS_BENCH_RNTUPLE=/data/HG00154.ram
   export RAMTOOLS_BENCH_REGIONS="chr1:1000000-1001000,chr1:1-50000000,chr21"
   build/benchmark/region_query_benchmark

``RAMTOOLS_BENCH_REGIONS`` is optional; without it a fixed set of ``chr1``
and ``chr21`` regions is used. Each region is one benchmark case, and the
label of each case shows the region and how many records it returned, so a
zero is easy to spot. The binary refuses to run without
``RAMTOOLS_BENCH_RNTUPLE`` rather than measure an open that failed.

``scripts/run_benchmark.py`` runs the set and ``scripts/render_benchmark.py``
turns the JSON output into Markdown tables.

Reproducing a number
--------------------

A benchmark result is only a result if someone else can get it. Record
alongside every number: the input file and how it was made (``samtools
sort`` order or not), the compression code, the thread count, the ROOT
version, and the machine.

HG00154 against BAM and CRAM
----------------------------

The numbers in the README come from the HG00154 sample of the 1000 Genomes
Project: 196,040,370 records, 72.06 GB of SAM, coordinate sorted, GRCh37
reference names. They were measured with ``/usr/bin/time`` and the tools
themselves rather than with the binaries above, on an Intel i3-7020U (2
cores, 4 threads, 12 GB of memory) with input and output on one hard disk,
against samtools 1.13 on 4 threads.

=============================================  ===========  =========  =========
Output                                         Wall         CPU (s)    Size (GB)
=============================================  ===========  =========  =========
BAM, ``samtools view -@ 4 -b``                 24:08        3,589      15.22
CRAM, ``samtools view -@ 4 -C -T ref.fasta``   12:34        1,817      7.77
RAM, ZSTD 5, one thread                        1:30:06      5,258      11.43
RAM, ZSTD 5, ``-threads 4``                    32:59        7,664      11.43
RAM, ZSTD 3, ``-threads 4``                    17:09        3,513      11.74
RAM, ZSTD 1, ``-threads 4``                    14:32        1,273      12.27
=============================================  ===========  =========  =========

CPU is user plus system time. A CRAM also needs its reference FASTA, 3.15 GB,
to be read. Every output holds the same records.

Region queries, best of three with a warm cache, counted with ``samtools
view -c`` and ``ramdump -c``:

================================  =========  =======  ========  =======
Region                            Records    BAM (s)  CRAM (s)  RAM (s)
================================  =========  =======  ========  =======
``1:1000000-1000100``             2          0.15     0.10      0.48
``13:32889611-32973805`` (BRCA2)  6,057      0.10     0.09      0.39
``1:10000000-20000000``           582,985    0.80     1.75      0.40
``2:1-100000000``                 7,048,385  8.59     21.43     0.85
================================  =========  =======  ========  =======

About 0.4 s of every RAM query is starting ROOT and opening the file; after
that a count reads about 15 million records per second, since it needs only
the ``refid``, ``pos`` and ``cigar`` columns.
