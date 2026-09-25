Querying and dumping
====================

``ramntupleview`` counts, ``ramdump`` prints. Both take the same region
syntax and use the same overlap rule, and both agree with ``samtools view``
on what a region contains, with the one exception noted below.

Regions
-------

A region is ``rname``, ``rname:pos`` or ``rname:start-end``. Positions are
1-based and inclusive, as in samtools. ``rname`` alone means the whole
reference; no region at all, or ``*``, means the whole file.

``rname:pos`` is the one form that differs from samtools: here it is the
single base ``pos``, while ``samtools view`` reads it as ``pos`` to the end
of the reference. Write ``rname:pos-end`` when you compare the two.

Reference names may contain colons (GRCh38 has ``HLA-A*01:01:01:01``). The
tools first try the whole string as a name, then split on the last colon,
which is the order samtools uses.

What "overlapping" means
------------------------

A record is in a region when it is placed on that reference and the span it
covers, from POS to POS plus the reference length of its CIGAR minus one,
meets the region. So a 100-base read starting 50 bases before the region is
included, a read whose CIGAR is ``50S50M`` covers 50 reference bases, not
100, and a placed record with no CIGAR counts as covering the one base it
sits on. That is htslib's ``bam_endpos`` rule.

Secondary and supplementary alignments are included, as ``samtools view``
includes them. Unmapped records are included if they carry a position (a
mate placed next to its partner) and excluded if they do not. Use ``-F`` on
``ramdump`` to filter them out, as you would with samtools.

ramntupleview
-------------

.. code-block:: bash

   ramntupleview <file.ram> [region]

Prints the time the query took and the count. On a sorted file the query
finds its first row by binary search over the ``refid`` and ``pos`` columns,
starting the longest reference span in the file before the region so that
nothing that starts earlier and reaches in is missed, and stops at the first
record past it. On an unsorted file it tests every record.

ramdump
-------

.. code-block:: bash

   ramdump [options] <file.ram> [region]

=================  ==========================================================
``-h``             include the header
``-H``             print the header only
``-c``             print only the number of matching records
``-f INT``         keep records with all of these FLAG bits set
``-F INT``         drop records with any of these FLAG bits set
``-o FILE``        write to FILE instead of standard output
=================  ==========================================================

FLAG values may be decimal or hexadecimal (``0x900``). Records are printed in
file order, exactly as they went in: the same columns, the same tags in the
same order, the same header lines.

Checking against samtools
-------------------------

Every result here can be verified. For a coordinate-sorted BAM with an
index:

.. code-block:: bash

   # counts
   ramdump -c reads.ram chr1:1-1000000
   samtools view -c reads.bam chr1:1-1000000

   # primary alignments only
   ramdump -c -F 0x904 reads.ram chr1
   samtools view -c -F 0x904 reads.bam chr1

   # the records themselves
   diff <(ramdump reads.ram chr1) <(samtools view reads.bam chr1)

   # the whole file, header included
   diff <(ramdump -h reads.ram) <(samtools view -h reads.bam)

The only expected difference is a ``@PG`` line that samtools appends to
record its own invocation.

Files converted before the header was stored as a single key print records
but no header. Convert them again.
