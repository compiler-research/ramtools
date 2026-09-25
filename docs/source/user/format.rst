The RAM file format
===================

A RAM file is a ROOT file with two RNTuples and one key. You can open it
with any ROOT build that has RNTuple, and nothing in it is specific to
RAMTools except the record encoding described below.

======================  =====================================================
Object                  Contents
======================  =====================================================
``RAM`` (RNTuple)       one entry per alignment record, field ``record``
``METADATA`` (RNTuple)  one entry: name tables, longest span, sort flag
``headers`` (TList)     the SAM header, one ``TNamed`` per line
======================  =====================================================

The record
----------

``record`` is a ``RAMNTupleRecord``. RNTuple stores each member as its own
column, so a query that needs only ``refid``, ``pos`` and ``cigar`` reads
only those three.

.. list-table::
   :header-rows: 1
   :widths: 24 18 58

   * - Member
     - Type
     - Meaning
   * - ``qname``
     - string
     - read name, verbatim
   * - ``flag``
     - uint16
     - SAM FLAG
   * - ``refid``
     - int32
     - index into ``rname_refs``; -1 for ``*``
   * - ``pos``
     - int32
     - 0-based leftmost position; -1 if unplaced
   * - ``mapq``
     - uint8
     - mapping quality
   * - ``cigar``
     - vector<uint32>
     - packed CIGAR, see below
   * - ``refnext``
     - int32
     - index into ``rnext_refs``; -1 for ``*``
   * - ``pnext``
     - int32
     - 0-based mate position; -1 if none
   * - ``tlen``
     - int32
     - template length
   * - ``seq``
     - string
     - packed sequence, see below
   * - ``qual``
     - string
     - quality, encoding per ``compression_flags``
   * - ``tags``
     - vector<string>
     - optional fields verbatim, e.g. ``NM:i:0``
   * - ``compression_flags``
     - uint32
     - how ``qual`` is stored

Positions are stored 0-based and converted back to SAM's 1-based form on the
way out.

**Reference names.** ``refid`` and ``refnext`` index two string tables held
in ``METADATA``. ``rname_refs`` is filled from the ``@SQ`` header lines in
order, then extended by any name seen in a record. ``rnext_refs`` is a
separate table so that ``=`` (same reference as the read) can be stored as
its own entry and reproduced as ``=``.

**CIGAR.** Each operation is one ``uint32``: ``length << 4 | op``, with
``op`` numbered ``M I D N S H P = X`` from 0 to 8, the BAM encoding. The
length has 28 bits, so an operation longer than 268,435,455 bases is
rejected at conversion. An empty vector is SAM's ``*``.

**Sequence.** ``seq`` starts with a 4-byte little-endian length followed by
the bases packed two per byte, 4 bits each, in htslib's order
``=ACMGRSVTWYHKDBN`` (so ``A`` is 1, ``C`` 2, ``G`` 4, ``T`` 8, ``N`` 15).
Lower-case input is stored as the same code. Any other byte is stored as
``N``. An empty ``seq`` string is SAM's ``*``.

**Quality.** ``compression_flags`` says which of three encodings was used:

- ``kPhred33`` (bit 14): the QUAL column stored verbatim.
- ``kIlluminaBinning`` (bit 15): one byte per base holding the binned Phred
  value, 0 to 40, not the ASCII character. Reading adds 33. An empty string
  is ``*``.
- ``kDrop`` (bit 16): nothing stored; reading returns ``*``.

Metadata
--------

``METADATA`` has a single entry with four fields:

- ``rname_refs`` and ``rnext_refs``: the two name tables.
- ``max_ref_span``: the longest reference span of any record in the file.
  A region query starts its search that far before the region, so no record
  that begins before the region and reaches into it can be skipped. A value
  of 0 means the file predates this field, and queries start at the
  reference's first record.
- ``coordinate_sorted``: whether the records are in coordinate order. Files
  written before this field are read as sorted, which is what they always
  were assumed to be.

Finding a region
----------------

There is no separate index: in a sorted file the ``refid`` and ``pos``
columns are one. To answer ``rname:start-end`` a query binary-searches
those two columns for the first row at or after ``start - max_ref_span`` on
that reference, reads ``refid``, ``pos`` and ``cigar`` forward from there,
and stops as soon as ``refid`` changes or ``pos`` passes ``end``. Records in
between are tested with the overlap rule in :doc:`querying`. The search
takes about log2(rows) probes, 28 for 196 million records, each reading a
page of ``refid`` and ``pos``. When ``coordinate_sorted`` is false the query
skips the search and the early stop and tests every record instead.

Files written before the search replaced it also hold an ``INDEX`` ntuple
with a sparse position index. Current tools ignore it.

The header
----------

``headers`` is a ``TList`` written as a single key. Each element is a
``TNamed`` whose name is the record type (``@HD``, ``@SQ``, ``@PG``) and
whose title is everything after the first tab. ``ramdump -h`` joins them
back with a tab, in order. Lines with nothing after the tag (a bare
``@CO``) get no tab.

Reading a RAM file from your own code
-------------------------------------

The tuples are plain RNTuples, so this works from ROOT with only the
``ramcore`` library loaded for the record dictionary:

.. code-block:: cpp

   auto reader = ROOT::RNTupleReader::Open("RAM", "reads.ram");
   auto pos = reader->GetView<int32_t>("record.pos");
   auto cigar = reader->GetView<std::vector<uint32_t>>("record.cigar");
   for (auto i : reader->GetEntryRange()) {
      // pos(i), cigar(i)
   }

``RAMNTupleRecord::OpenRAMFile`` does the same and also loads the name
tables, the longest span and the sort flag, which the getters for names and
the region scan need.
