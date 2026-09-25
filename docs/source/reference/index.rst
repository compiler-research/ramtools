API reference
=============

Everything below is in the ``ramcore`` shared library and declared in the
headers under ``inc/``. The command-line tools are thin wrappers over these
functions, so anything a tool does can be done from C++.

.. toctree::
   :maxdepth: 1

   record
   tables
   parser
   converters
   query

Before you start
----------------

**Shared state.** The reference name tables, the longest span and the sort
flag are static members of ``RAMNTupleRecord``, one set per process. ``InitializeRefs()`` resets the per-file parts; only the
writers and ``OpenRAMFile()`` call it, and constructing a record never
does. Two files cannot be open for querying at the same time in one
process. The rule and the bug it prevents are described in
:doc:`../dev/architecture`.

**Values in, values out.** Setters take SAM values (1-based positions, text
CIGAR, plain bases, Phred+33 quality) and store the encoded form. Getters
decode. You only meet the packed representation if you read the columns
through RNTuple directly, and :doc:`../user/format` describes it.

A conversion and a query
------------------------

.. code-block:: cpp

   #include "ramcore/SamToNTuple.h"
   #include "ramcore/RAMNTupleView.h"
   #include "rntuple/RAMNTupleRecord.h"

   // SAM -> RAM, ZSTD level 5, quality kept verbatim, on four threads
   if (!samtoramntuple("reads.sam", "reads.ram", 505, RAMNTupleRecord::kPhred33, /*threads=*/4))
      return 1;

   // open: loads the name tables, the longest span and the sort flag
   auto reader = RAMNTupleRecord::OpenRAMFile("reads.ram");

   // walk a region; the callback receives each matching row
   auto view = reader->GetView<RAMNTupleRecord>("record");
   ramntuplescan(*reader, "chr1:10000-20000", [&](Long64_t row) {
      const RAMNTupleRecord &rec = view(row);
      // rec.GetQNAME(), rec.GetPOS(), rec.GetCIGAR(), rec.GetSEQ() ...
   });
