Record and encoders
===================

RAMNTupleRecord
---------------

.. doxygenclass:: RAMNTupleRecord

Columns
~~~~~~~

The public data members are the columns of the ``RAM`` ntuple. Read them
directly only when you need the raw encoding; the getters below decode.

.. list-table::
   :header-rows: 1
   :widths: 26 20 54

   * - Member
     - Type
     - Holds
   * - ``qname``
     - ``std::string``
     - read name, verbatim
   * - ``flag``
     - ``uint16_t``
     - SAM FLAG
   * - ``refid``
     - ``int32_t``
     - index into the RNAME table; -1 for ``*``
   * - ``pos``
     - ``int32_t``
     - 0-based leftmost position; -1 when unplaced
   * - ``mapq``
     - ``uint8_t``
     - mapping quality
   * - ``cigar``
     - ``std::vector<uint32_t>``
     - operations packed as ``(length << 4) | op``
   * - ``refnext``
     - ``int32_t``
     - index into the RNEXT table; -1 for ``*``
   * - ``pnext``
     - ``int32_t``
     - 0-based mate position; -1 when none
   * - ``tlen``
     - ``int32_t``
     - template length
   * - ``seq``
     - ``std::string``
     - 4-byte little-endian length, then bases packed two per byte
   * - ``qual``
     - ``std::string``
     - quality, encoded per ``compression_flags``
   * - ``tags``
     - ``std::vector<std::string>``
     - optional fields verbatim, e.g. ``NM:i:0``
   * - ``compression_flags``
     - ``uint32_t``
     - one of ``EQualCompressionBits``

Quality policy
~~~~~~~~~~~~~~

.. doxygenenum:: RAMNTupleRecord::EQualCompressionBits

.. doxygenfunction:: RAMNTupleRecord::SetCompressionMode

.. doxygenfunction:: RAMNTupleRecord::SetBit

.. doxygenfunction:: RAMNTupleRecord::TestBit

Setters
~~~~~~~

All take SAM values: 1-based positions, a text CIGAR, plain bases, Phred+33
quality. ``SetREFID`` and ``SetREFNEXT`` are the names the converters use
for ``SetRNAME`` and ``SetRNEXT``; ``SetOPT`` and ``ResetNOPT`` are the
same as ``AddTag`` and ``ClearTags``.

.. doxygenfunction:: RAMNTupleRecord::SetQNAME

.. doxygenfunction:: RAMNTupleRecord::SetFLAG

.. doxygenfunction:: RAMNTupleRecord::SetRNAME

.. doxygenfunction:: RAMNTupleRecord::SetPOS

.. doxygenfunction:: RAMNTupleRecord::SetMAPQ

.. doxygenfunction:: RAMNTupleRecord::SetCIGAR

.. doxygenfunction:: RAMNTupleRecord::SetRNEXT

.. doxygenfunction:: RAMNTupleRecord::SetPNEXT

.. doxygenfunction:: RAMNTupleRecord::SetTLEN

.. doxygenfunction:: RAMNTupleRecord::SetSEQ

.. doxygenfunction:: RAMNTupleRecord::SetQUAL

.. doxygenfunction:: RAMNTupleRecord::AddTag

.. doxygenfunction:: RAMNTupleRecord::ClearTags

Getters
~~~~~~~

All give SAM values back. ``GetRNAME`` and ``GetRNEXT`` resolve the ids
through the shared name tables, so they need the file's metadata loaded
(``OpenRAMFile`` does that).

.. doxygenfunction:: RAMNTupleRecord::GetQNAME

.. doxygenfunction:: RAMNTupleRecord::GetFLAG

.. doxygenfunction:: RAMNTupleRecord::GetRNAME

.. doxygenfunction:: RAMNTupleRecord::GetREFID

.. doxygenfunction:: RAMNTupleRecord::GetPOS

.. doxygenfunction:: RAMNTupleRecord::GetMAPQ

.. doxygenfunction:: RAMNTupleRecord::GetCIGAR

.. doxygenfunction:: RAMNTupleRecord::GetRNEXT

.. doxygenfunction:: RAMNTupleRecord::GetREFNEXT

.. doxygenfunction:: RAMNTupleRecord::GetPNEXT

.. doxygenfunction:: RAMNTupleRecord::GetTLEN

.. doxygenfunction:: RAMNTupleRecord::GetSEQ

.. doxygenfunction:: RAMNTupleRecord::GetQUAL

.. doxygenfunction:: RAMNTupleRecord::GetTags

.. doxygenfunction:: RAMNTupleRecord::GetSEQLEN

.. doxygenfunction:: RAMNTupleRecord::GetRefSpan

.. doxygenfunction:: RAMNTupleRecord::GetNCIGAROP

.. doxygenfunction:: RAMNTupleRecord::GetCIGAROPLEN

.. doxygenfunction:: RAMNTupleRecord::GetCIGAROP

.. doxygenfunction:: RAMNTupleRecord::Print

Shared state
~~~~~~~~~~~~

.. doxygenfunction:: RAMNTupleRecord::InitializeRefs

.. doxygenfunction:: RAMNTupleRecord::GetRnameRefs

.. doxygenfunction:: RAMNTupleRecord::GetRnextRefs

.. doxygenfunction:: RAMNTupleRecord::GetMaxRefSpan

.. doxygenfunction:: RAMNTupleRecord::NoteRefSpan

.. doxygenfunction:: RAMNTupleRecord::IsCoordinateSorted

.. doxygenfunction:: RAMNTupleRecord::SetCoordinateSorted

.. doxygenfunction:: RAMNTupleRecord::NotePlacement

Files
~~~~~

.. doxygenfunction:: RAMNTupleRecord::OpenRAMFile

.. doxygenfunction:: RAMNTupleRecord::MakeModel

.. doxygenfunction:: RAMNTupleRecord::WriteAllRefs

.. doxygenfunction:: RAMNTupleRecord::ReadAllRefs

Encoders
--------

The functions behind the setters and getters, in namespace
``RAMNTupleUtils``. Use them when you read columns directly.

.. doxygenfunction:: RAMNTupleUtils::EncodeSequence

.. doxygenfunction:: RAMNTupleUtils::DecodeSequence

.. doxygenfunction:: RAMNTupleUtils::EncodeQuality

.. doxygenfunction:: RAMNTupleUtils::DecodeQuality

.. doxygenfunction:: RAMNTupleUtils::ParseCIGAR

.. doxygenfunction:: RAMNTupleUtils::FormatCIGAR

CIGAR operation codes
---------------------

``ramcore/CigarOps.h`` defines the codes used inside a packed CIGAR, in
BAM's order.

.. list-table::
   :header-rows: 1
   :widths: 30 10 60

   * - Constant
     - Code
     - Operation
   * - ``RAM_CIGAR_M``
     - 0
     - alignment match, consumes reference
   * - ``RAM_CIGAR_I``
     - 1
     - insertion
   * - ``RAM_CIGAR_D``
     - 2
     - deletion, consumes reference
   * - ``RAM_CIGAR_N``
     - 3
     - skipped region, consumes reference
   * - ``RAM_CIGAR_S``
     - 4
     - soft clip
   * - ``RAM_CIGAR_H``
     - 5
     - hard clip
   * - ``RAM_CIGAR_P``
     - 6
     - padding
   * - ``RAM_CIGAR_EQUAL``
     - 7
     - sequence match, consumes reference
   * - ``RAM_CIGAR_X``
     - 8
     - sequence mismatch, consumes reference
