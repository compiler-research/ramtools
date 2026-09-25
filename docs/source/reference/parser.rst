SAM parser
==========

``ramcore::SamParser`` is the only path from SAM text into the writers. It
validates every mandatory column, so a converter never sees a malformed
line. ``ParseRecord()`` parses one line and is what the threaded converter's
workers call; ``ParseFile()`` reads a whole file and delivers records
through callbacks, which the split converter uses. The BAM path does not use it;
``bamtoramntuple`` formats htslib's fields to the same values instead.

SamRecord
---------

.. doxygenstruct:: ramcore::SamRecord
   :members:
   :undoc-members:

SamParser
---------

.. doxygenclass:: ramcore::SamParser

.. doxygenfunction:: ramcore::SamParser::ParseRecord

.. doxygenfunction:: ramcore::SamParser::ParseFile

.. doxygenfunction:: ramcore::SamParser::GetLinesProcessed

.. doxygenfunction:: ramcore::SamParser::GetRecordsProcessed

Helpers
-------

.. doxygenfunction:: ramcore::StripCRLF
