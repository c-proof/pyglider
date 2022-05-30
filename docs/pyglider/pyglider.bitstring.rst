======================
``pyglider.bitstring``
======================

.. automodule:: pyglider.bitstring

   .. contents::
      :local:

.. currentmodule:: pyglider.bitstring


Functions
=========

- :py:func:`pack`:
  Pack the values according to the format string and return a new BitStream.


.. autofunction:: pack


Classes
=======

- :py:class:`ConstBitArray`:
  A container holding an immutable sequence of bits.

- :py:class:`ConstBitStream`:
  A container or stream holding an immutable sequence of bits.

- :py:class:`BitStream`:
  A container or stream holding a mutable sequence of bits

- :py:class:`BitArray`:
  A container holding a mutable sequence of bits.

- :py:class:`Bits`:
  A container holding an immutable sequence of bits.

- :py:class:`BitString`:
  A container or stream holding a mutable sequence of bits


.. autoclass:: ConstBitArray
   :members:

   .. rubric:: Inheritance
   .. inheritance-diagram:: ConstBitArray
      :parts: 1

.. autoclass:: ConstBitStream
   :members:

   .. rubric:: Inheritance
   .. inheritance-diagram:: ConstBitStream
      :parts: 1

.. autoclass:: BitStream
   :members:

   .. rubric:: Inheritance
   .. inheritance-diagram:: BitStream
      :parts: 1

.. autoclass:: BitArray
   :members:

   .. rubric:: Inheritance
   .. inheritance-diagram:: BitArray
      :parts: 1

.. autoclass:: Bits
   :members:

   .. rubric:: Inheritance
   .. inheritance-diagram:: Bits
      :parts: 1

.. autoclass:: BitString
   :members:

   .. rubric:: Inheritance
   .. inheritance-diagram:: BitString
      :parts: 1


Exceptions
==========

- :py:exc:`Error`:
  Base class for errors in the bitstring module.

- :py:exc:`ReadError`:
  Reading or peeking past the end of a bitstring.

- :py:exc:`InterpretError`:
  Inappropriate interpretation of binary data.

- :py:exc:`ByteAlignError`:
  Whole-byte position or length needed.

- :py:exc:`CreationError`:
  Inappropriate argument during bitstring creation.


.. autoexception:: Error

   .. rubric:: Inheritance
   .. inheritance-diagram:: Error
      :parts: 1

.. autoexception:: ReadError

   .. rubric:: Inheritance
   .. inheritance-diagram:: ReadError
      :parts: 1

.. autoexception:: InterpretError

   .. rubric:: Inheritance
   .. inheritance-diagram:: InterpretError
      :parts: 1

.. autoexception:: ByteAlignError

   .. rubric:: Inheritance
   .. inheritance-diagram:: ByteAlignError
      :parts: 1

.. autoexception:: CreationError

   .. rubric:: Inheritance
   .. inheritance-diagram:: CreationError
      :parts: 1


Variables
=========

- :py:data:`bytealigned`

.. autodata:: bytealigned
   :annotation:

   .. code-block:: text

      False
