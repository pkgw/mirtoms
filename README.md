mirtoms
=======

This program converts MIRIAD datasets to MeasurementSet datasets. It is
heavily derived from Peter Teuben's `carmafiller`, which in turn came from
`bimafiller`, which in turn came from `uvfitsfiller`.

This program is customized to work with simple ATA datasets. Key differences
with `carmafiller` are:

  * it is [less chatty](http://newton.cx/~peter/2012/12/silence-is-golden/);
  * it includes appropriate hardcoded info for the ATA;
  * groups of related polarization records are written properly, so
    CASA Stokes processing can work;
  * the starting SCAN_NUM can be set on the command-line, so that multiple
    generated MSes can be concatenated without warnings;
  * extra information is included to allow reverse-mapping from the MS
    back to MIRIAD data.

These changes make it so that:

  * you can export a bunch of MIRIAD datasets, concatenate them, and image
    them with the CASA imager, with propert Stokes processing;
  * you can run André Offringa’s [aoflagger](http://aoflagger.sourceforge.net)
    on your MIRIAD datasets and then import the flags back into MIRIAD,
    using the included `mirmsflagextract` tool.

The `mirtoms` tool has severe limitations and will only work with very simple
MIRIAD datasets. It also probably gets various details wrong that will bite
you in the millimeter regime, but not the centimeter regime where I work.
`carmafiller` is more robust in both ways.


Installation
============

The `Makefile` is customized to my machine. Sorry. Compilation is
extremely straightforward if you can get the simultaneous linking with both
MIRIAD and CASA libraries to work. If you're lucky, you should be able to
build with just:

    make MIR=/path/to/miriad CASACORE=/path/to/casacore


License and Copyright Information
=================================

The code is licensed under the GNU GPL, version 2 or later.

This file is copyright 2013 Peter Williams. It is free documentation; the
copyright holder gives unlimited permission to copy, distribute, and modify
it.
