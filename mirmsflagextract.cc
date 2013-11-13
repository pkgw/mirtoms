/* mirmsflagextract: copy flags from a mirtoms MS back into the MIRIAD dataset
   Copyright 2013 Peter Williams
   Licensed under the GNU GPL version 2 or later.

   Derived from mirtoms.cc, which from Peter Teuben's carmafiller, which in
   turn came from bimafiller, and before that, uvfitsfiller.
*/

/*
  A *very* incomplete list of deficiencies and TODOs:
  - we ignore "wide" channels
  - option to 'or' in flags instead of straight copying
*/

#include <casa/aips.h>
#include <casa/stdio.h>
#include <casa/iostream.h>
#include <casa/OS/File.h>
#include <casa/Utilities/GenSort.h>
#include <casa/Arrays/Cube.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/ArrayUtil.h>
#include <casa/Arrays/ArrayLogical.h>
#include <casa/Arrays/MatrixMath.h>
#include <casa/Inputs/Input.h>
#include <casa/namespace.h>

#include <measures/Measures.h>
#include <measures/Measures/MPosition.h>
#include <measures/Measures/MeasData.h>
#include <measures/Measures/Stokes.h>
#include <tables/Tables.h>
#include <tables/Tables/TableInfo.h>
#include <ms/MeasurementSets.h>

#include <miriad-c/maxdimc.h>
#include <miriad-c/miriad.h>


#define MYMAXCHAN 8192 // I feel so dirty.
#define WARN(message) (cerr << "warning: " << message << endl);

const char *MIR_REC_COL = "MIRIAD_RECNUM";


void
extract_flags (String& mspath, String& vispath)
{
    if (sizeof (double) != sizeof (Double))
	WARN ("sizeof(Double) != sizeof(double); mirmsflagextract will probably fail");
    if (sizeof (int) != sizeof (Int))
	WARN ("sizeof(Int) != sizeof(int); mirmsflagextract will probably fail");

    // Open the MIRIAD dataset.

    int mirhandle;

    uvopen_c (&mirhandle, vispath.chars (), "old");
    uvset_c (mirhandle, "preamble", "uvw/time/baseline", 0, 0.0, 0.0, 0.0);

    // Open the CASA dataset.

    MeasurementSet ms (mspath, Table::Old);
    MSColumns msc (ms);
    ScalarColumn<Int> mirreccol (ms, MIR_REC_COL);

    /* Start charging through. CASA stores multiple polarization records in one
       logical row, while MIRIAD separates out the records. So we read a CASA row
       first, save the data, and apply it to 1-4 MIRIAD records. Our iteration is
       all driven by the MIRIAD dataset, though.
    */

    Int recnum = 0, row = 0;
    int polsleft = 0;
    double preamble[5];
    float data[2 * MYMAXCHAN]; // complex, so 2 floats per channel
    int flags[MYMAXCHAN];

    while (1) {
	/* As far as I know, we need to actually read the UV data and friends,
	   though we don't actually use them for anything... */
	Int nread;
	uvread_c (mirhandle, preamble, data, flags, MYMAXCHAN, &nread);
	if (nread <= 0)
	    break;

	if (polsleft == 0) {
	    // We just started a new simultaneous polarization record.
	    uvrdvr_c (mirhandle, H_INT, "npol", (char *) &polsleft, NULL, 1);

	    Int msrecnum = mirreccol.get (row);
	    if (msrecnum != recnum)
		throw AipsError ("synchronization failure at MIRIAD record #" +
				 String::toString (recnum) +
				 " (0-based) and CASA row #" +
				 String::toString (row) +
				 " (0-based): should be at MIRIAD record #" +
				 String::toString (msrecnum));

	    row++; // next time, next row.
	}

	int mirpol;
	uvrdvr_c (mirhandle, H_INT, "pol", (char *) &mirpol, NULL, 1);

	recnum++;
	polsleft--;
    }

    // Wrap up.

    uvclose_c (mirhandle);
}


int
main (int argc, char **argv)
{
    try {
	Input inp (1);
	inp.version ("");
	inp.create ("vis", "", "path of MIRIAD dataset to modify", "string");
	inp.create ("ms", "", "path of MeasurementSet dataset with flags", "string");
	inp.readArguments (argc, argv);

	String vis (inp.getString ("vis"));
	if (vis == "")
	    throw AipsError ("no MIRIAD input path (vis=) given");
	if (! File (vis).isDirectory ())
	    throw AipsError ("MIRIAD input path (vis=) does not refer to a directory");

	String ms (inp.getString ("ms"));
	if (ms == "")
	    throw AipsError ("no MS input path (ms=) given");
	if (! File (ms).isDirectory ())
	    throw AipsError ("MS input path (ms=) does not refer to a directory");

	extract_flags (ms, vis);
    } catch (AipsError x) {
	cerr << "error: " << x.getMesg () << endl;
	return 1;
    }

    return 0;
}
