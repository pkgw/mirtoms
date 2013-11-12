/* pwcarmafiller: hacked MIRIAD to MeasurementSet conversion
   adapted from Peter Teuben's carmafiller by Peter Williams.

   Copyright 2010-2013 Peter Williams, Peter Teuben
   Earlier versions copyright 1997, 2000, 2001, 2002
     Associated Universities, Inc. Washington DC, USA.

   Licensed under the GNU GPL version 2 or later.
*/

/*
Based off of carmafiller, which came from bimafiller, and before that,
uvfitsfiller.

The big problem with this program right now is that we really need to make two
passes through the dataset to build up lists of all of the pol'n configs,
pointings, etc., before we can actually start writing everything out.
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


#ifndef MAXFIELD
# define MAXFIELD 256 // TODO: kill this hardcoding.
#endif

#define WARN(message) (cerr << "warning: " << message << endl);


typedef struct window {
    // CASA defines everything mid-band, mid-interval
    int    nspect;                   // number of valid windows (<=MAXWIN, typically 16)
    int    nschan[MAXWIN+MAXWIDE];   // number of channels in a window
    int    ischan[MAXWIN+MAXWIDE];   // starting channel of a window (1-based)
    double sdf[MAXWIN+MAXWIDE];      // channel separation
    double sfreq[MAXWIN+MAXWIDE];    // frequency of first channel in window (doppler changes)
    double restfreq[MAXWIN+MAXWIDE]; // rest freq, if appropriate
    char   code[MAXWIN+MAXWIDE];     // code to CASA identification (N, W or S; S not used anymore)

    int    nwide;                    // number of wide band channels
    float  wfreq[MAXWIDE];           // freq
    float  wwidth[MAXWIDE];          // width
} WINDOW;


class CarmaFiller {
public:
    CarmaFiller (String& infile, Int debug_level=0, Bool apply_tsys=False);

    void checkInput ();
    void setupMeasurementSet (const String& ms_path);
    void fillObsTables ();
    void fillMSMainTable (Int snumbase);
    void fillAntennaTable ();
    void fillSyscalTable ();
    void fillSpectralWindowTable ();
    void fillFieldTable ();
    void fillSourceTable ();
    void fillFeedTable ();
    void fixEpochReferences ();

private:
    void setup_tracking ();
    void track_updates ();
    void init_window_info ();

    bool uv_hasvar (const char *varname);
    char *uv_getstr (const char *varname);
    int uv_getint (const char *varname);
    float uv_getfloat (const char *varname);
    double uv_getdouble (const char *varname);
    void uv_getfloats (const char *varname, float *dest, int count);
    void uv_getdoubles (const char *varname, double *dest, int count);

    String infile_p;
    Int uv_handle_p;
    MeasurementSet ms_p;
    MSColumns *msc_p;
    Int debug_level;
    String telescope_name, project_name, object_p, telescope_p,
	observer_name;
    Vector<Int> nPixel_p, corrType_p, corrIndex_p;
    Matrix<Int> corrProduct_p;
    Double epoch_p;
    MDirection::Types epochRef_p;
    Int num_arrays;
    Block<Int> nAnt_p;
    Block<Vector<Double> > receptorAngle_p;
    Vector<Double> arrayXYZ_p; // 3 elements
    Vector<Double> ras_p, decs_p; // ra/dec for source list (source_p)
    Vector<String> source_p, purpose_p; // list of source names (?? object_p ??)

    // the following variables are for miriad, hence not Double/Int/Float

    double preamble[5], first_time;
    int ifield, nfield, npoint, nsource;     // both dra/ddec should become Vector's
    float dra[MAXFIELD], ddec[MAXFIELD];       // offset in radians
    double ra[MAXFIELD], dec[MAXFIELD];
    int field[MAXFIELD];                     // source index
    int sid_p[MAXFIELD];
    float dra_p, ddec_p;
    int pol_p;

    Vector<Int> polmapping;
    Int nants_p, nchan_p, nwide_p, npol_p;
    Double antpos[3*MAXANT];
    double longitude;
    Double ra_p, dec_p;       // current pointing center RA,DEC at EPOCH
    Float inttime_p;
    Double freq_p;            // rest frequency of the primary line
    Int mount_p;
    Double time_p;            // current MJD time

    WINDOW win;
    Bool apply_tsys;    /* tsys weights */

    float data[2*MAXCHAN], wdata[2*MAXCHAN];	// 2*MAXCHAN since (Re,Im) pairs complex numbers
    int flags[MAXCHAN], wflags[MAXCHAN];
    float systemp[MAXANT*MAXWIDE];
};


CarmaFiller::CarmaFiller (String& infile, Int debug_level, Bool apply_tsys)
{
    num_arrays = 0;
    nfield = 0;
    npoint = 0;

    infile_p = infile;
    this->debug_level = debug_level;
    this->apply_tsys = apply_tsys;

    if (sizeof (double) != sizeof (Double))
	WARN ("sizeof(Double) != sizeof(double); pwcarmafiller will probably fail");
    if (sizeof (int) != sizeof (Int))
	WARN ("sizeof(Int) != sizeof(int); pwcarmafiller will probably fail");

    uvopen_c (&uv_handle_p, infile_p.chars (), "old");
    uvset_c (uv_handle_p, "preamble", "uvw/time/baseline", 0, 0.0, 0.0, 0.0);
    setup_tracking ();
}


bool
CarmaFiller::uv_hasvar (const char *varname)
{
    /* Also tests whether the variable has been updated if it's being tracked. */
    int vupd, vlen;
    char vtype[10];

    uvprobvr_c (uv_handle_p, varname, vtype, &vlen, &vupd);
    return vupd;
}

char *
CarmaFiller::uv_getstr (const char *varname)
{
    char *value = new char[64];
    // note: can't use sizeof(*value) since the size parameter is an int. Boo.
    uvgetvr_c (uv_handle_p, H_BYTE, varname, value, 64);
    return value;
}

int
CarmaFiller::uv_getint (const char *varname)
{
    int value;
    uvgetvr_c (uv_handle_p, H_INT, varname, (char *) &value, 1);
    // XXX: error checking!
    return value;
}

float
CarmaFiller::uv_getfloat (const char *varname)
{
    float value;
    uvgetvr_c (uv_handle_p, H_REAL, varname, (char *) &value, 1);
    return value;
}

double
CarmaFiller::uv_getdouble (const char *varname)
{
    double value;
    uvgetvr_c (uv_handle_p, H_DBLE, varname, (char *) &value, 1);
    return value;
}

void
CarmaFiller::uv_getfloats (const char *varname, float *dest, int count)
{
    uvgetvr_c (uv_handle_p, H_REAL, varname, (char *) dest, count);
}

void
CarmaFiller::uv_getdoubles (const char *varname, double *dest, int count)
{
    uvgetvr_c (uv_handle_p, H_DBLE, varname, (char *) dest, count);
}


void
CarmaFiller::checkInput ()
{
    Int i, nread, nwread;

    uvread_c (uv_handle_p, preamble, data, flags, MAXCHAN, &nread);
    uvwread_c (uv_handle_p, wdata, wflags, MAXCHAN, &nwread);
    if (nread <= 0 && nwread <= 0)
	throw AipsError ("no UV data present");

    init_window_info ();

    if (win.nspect > 0) {
	nchan_p = nread;
	nwide_p = nwread;
    } else {
	nchan_p = nread;
	nwide_p = 0;
    }

    // Get the initial array configuration
    nants_p = uv_getint ("nants");
    uv_getdoubles ("antpos", antpos, 3 * nants_p);
    longitude = uv_getdouble ("longitu");

    // Note: systemp is stored systemp[nants][nwin] in C notation
    if (win.nspect > 0)
	uv_getfloats ("systemp", systemp, nants_p * win.nspect);
    else
	uv_getfloats ("wsystemp", systemp, nants_p);

    if (win.nspect > 0)
	uv_getdoubles ("restfreq", win.restfreq, win.nspect);

    if (uv_hasvar ("project"))
	project_name = uv_getstr ("project");
    else
	project_name = "unknown";

    object_p = uv_getstr ("source");
    telescope_name = uv_getstr ("telescop");

    if (uv_hasvar ("observer"))
	observer_name = uv_getstr ("observer");
    else
	observer_name = "unknown";

    mount_p = 0;

    epoch_p = uv_getfloat ("epoch");
    epochRef_p = MDirection::J2000;
    if (nearAbs (epoch_p, 1950.0, 0.01))
	epochRef_p = MDirection::B1950;

    // TODO: these should all be handled on-the-fly.
    npol_p = uv_getint ("npol");
    pol_p = uv_getint ("pol");
    inttime_p = uv_getfloat ("inttime");
    freq_p = uv_getdouble ("freq") * 1e9; // GHz -> Hz

    ra_p = uv_getdouble ("ra");
    dec_p = uv_getdouble ("dec");

    if (hexists_c (uv_handle_p, "gains"))
	WARN ("a gains table is present, but this tool cannot apply them");
    if (hexists_c (uv_handle_p, "bandpass"))
	WARN ("a bandpass table is present, but this tool cannot apply them");
    if (hexists_c (uv_handle_p, "leakage"))
	WARN ("a leakage table is present, but this tool cannot apply them");

    uvrewind_c (uv_handle_p);

    // XXX: hardcoding assumption of full-Stokes XY pol
    npol_p = 4;
    corrType_p.resize (npol_p);
    corrType_p(0) = Stokes::XX;
    corrType_p(1) = Stokes::XY;
    corrType_p(2) = Stokes::YX;
    corrType_p(3) = Stokes::YY;
    polmapping.resize (13);
    polmapping = -1;
    polmapping(-5 + 8) = 0;
    polmapping(-6 + 8) = 3;
    polmapping(-7 + 8) = 1;
    polmapping(-8 + 8) = 2;

    corrProduct_p.resize (2, npol_p);
    corrProduct_p = 0;

    for (i = 0; i < npol_p; i++) {
	Fallible<Int> receptor = Stokes::receptor1 (Stokes::type (corrType_p(i)));
	if (receptor.isValid ())
	    corrProduct_p(0,i) = receptor;

	receptor = Stokes::receptor2 (Stokes::type (corrType_p(i)));
	if (receptor.isValid ())
	    corrProduct_p(1,i) = receptor;
    }
}


void
CarmaFiller::setupMeasurementSet (const String& ms_path)
{
    // Begin cargo-cult programming.

    TableDesc td = MS::requiredTableDesc ();
    MS::addColumnToDesc (td, MS::DATA, 2);
    td.removeColumn (MS::columnName (MS::FLAG));
    MS::addColumnToDesc (td, MS::FLAG, 2);

    td.defineHypercolumn ("TiledData", 3, stringToVector (MS::columnName (MS::DATA)));
    td.defineHypercolumn ("TiledFlag", 3, stringToVector (MS::columnName (MS::FLAG)));
    td.defineHypercolumn ("TiledUVW", 2, stringToVector (MS::columnName (MS::UVW)));

    SetupNewTable newtab (ms_path, td, Table::New);
    IncrementalStMan incrStMan ("ISMData");
    newtab.bindAll (incrStMan, True);

    Int tileSize = nchan_p / 10 + 1;

    TiledShapeStMan tiledStMan1 ("TiledData",IPosition (3, npol_p, tileSize,
							16384 / npol_p / tileSize));
    TiledShapeStMan tiledStMan1f ("TiledFlag", IPosition (3, npol_p, tileSize,
							  16384 / npol_p / tileSize));
    TiledColumnStMan tiledStMan3 ("TiledUVW", IPosition (2, 3, 1024));

    newtab.bindColumn (MS::columnName (MS::DATA), tiledStMan1);
    newtab.bindColumn (MS::columnName (MS::FLAG), tiledStMan1f);
    newtab.bindColumn (MS::columnName (MS::UVW), tiledStMan3);

    TableLock lock (TableLock::PermanentLocking);
    MeasurementSet ms (newtab, lock);
    Table::TableOption option = Table::New;

    ms.createDefaultSubtables (option);

    ms.spectralWindow ().addColumn (ArrayColumnDesc<Int>(
					MSSpectralWindow::columnName(MSSpectralWindow::ASSOC_SPW_ID),
					MSSpectralWindow::columnStandardComment(MSSpectralWindow::ASSOC_SPW_ID)));

    ms.spectralWindow ().addColumn (ArrayColumnDesc<String>(
					MSSpectralWindow::columnName(MSSpectralWindow::ASSOC_NATURE),
					MSSpectralWindow::columnStandardComment(MSSpectralWindow::ASSOC_NATURE)));

    ms.spectralWindow ().addColumn (ScalarColumnDesc<Int>(
					MSSpectralWindow::columnName(MSSpectralWindow::DOPPLER_ID),
					MSSpectralWindow::columnStandardComment(MSSpectralWindow::DOPPLER_ID)));

    // the SOURCE table; 1 extra optional column needed
    TableDesc sourceDesc = MSSource::requiredTableDesc ();
    MSSource::addColumnToDesc (sourceDesc, MSSourceEnums::REST_FREQUENCY, 1);
    SetupNewTable sourceSetup (ms.sourceTableName (), sourceDesc, option);
    ms.rwKeywordSet ().defineTable (MS::keywordName(MS::SOURCE), Table (sourceSetup));

    // the DOPPLER table; no optional columns needed
    TableDesc dopplerDesc = MSDoppler::requiredTableDesc ();
    SetupNewTable dopplerSetup (ms.dopplerTableName (),dopplerDesc, option);
    ms.rwKeywordSet ().defineTable (MS::keywordName(MS::DOPPLER), Table (dopplerSetup));

    // the SYSCAL table; 1 optional column needed
    TableDesc syscalDesc = MSSysCal::requiredTableDesc ();
    MSSysCal::addColumnToDesc (syscalDesc, MSSysCalEnums::TSYS, 1);
    SetupNewTable syscalSetup (ms.sysCalTableName (), syscalDesc, option);
    ms.rwKeywordSet ().defineTable (MS::keywordName (MS::SYSCAL), Table (syscalSetup));

    ms.initRefs(); // update the references to the subtable keywords

    {
	// Set the TableInfo. I'm assuming that the braces are to trigger a destructor.
	TableInfo& info (ms.tableInfo ());
	info.setType (TableInfo::type (TableInfo::MEASUREMENTSET));
	info.setSubType (String ("MIRIAD"));
	info.readmeAddLine ("made with pwcarmafiller");
    }

    ms_p = ms;
    msc_p = new MSColumns (ms_p);
}


void
CarmaFiller::fillObsTables ()
{
    ms_p.observation ().addRow ();
    MSObservationColumns msObsCol (ms_p.observation ());

    msObsCol.telescopeName ().put (0, telescope_name);
    msObsCol.observer ().put (0, observer_name);
    msObsCol.project ().put (0, project_name);

    MSHistoryColumns msHisCol (ms_p.history ());

    Int row = 0;
    char hline[8192]; // sigh, magic buffer sizes

    hisopen_c (uv_handle_p, "read");

    while (1) {
	int heof;

	hisread_c (uv_handle_p, hline, sizeof (hline), &heof);
	if (heof)
	    break;

	ms_p.history ().addRow ();
	msHisCol.observationId ().put (row, 0);
	msHisCol.priority ().put (row, "NORMAL");
	msHisCol.origin ().put (row, "CarmaFiller::fillObsTables");
	msHisCol.application ().put (row, "pwcarmafiller");
	msHisCol.cliCommand ().put (row, Vector<String> (0));
	msHisCol.message ().put (row, hline);
	row++;
    }

    hisclose_c (uv_handle_p);
}


void
CarmaFiller::fillMSMainTable (Int snumbase)
{
    MSColumns& msc (*msc_p);
    Int nCorr = npol_p;
    Int nChan = nchan_p;
    Int nCat  = 3; // number of flagging categories
    Int iscan = snumbase;
    Int ifield_old;

    Matrix<Complex> vis (nCorr, nChan);
    Vector<Float>   sigma (nCorr);
    Vector<String>  cat (nCat);
    cat(0) = "FLAG_CMD";
    cat(1) = "ORIGINAL";
    cat(2) = "USER";
    msc.flagCategory ().rwKeywordSet ().define ("CATEGORY", cat);
    Cube<Bool> flagCat (nCorr, nChan, nCat, False);
    Matrix<Bool> flag = flagCat.xyPlane (0); // references flagCat's storage
    Vector<Float> w1 (nCorr), w2 (nCorr);

    uvrewind_c (uv_handle_p); // may not be necessary anymore? Can't hurt ...

    nAnt_p.resize (1);
    nAnt_p[0] = 0;

    receptorAngle_p.resize (1);
    Int recnum, row = -1;
    int polsleft = 0;
    Double interval;
    Bool lastRowFlag = False;
    int nread, nwread;
    Int ant1, ant2;
    Double time;
    Vector<Double> uvw (3);

    for (recnum = 0; ; recnum++) {
	uvread_c (uv_handle_p, preamble, data, flags, MAXCHAN, &nread);
	if (nread <= 0)
	    break;

	if (win.nspect > 0)
	    uvwread_c (uv_handle_p, wdata, wflags, MAXCHAN, &nwread);
	else
	    nwread = 0;

	if (nread != nchan_p)
	    throw AipsError ("cannot handle nchan changing from " + String (nchan_p) +
			     " to " + String (nread));

	if (nwread != nwide_p)
	    throw AipsError ("cannot handle nwide changing from " + String (nwide_p) +
			     " to " + String (nwread));

	if (polsleft == 0) {
	    // starting a new simultaneous polarization record
	    uvrdvr_c (uv_handle_p, H_INT, "npol", (char *) &polsleft, NULL, 1);

	    int baseline = (int) preamble[4];
	    // XXX: we're not handling the MIRIAD >256-ant convention
	    ant1 = baseline / 256;
	    ant2 = baseline - ant1 * 256;

	    // get time in MJD seconds ; input was in JD
	    time = (preamble[3] - 2400000.5) * C::day;
	    time_p = time;
	    interval = inttime_p;

	    if (uvupdate_c (uv_handle_p))
		track_updates (); // something important changed.

	    // CARMA stuff for different "arrays" in the MS
	    nAnt_p[num_arrays-1] = max (nAnt_p[num_arrays-1], ant1);
	    nAnt_p[num_arrays-1] = max (nAnt_p[num_arrays-1], ant2);

	    // change antenna numbering convention from MIRIAD to CASA.
	    ant1--;
	    ant2--;

	    // convert ns -> m and to CASA/AIPS sign convention
	    uvw(0) = preamble[0];
	    uvw(1) = preamble[1];
	    uvw(2) = preamble[2];
	    uvw *= -1e-9 * C::c;

	    flag = 1; // clear all, in case current npol != nCorr
	    vis = 0;
	}

	int mirpol;
	Int casapolidx;
	uvrdvr_c (uv_handle_p, H_INT, "pol", (char *) &mirpol, NULL, 1);
	casapolidx = polmapping(mirpol + 8);

	if (casapolidx < 0)
	    throw AipsError ("unexpected MIRIAD polarization " + mirpol);

	Int count = 0;

	for (Int chan = 0; chan < nChan; chan++) {
	    // MIRIAD uses ant1->ant2; FITS/AIPS/CASA use ant2->ant1
	    // Along with negating UVW, we need to conjugate the visibility.
	    Bool visFlag = (flags[count/2] == 0) ? False : True;
	    Float visReal = +data[count]; count++;
	    Float visImag = -data[count]; count++;
	    Float wt = 1.0;

	    if (!visFlag)
		wt = -wt;

	    flag(casapolidx,chan) = (wt <= 0);
	    vis(casapolidx,chan) = Complex (visReal, visImag);
	}

	polsleft--;

	if (polsleft == 0 && !allTrue (flag)) {
	    // Done with this set of simultaneous pols, and not all flagged.
	    // CASA "IF" is our "spectral window" concept.

	    for (Int ifno = 0; ifno < win.nspect; ifno++) {
		// IFs go to separate rows in the MS, pol's do not!
		ms_p.addRow ();
		row++;

		if (row == 0) {
		    // first fill in values for all the unused columns
		    ifield_old = ifield;
		    msc.feed1 ().put (row, 0);
		    msc.feed2 ().put (row, 0);
		    msc.flagRow ().put (row, False);
		    lastRowFlag = False;
		    msc.scanNumber ().put (row, iscan);
		    msc.processorId ().put (row, -1);
		    msc.observationId ().put (row, 0);
		    msc.stateId ().put (row, -1);

		    if (!apply_tsys) {
			Vector<Float> tmp (nCorr);
			tmp = 1.0;
			msc.weight ().put (row, tmp);
			msc.sigma ().put (row, tmp);
		    }
		}

		msc.exposure ().put (row, interval);
		msc.interval ().put (row, interval);

		// Copying all of the data!

		Matrix<Complex> tvis (nCorr, win.nschan[ifno]);
		Cube<Bool> tflagCat (nCorr, win.nschan[ifno], nCat, False);
		Matrix<Bool> tflag = tflagCat.xyPlane (0); // references flagCat's storage

		Int woffset = win.ischan[ifno] - 1;
		Int wsize = win.nschan[ifno];
		for (int j = 0; j < nCorr; j++) {
		    for (Int i = 0; i < wsize; i++) {
			tvis(j,i) = vis(j,i+woffset);
			tflag(j,i) = flag(j,i+woffset);
		    }
		}

		msc.data ().put(row, tvis);
		msc.flag ().put(row, tflag);
		msc.flagCategory ().put (row, tflagCat);

		Bool rowFlag = allEQ (flag, True);
		if (rowFlag != lastRowFlag) {
		    msc.flagRow ().put (row, rowFlag);
		    lastRowFlag = rowFlag;
		}

		msc.antenna1 ().put (row, ant1);
		msc.antenna2 ().put (row, ant2);
		msc.time ().put (row, time); // note: CARMA timing convention changed in 2009
		msc.timeCentroid ().put (row, time);

		if (apply_tsys) {
		    w2 = 1.0; // "i use this as a 'version' id  to test FC refresh bugs :-)"

		    if (systemp[ant1] == 0 || systemp[ant2] == 0)
			w1 = 0.0;
		    else
			w1 = 1.0 / sqrt ((double) (systemp[ant1] * systemp[ant2]));

		    msc.weight ().put (row, w1);
		    msc.sigma ().put (row, w2);
		}

		msc.uvw ().put (row, uvw);
		msc.arrayId ().put (row, num_arrays - 1);
		msc.dataDescId ().put (row, ifno);
		msc.fieldId ().put (row, ifield);

		if (ifield_old != ifield)
		    iscan++;

		ifield_old = ifield;
		msc.scanNumber ().put (row, iscan);
	    }
	}
    }

    cout << infile_p << ": Processed " << recnum << " visibilities."
	 << endl;
    cout << "Found " << npoint << " pointings with "
	 << nfield << " unique source/fields "
	 << source_p.nelements() << " sources and "
	 << num_arrays << " arrays."
	 << endl;
}


void
CarmaFiller::fillAntennaTable ()
{
    // TODO: much of this should be tabulated, or preferably not hardcoded at all ...
    arrayXYZ_p.resize (3);

    if (telescope_name == "HATCREEK" || telescope_name == "BIMA") {
	arrayXYZ_p(0) = -2523862.04;
	arrayXYZ_p(1) = -4123592.80;
	arrayXYZ_p(2) =  4147750.37;
    } else if (telescope_name == "ATA") {
	// XXX: not 100% sure that this should be identical to HATCREEK/BIMA.
	arrayXYZ_p(0) = -2523862.04;
	arrayXYZ_p(1) = -4123592.80;
	arrayXYZ_p(2) =  4147750.37;
    } else if (telescope_name == "ATCA") {
	arrayXYZ_p(0) = -4750915.84;
	arrayXYZ_p(1) =  2792906.18;
	arrayXYZ_p(2) = -3200483.75;
    } else if (telescope_name == "OVRO" || telescope_name == "CARMA") {
	arrayXYZ_p(0) = -2397389.65197;
	arrayXYZ_p(1) = -4482068.56252;
	arrayXYZ_p(2) =  3843528.41479;
    } else {
	WARN ("no hardcoded array position available for name " << telescope_name);
	WARN ("assumed center of zero is grievously wrong");
	arrayXYZ_p = 0.0;
    }

    // Should use antdiam if available ...

    Float diameter = 25;
    if (telescope_name == "ATCA")
	diameter = 22;
    else if (telescope_name == "HATCREEK")
	diameter = 6;
    else if (telescope_name == "BIMA")
	diameter = 6;
    else if (telescope_name == "ATA")
	diameter = 6.1;
    else if (telescope_name == "CARMA")
	diameter = 8;
    else if (telescope_name == "OVRO")
	diameter = 10;
    else
	WARN ("no hardcoded antenna diameter for " << telescope_name <<
	      "; assuming " << diameter);

    if (nants_p == 15 && telescope_name == "OVRO") {
	WARN ("assuming CARMA-15 array (6 OVRO, 9 BIMA)");
	telescope_name = "CARMA";
    } else if (nants_p == 23 && telescope_name == "OVRO") {
	WARN ("assuming CARMA-23 array (6 OVRO, 9 BIMA, 8 SZA)");
	telescope_name = "CARMA";
    }

    Matrix<Double> posRot = Rot3D (2, longitude);

    MSAntennaColumns& ant (msc_p->antenna ());
    Vector<Double> antXYZ(3);

    if (num_arrays == 0)
	ant.setPositionRef (MPosition::ITRF);

    Int row = ms_p.antenna ().nrow () - 1;

    for (Int i = 0; i < nants_p; i++) {
	ms_p.antenna().addRow ();
	row++;

	ant.dishDiameter ().put (row, diameter);

	antXYZ(0) = antpos[i];
	antXYZ(1) = antpos[i + nants_p];
	antXYZ(2) = antpos[i + nants_p * 2];
	antXYZ *= 1e-9 * C::c; // ns -> m

	// should track the "mount" UV variable
	String mount;
	switch (mount_p) {
	case 0: mount = "ALT-AZ"; break;
	case 1: mount = "EQUATORIAL"; break;
	case 2: mount = "X-Y"; break;
	case 3: mount = "ORBITING"; break;
	case 4: mount = "BIZARRE"; break;
	default: mount = "UNKNOWN"; break;
	}

	ant.mount ().put (row, mount);
	ant.flagRow ().put (row, False);
	ant.name ().put (row, String::toString (i+1));
	ant.station ().put (row, "ANT" + String::toString (i+1));
	ant.type ().put (row, "GROUND-BASED");

	// Store absolute positions, with all offsets 0

	Vector<Double> offsets (3);
	offsets = 0.;
	antXYZ = product (posRot, antXYZ);
	ant.position ().put (row, antXYZ + arrayXYZ_p);
	ant.offset ().put (row, offsets);
    }

    num_arrays++;
    nAnt_p.resize (num_arrays);
    nAnt_p[num_arrays - 1] = 0;

    if (num_arrays > 1)
	return;

    // now do some things which only need to happen the first time around
    // store these items in non-standard keywords for now
    ant.name ().rwKeywordSet ().define ("ARRAY_NAME", telescope_name);
    ant.position ().rwKeywordSet ().define ("ARRAY_POSITION", arrayXYZ_p);
}


void
CarmaFiller::fillSyscalTable ()
{
    MSSysCalColumns& msSys (msc_p->sysCal ());
    Vector<Float> tsys(1);
    Int row = -1;

    // Note that we're using only one value for each receptor, since MIRIAD
    // has weak support for differing values (cf. xtsys and ytsys variables).

    for (Int i = 0; i < nants_p; i++) {
	ms_p.sysCal ().addRow ();
	row++;

	msSys.antennaId ().put (row, i);
	msSys.feedId ().put (row, 0);
	msSys.spectralWindowId ().put (row, -1);
	msSys.time ().put (row, time_p);
	msSys.interval ().put (row, -1.0);
	tsys(0) = systemp[i];
	msSys.tsys ().put (row, tsys);
    }
}


void
CarmaFiller::fillSpectralWindowTable ()
{
    MSSpWindowColumns& msSpW (msc_p->spectralWindow ());
    MSDataDescColumns& msDD (msc_p->dataDescription ());
    MSPolarizationColumns& msPol (msc_p->polarization ());
    MSDopplerColumns& msDop (msc_p->doppler ());

    MDirection::Types dirtype = epochRef_p;
    MEpoch ep (Quantity (time_p, "s"), MEpoch::UTC);
    MPosition obspos (MVPosition (arrayXYZ_p), MPosition::ITRF);
    MDirection dir (Quantity (ra_p, "rad"), Quantity (dec_p, "rad"), dirtype);
    MeasFrame frame (ep, obspos, dir);

    MFrequency::Types freqsys_p = MFrequency::LSRK;
    MFrequency::Convert tolsr (MFrequency::TOPO, MFrequency::Ref (freqsys_p, frame));

    // We currently handle only one polarization setup.
    ms_p.polarization ().addRow ();
    msPol.numCorr ().put (0, npol_p);
    msPol.corrType ().put (0, corrType_p);
    msPol.corrProduct ().put (0, corrProduct_p);
    msPol.flagRow ().put (0, False);

    for (Int i = 0; i < win.nspect; i++) {
	ms_p.doppler ().addRow ();
	msDop.dopplerId ().put (i, i);
	msDop.sourceId ().put (i, -1); // applies to all sources.
	msDop.transitionId ().put (i, -1);
	msDop.velDefMeas ().put (i, MDoppler (Quantity (0), MDoppler::RADIO));
    }

    for (Int i = 0; i < win.nspect; i++) {
	Int n = win.nschan[i];
	Vector<Double> f(n), w(n);

	ms_p.spectralWindow ().addRow ();
	ms_p.dataDescription ().addRow ();

	msDD.spectralWindowId ().put (i, i);
	msDD.polarizationId ().put (i, 0);
	msDD.flagRow ().put (i, False);

	msSpW.numChan ().put (i, win.nschan[i]);

	Double BW = 0.0;
	Double fwin = win.sfreq[i] * 1e9; // GHz -> Hz; a lot more of this on the way
	fwin = tolsr (fwin).getValue ().getValue ();

	for (Int j = 0; j < win.nschan[i]; j++) {
	    f(j) = fwin + j * win.sdf[i] * 1e9;
	    w(j) = abs (win.sdf[i] * 1e9);
	    BW += w(j);
	}

	msSpW.chanFreq ().put (i, f);
	if (i < win.nspect)
	    msSpW.refFrequency ().put (i, win.restfreq[i] * 1e9);
	else
	    msSpW.refFrequency ().put (i, freq_p);

	msSpW.resolution ().put (i, w);
	msSpW.chanWidth ().put (i, w);
	msSpW.effectiveBW ().put (i, w);
	msSpW.totalBandwidth ().put (i, BW);
	msSpW.ifConvChain ().put (i, 0);
	msSpW.measFreqRef ().put (i, freqsys_p);
	if (i < win.nspect)
	    msSpW.dopplerId ().put (i, i); // CARMA has only 1 ref freq line
	else
	    msSpW.dopplerId ().put (i, -1); // no ref

	if (win.sdf[i] > 0)
	    msSpW.netSideband ().put (i, 1);
	else if (win.sdf[i] < 0)
	    msSpW.netSideband ().put (i, -1);
	else
	    msSpW.netSideband ().put (i, 0);

	switch (win.code[i]) {
	case 'N':
	    msSpW.freqGroup ().put (i, 1);
	    msSpW.freqGroupName ().put (i, "MULTI-CHANNEL-DATA");
	    break;
	case 'S':
	    msSpW.freqGroup ().put (i, 2);
	    msSpW.freqGroupName ().put (i, "MULTI-CHANNEL-AVG");
	    break;
	default:
	    throw AipsError ("bad code for a spectral window");
	}
    }
}


void
CarmaFiller::fillFieldTable ()
{
    msc_p->setDirectionRef (epochRef_p);

    MSFieldColumns& msField (msc_p->field ());

    Vector<Double> radec(2), pm(2);
    Vector<MDirection> radecMeas(1);
    Double cosdec;

    pm = 0; // We don't store proper motion.

    if (nfield == 0) {
	// if no pointings found, say there is 1
	WARN ("no dra/ddec pointings found; creating one");
	nfield = npoint = 1;
	dra[0] = ddec[0] = 0.0;
    }

    for (Int fld = 0; fld < nfield; fld++) {
	int sid = sid_p[fld];

	ms_p.field ().addRow ();
	msField.sourceId ().put (fld, sid - 1);
	msField.name ().put (fld, source_p[field[fld]]);
	msField.code ().put (fld, purpose_p[field[fld]]);
	msField.numPoly ().put(fld, 0);

	cosdec = cos (dec[fld]);
	radec(0) = ra[fld] + dra[fld] / cosdec;
	radec(1) = dec[fld] + ddec[fld];

	radecMeas (0).set (MVDirection (radec(0), radec(1)), MDirection::Ref (epochRef_p));

	msField.delayDirMeasCol ().put (fld, radecMeas);
	msField.phaseDirMeasCol ().put (fld, radecMeas);
	msField.referenceDirMeasCol ().put (fld, radecMeas);

	// Need to convert epoch in years to MJD time. We're assuming UTC here
	// (and TAI elsewhere!)
	if (nearAbs (epoch_p, 2000.0, 0.01))
	    msField.time ().put (fld, MeasData::MJD2000 * C::day);
	else if (nearAbs (epoch_p, 1950.0, 0.01))
	    msField.time ().put (fld, MeasData::MJDB1950 * C::day);
	else
	    WARN ("cannot handle epoch " << epoch_p);
    }
}


void
CarmaFiller::fillSourceTable ()
{
    MSSourceColumns& msSource (msc_p->source ());
    Int srcidx = -1;
    Vector<Double> radec(2);

    for (uInt i = 0; i < source_p.nelements (); i++) {
	uInt j;

	for (j = 0; j < i; j++)
	    if (source_p[i] == source_p[j])
		break;

	if (j < i)
	    continue; // duplicate source

	srcidx++;
	ms_p.source ().addRow ();

	radec(0) = ras_p[i];
	radec(1) = decs_p[i];

	msSource.sourceId ().put (srcidx, srcidx);
	msSource.name ().put (srcidx, source_p[i]);
	// "FIX it due to a bug in MS2 code (6feb2001)":
	msSource.spectralWindowId ().put (srcidx, 0);
	msSource.direction ().put (srcidx, radec);

	if (win.nspect > 0) {
	    Vector<Double> restFreq(win.nspect);
	    for (Int i = 0; i < win.nspect; i++)
		restFreq(i) = win.restfreq[i] * 1e9;

	    msSource.numLines ().put (srcidx, win.nspect);
	    msSource.restFrequency ().put (srcidx, restFreq);
	}

	// valid at all times:
	msSource.time ().put (srcidx, 0.0);
	msSource.interval ().put (srcidx, 0);
    }
}


void
CarmaFiller::fillFeedTable ()
{
    MSFeedColumns msfc (ms_p.feed ());
    MSPolarizationColumns& msPolC (msc_p->polarization ());

    Int numCorr = msPolC.numCorr ()(0);
    Vector<String> rec_type(2);
    rec_type = "";

    if (corrType_p(0) >= Stokes::RR && corrType_p(numCorr-1) <= Stokes::LL) {
	rec_type(0) = "R";
	rec_type(1) = "L";
    }

    if (corrType_p(0) >= Stokes::XX && corrType_p(numCorr-1) <= Stokes::YY) {
	rec_type(0) = "X";
	rec_type(1) = "Y";
    }

    Matrix<Complex> polResponse(2,2);
    polResponse = 0.;
    polResponse(0,0) = polResponse(1,1) = 1.;

    Matrix<Double> offset(2,2);
    offset = 0.;

    Vector<Double> position(3);
    position = 0.;

    Vector<Double> ra(2);
    ra = 0.0;

    Int row = -1;

    for (Int arr = 0; arr < (Int) nAnt_p.nelements (); arr++) {
	for (Int ant = 0; ant < nAnt_p[arr]; ant++) {
	    ms_p.feed ().addRow ();
	    row++;

	    msfc.antennaId ().put (row, ant);
	    msfc.beamId ().put (row, -1);
	    msfc.feedId ().put (row, 0);
	    msfc.interval ().put (row, DBL_MAX);
	    msfc.spectralWindowId ().put (row, -1);
	    msfc.time ().put (row, 0.);
	    msfc.numReceptors ().put (row, 2);
	    msfc.beamOffset ().put (row, offset);
	    msfc.polarizationType ().put (row, rec_type);
	    msfc.polResponse ().put (row, polResponse);
	    msfc.position ().put (row, position);
	    msfc.receptorAngle ().put (row, ra);
	}
    }
}


void
CarmaFiller::fixEpochReferences ()
{
    String time_ref("TAI"); // hardcoded for now.

    if (time_ref == "IAT")
	time_ref = "TAI";

    if (time_ref == "UTC" || time_ref == "TAI") {
	String key ("MEASURE_REFERENCE");
	MSColumns msc (ms_p);

	msc.time ().rwKeywordSet ().define (key, time_ref);
	msc.feed ().time ().rwKeywordSet ().define (key, time_ref);
	msc.field ().time ().rwKeywordSet ().define (key, time_ref);
    } else if (time_ref != "")
	WARN ("unhandled time reference system " << time_ref);
}


void
CarmaFiller::setup_tracking ()
{
    uvtrack_c (uv_handle_p, "nschan", "u");
    uvtrack_c (uv_handle_p, "nspect", "u");
    uvtrack_c (uv_handle_p, "ischan", "u");
    uvtrack_c (uv_handle_p, "sdf", "u");
    uvtrack_c (uv_handle_p, "sfreq", "u");
    uvtrack_c (uv_handle_p, "restfreq", "u");
    uvtrack_c (uv_handle_p, "freq", "u");
    uvtrack_c (uv_handle_p, "nwide", "u");
    uvtrack_c (uv_handle_p, "wfreq", "u");
    uvtrack_c (uv_handle_p, "wwidth", "u");
    uvtrack_c (uv_handle_p, "antpos", "u");
    uvtrack_c (uv_handle_p, "dra", "u");
    uvtrack_c (uv_handle_p, "ddec", "u");
    uvtrack_c (uv_handle_p, "ra", "u");
    uvtrack_c (uv_handle_p, "dec", "u");
    uvtrack_c (uv_handle_p, "inttime", "u");
}


void
CarmaFiller::track_updates ()
{
    // "uv_hasvar" is a misnomer here. It returns true if the variable has been updated.

    if (uv_hasvar ("inttime"))
	inttime_p = uv_getfloat ("inttime");

    if (uv_hasvar ("antpos")) {
	nants_p = uv_getint ("nants");
	uv_getdoubles ("antpos", antpos, 3 * nants_p);
    }

    if (win.nspect > 0) {
	if (uv_hasvar ("systemp"))
	    uv_getfloats ("systemp", systemp, nants_p * win.nspect);
    } else {
	if (uv_hasvar ("wsystemp"))
	    uv_getfloats ("wsystemp", systemp, nants_p);
    }

    int source_updated = uv_hasvar ("source");

    if (source_updated) {
	object_p = uv_getstr ("source");

	// This leads to duplicate values; we strip them out later.
	source_p.resize (source_p.nelements () + 1, True);
	source_p[source_p.nelements () - 1] = object_p;

	ras_p.resize (ras_p.nelements () + 1, True);
	decs_p.resize (decs_p.nelements () + 1, True);
	ras_p[ras_p.nelements() - 1] = 0.0;
	decs_p[decs_p.nelements() - 1] = 0.0;

	purpose_p.resize (purpose_p.nelements () + 1, True);
	purpose_p[purpose_p.nelements () - 1] = "S";
    }

    if (source_updated || uv_hasvar ("dra") || uv_hasvar ("ddec")) {
	int i, j, k;

	npoint++;
	ra_p = uv_getdouble ("ra");
	dec_p = uv_getdouble ("dec");
	dra_p = ddec_p = 0.;
	object_p = uv_getstr ("source");

	for (i = 0, j = -1;  i < (int) source_p.nelements (); i++) {
	    if (source_p[i] == object_p) {
		j = i;
		break;
	    }
	}

	for (i = 0, k = -1; i < nfield; i++) {
	    if (dra[i] == dra_p && ddec[i] == ddec_p && field[i] == j) {
		k = i;
		break;
	    }
	}

	if (k >= 0) {
	    // This source/field combination is already known.
	    ifield = k;
	} else {
	    ifield = nfield;
	    nfield++;

	    if (nfield >= MAXFIELD)
		throw AipsError ("cannot handle more than " + String::toString (MAXFIELD) + " fields");

	    ra[ifield] = ra_p;
	    dec[ifield] = dec_p;
	    dra[ifield] = dra_p;
	    ddec[ifield] = ddec_p;
	    field[ifield] = j;
	    sid_p[ifield] = j + 1;

	    if (dra_p == 0.0 && ddec_p == 0.0) {
		// Store ra/dec for SOURCE table as well.
		ras_p[j] = ra_p;
		decs_p[j] = dec_p;
	    }
	}
    }
}


void
CarmaFiller::init_window_info ()
{
    /* "this is also a nasty routine. It makes assumptions on a relationship
       between narrow and window averages which normally exists for CARMA
       telescope data, but which in principle can be modified by uvcat/uvaver
       and possibly break this routine... (there has been some talk at the
       site to write subsets of the full data, which could break this
       routine)"
    */

    int nchan, nspect, nwide;

    if (uv_hasvar ("nchan"))
	uvrdvr_c (uv_handle_p, H_INT, "nchan", (char *) &nchan, NULL, 1);
    else
	nchan = 0;

    if (uv_hasvar ("nspect"))
	uvrdvr_c (uv_handle_p, H_INT, "nspect", (char *) &nspect, NULL, 1);
    else
	nspect = 0;

    win.nspect = nspect;

    if (uv_hasvar ("nwide"))
	uvrdvr_c (uv_handle_p, H_INT, "nwide", (char *) &nwide, NULL, 1);
    else
	nwide = 0;

    win.nwide = nwide;

    if (nspect > 0 && nspect <= MAXWIN) {
	if (uv_hasvar ("ischan"))
	    uvgetvr_c (uv_handle_p, H_INT, "ischan", (char *) win.ischan, nspect);
	else if (nspect == 1)
	    win.ischan[0] = 1;
	else
	    throw AipsError ("missing ischan");

	if (uv_hasvar ("nschan"))
	    uvgetvr_c (uv_handle_p, H_INT, "nschan", (char *) win.nschan, nspect);
	else if (nspect == 1)
	    win.nschan[0] = nchan_p;
	else
	    throw AipsError ("missing nschan");

	if (uv_hasvar ("restfreq"))
	    uv_getdoubles ("restfreq", win.restfreq, nspect);
	else
	    throw AipsError ("missing restfreq");

	if (uv_hasvar ("sdf"))
	    uv_getdoubles ("sdf", win.sdf, nspect);
	else if (nspect > 1)
	    throw AipsError ("missing sdf");

	if (uv_hasvar ("sfreq"))
	    uv_getdoubles ("sfreq", win.sfreq, nspect);
	else
	    throw AipsError ("missing sfreq");
    }

    if (nwide > 0 && nwide <= MAXWIDE) {
	if (uv_hasvar ("wfreq"))
	    uv_getfloats ("wfreq", win.wfreq, nwide);
	if (uv_hasvar ("wwidth"))
	    uv_getfloats ("wwidth", win.wwidth, nwide);
    }

    // cidx points into the combined win.xxx[] elements
    int cidx = 0;

    for (int i = 0; i < nspect; i++) {
	win.code[cidx] = 'N';
	cidx++;
    }

    for (int i = 0; i < nwide; i++) {
	Int side = win.sdf[i] < 0 ? -1 : 1;

	win.code[cidx] = 'S';
	win.ischan[cidx] = nchan + i + 1;
	win.nschan[cidx] = 1;
	win.sfreq[cidx] = win.wfreq[i];
	win.sdf[cidx] = side * win.wwidth[i];
	win.restfreq[cidx] = -1.0; // no meaning
	cidx++;
    }
}


int
main (int argc, char **argv)
{
    try {
	Input inp (1);
	inp.version ("");
	inp.create ("vis", "", "Name of CARMA dataset name", "string");
	inp.create ("ms", "", "Name of MeasurementSet", "string");
	inp.create ("tsys", "False", "Fill WEIGHT from Tsys in data?", "bool");
	inp.create ("snumbase", "0", "Starting SCAN_NUMBER value", "int");
	inp.create ("polmode", "0", "(deprecated; ignored)", "int");
	inp.readArguments (argc, argv);

	String vis (inp.getString ("vis"));
	if (vis == "")
	    throw AipsError ("no input path (vis=) given");
	if (! File (vis).isDirectory ())
	    throw AipsError ("input path (vis=) does not refer to a directory");

	String ms (inp.getString ("ms"));
	if (ms == "")
	    ms = vis.before ('.') + ".ms";

	Bool apply_tsys = inp.getBool ("tsys");
	Int snumbase = inp.getInt ("snumbase");

	// I don't understand what's going on here:
	int debug = -1;
	while (inp.debug (debug + 1))
	    debug++;

	CarmaFiller cf (vis, debug, apply_tsys);

	cf.checkInput ();
	cf.setupMeasurementSet (ms);
	cf.fillObsTables ();
	cf.fillAntennaTable ();
	cf.fillMSMainTable (snumbase);
	cf.fillSyscalTable ();
	cf.fillSpectralWindowTable ();
	cf.fillFieldTable ();
	cf.fillSourceTable ();
	cf.fillFeedTable ();
	cf.fixEpochReferences ();
    } catch (AipsError x) {
	cerr << "error: " << x.getMesg () << endl;
	return 1;
    }

    return 0;
}
