__author__ = 'kw'

activate_this = '/Users/admin/Desktop/MainBeltComets/bin/activate_this.py'
execfile(activate_this, dict(__file__=activate_this))
import pandas as pd
import numpy as np
import vos
import argparse
import os
from astropy.io import ascii
import sep

import sep_phot
import ossos_scripts.wcs as wcs
# from ossos_scripts import psf

import sys

sys.path.append('User/admin/Desktop/OSSOS/MOP/src/ossos-pipeline/ossos')
from ossos import daophot
from ossos import storage

sys.path.append('/Users/admin/Desktop/PythonPhot/PythonPhot')
from PythonPhot import getpsf, aper

_VOS_DIR = 'vos:kawebb/postage_stamps/all'
_DIR_PATH_BASE = os.path.dirname(os.path.abspath(__file__))
_DIR_PATH = '{}/asteroid_families'.format(_DIR_PATH_BASE)
_STAMPS_DIR = '{}/all/all_stamps'.format(_DIR_PATH)
_OSSOS_PATH_BASE = 'vos:OSSOS/dbimages'

_CCD_HEADER = 'EXTNAME'
_ZMAG_HEADER = 'PHOTZP'

_APERTURE = 10.0
_THRESHOLD = 3.0
_MAG_LIM = 20.0  # minimum brightness of catalog star to compare psf with
client = vos.Client()


def main():
    parser = argparse.ArgumentParser(
        description='For a given set of fits images: preforms photometry, identifies a specific object, and returns \
                        the orbital elements of the object')
    parser.add_argument("--family", '-f',
                        action="store",
                        default='434',
                        help="Asteroid family name. Usually the asteroid number of the largest member.")
    parser.add_argument("--object", '-o',
                        action='store',
                        default='25076',
                        help='The object to be identified.')
    parser.add_argument("--expnum", '-x',
                        action='store',
                        default='1735257p',
                        help='The expsure in which the object is being identified.')

    args = parser.parse_args()

    detect_mbc(args.family, args.object, args.expnum)


def detect_mbc(family_name, object_name, expnum):
    """
    Compare psf of asteroid with mean of stars to detect possible comae
    """
    phot_output_file = '{}/{}/phot_output/{}_output_test.txt'.format(_DIR_PATH, family_name, family_name)
    phot_output = pd.read_table(phot_output_file, sep='\t', dtype={'object': object})
    asteroid_id = phot_output.query('object == "{}" & expnum == "{}"'.format(object_name, expnum))

    catalog_stars, x_list, y_list, header, data = iden_catalog_stars(object_name, expnum, asteroid_id)

    #file_path = 'asteroid_families/all/all_stamps/25076_1735257p_27.654770_14.390550.fits'
    fits_file, file_path = copy_fits_to_local(object_name, expnum)

    # change apcor?
    pdump = daophot.phot(file_path, x_list, y_list, _APERTURE, exptime=header[_CCD_HEADER], zmag=header[_ZMAG_HEADER])
    pdump_table = pd.DataFrame(np.array(pdump))
    print pdump_table


    ossos_path = '{}/{}/{}'.format(_OSSOS_PATH_BASE, expnum.strip('p'), header[_CCD_HEADER])
    file_psf = '{}{}.psf.fits'.format(expnum, header[_CCD_HEADER].split('d')[1])
    file_psf_path = '{}/{}'.format(ossos_path, file_psf)

    storage.copy(file_psf_path, _STAMPS_DIR)

    # psf.build_psf(expnum, header['EXTNAME'])
    #build_psf(data, header, pdump_table)


def build_psf(data, header, pdump):

    gain = header['GAIN']
    zeropt = header['PHOTZP']

    try:
        bkg = sep.Background(data)  # , mask=mask, bw=64, bh=64, fw=3, fh=3) # optional parameters
    except ValueError:
        data = data.byteswap(True).newbyteorder()
        bkg = sep.Background(data)  # , mask=mask, bw=64, bh=64, fw=3, fh=3) # optional parameters
    sky = bkg.back()

    xpos = pdump['XCENTER'].values
    ypos = pdump['YCENTER'].values
    mag = pdump['MAG'].values

    # run aper to get mags and sky values for specified coords
    # mag, magerr, flux, fluxerr, sky, skyerr, badflag, outstr =  \
    #    aper.aper(data, xpos, ypos, phpadu=gain, apr=_APERTURE, zeropoint=zeropt)

    # image, xc, ys, apmag, sky, ronois, phpadu, idpsf, psfrad, fitrad, psfname
    gauss, psf, psfmag = getpsf.getpsf(data, xpos, ypos, mag, sky, ronois=1, phpadu=gain, idpsf=np.arange(len(xpos)),
                                       psfrad=_APERTURE, fitrad=3., psfname='output_psf.fits')

    print gauss
    print psf
    print psfmag


def copy_fits_to_local(object_name, expnum):
    """
    Creates local copy of fits file from VOSpace
    """

    assert storage.exists(_VOS_DIR)
    for fits_file in client.listdir(_VOS_DIR):  # images named with convention: object_expnum_RA_DEC.fits
        if fits_file.endswith('.fits'):
            objectname_file = fits_file.split('_')[0]
            expnum_file = fits_file.split('_')[1]

            if (expnum_file == expnum) and (objectname_file == object_name):
                print '-- Stamp found'
                file_path = '{}/{}'.format(_STAMPS_DIR, fits_file)
                storage.copy('{}/{}'.format(_VOS_DIR, fits_file), file_path)

                return fits_file, file_path


def iden_catalog_stars(object_name, expnum, asteroid_id):
    """
    Compare photometry output with Stephen's catalogue to identify stars to use for PSF comparison
    """

    septable, header, data = sep_phot.get_fits_data(object_name, expnum, _APERTURE, _THRESHOLD)
    table = sep_phot.append_table(septable, wcs.WCS(header), header['PHOTZP'], asteroid_id['ra'][0],
                                  asteroid_id['dec'][0])
    transients, catalogue, num_cat_objs, catalog_stars = sep_phot.compare_to_catalogue(table)

    bright_stars = catalog_stars.query('mag < {}'.format(_MAG_LIM))
    bright_stars.reset_index(drop=True, inplace=True)

    return bright_stars, bright_stars['x'].values, bright_stars['y'].values, header, data

'''
def xbuild_psf(expnum, header, filename='test.txt'):
    # The itime is set to 1 as the SGwyn ZP includes the exposure time.
    DEFAULT_INTEGRATION_TIME = 1
    DEFAULT_READ_NOISE = 4
    DATA_MAX = 30000
    DATA_MIN = -1000
    MAXIMUM_LINEAR_COUNTS = 30000
    BRIGHT_STAR_THRESHOLD = 4
    NOMINAL_FWHM = 4.0
    PSF_ORDER = 3
    SMALL_APERTURE_FRACTION = 1.1
    LARGE_APERTURE_FRACTION = 5.0

    # Width of sky annulus in units of fwhm
    DEFAULT_SKY_DANNULUS = 1

    ZEROPOINT_HEADER_KEYWORD = header["PHOTZP"]
    READ_NOISE_KEYWORD = header["RDNOISE"]
    UTC_OBS_HEADER_KEYWORD = header["UTC-OBS"]
    INTEGRATION_TIME_HEADER_KEYWORD = header["EXPTIME"]
    AIRMASS_HEADER_KEYWORD = header["AIRMASS"]
    GAIN_HEADER_KEYWORD = header["GAIN"]

    try:
        fwhm = storage.get_fwhm(expnum, header['EXTNAME'])  # possibly EXTVER ?
    except IOError as err:
        print 'ERROR: {}'.format(err)
        fwhm = NOMINAL_FWHM
    apmax = int(10 * LARGE_APERTURE_FRACTION * fwhm + 0.5) / 10.0
    apmin = int(10 * SMALL_APERTURE_FRACTION * fwhm + 0.5) / 10.0
    try:
        zeropoint = storage.get_zeropoint(expnum, header['EXTNAME'])
    except:
        print 'Zeropoint error, refer to psf.py and storage.py'
        # return build_zeropoint_used()

    from pyraf import iraf
    from iraf import noao, digiphot, daophot

    iraf.unlearn('datapars')
    iraf.unlearn('findpars')
    iraf.unlearn('centerpars')
    iraf.unlearn('fitskypars')
    iraf.unlearn('photpars')
    iraf.unlearn('daopars')
    iraf.unlearn('daofind')
    iraf.unlearn('phot')

    iraf.daophot.verify = iraf.no
    iraf.daophot.verbose = iraf.no
    iraf.daophot.update = iraf.no

    iraf.psf.showplots = iraf.no
    iraf.psf.interactive = iraf.no

    iraf.datapars.exposure = ""  # default: ""       Exposure time image header keyword
    iraf.datapars.datamin = 0  # default: INDEF    Minimum good data value
    iraf.datapars.itime = DEFAULT_INTEGRATION_TIME  # default: 1.       Exposure time
    iraf.datapars.gain = GAIN_HEADER_KEYWORD  # default: ""       CCD gain image header keyword
    iraf.datapars.airmass = AIRMASS_HEADER_KEYWORD  # default: ""       Airmass image header keyword
    iraf.datapars.obstime = UTC_OBS_HEADER_KEYWORD  # default: ""       Time of observation image header keyword
    iraf.datapars.ccdread = READ_NOISE_KEYWORD  # default: ""       CCD readout noise image header keyword
    iraf.datapars.readnoise = DEFAULT_READ_NOISE  # default: 0.0      CCD readout noise in electrons
    iraf.datapars.datamax = DATA_MAX  # default: INDEF    Maximum good data value
    iraf.datapars.datamin = DATA_MIN  # default: INDEF    Minimum good data value
    iraf.datapars.fwhmpsf = fwhm  # default: 2.5      FWHM of the PSF in scale units

    iraf.photpars.apertures = apmax  # default: "3."     List of aperture radii in scale units
    iraf.photpars.apertures = apmin
    iraf.photpars.zmag = zeropoint  # default: 25.0     Zero point of magnitude scale

    iraf.fitskypars.salgori = "centroid"  # default: "mode"   Sky fitting algorithm
    iraf.fitskypars.dannulus = DEFAULT_SKY_WIDTH * fwhm  # default: 10.      Width of sky annulus in scale units
    iraf.fitskypars.annulus = apmax + 2  # default: 10.0     Inner radius of sky annulus in scale units

    iraf.centerpars.calgorithm = "centroid"  # default: "none"   Centering algorithm
    iraf.centerpars.maxshift = 3  # default: 1.       Maximum center shift in scale units

    # iraf.daopars.varorder = 1                           # default: 0        Order of empirical component of psf model
    iraf.daopars.psfrad = apmax  # default: 11.0     Radius of psf model in scale units
    iraf.daopars.fitrad = apmax  # default: 3.0      Fitting radius in scale units
    # iraf.daopars.fitsky = 'yes'                         # default: no       Recompute group sky value during fit?
    # iraf.daopars.sannulus = 8                           # default: 0.0      Inner radius of sky fitting annulus in scale units
    # iraf.daopars.wsannulus = 5.                         # default: 11.0     Width of sky fitting annulus in scale units
    # iraf.daopars.proferr = 2.5                          # default: 5.0      Profile error in percent
    # iraf.daopars.clipexp = 4                            # default: 6        Bad data clipping exponent
    # iraf.daopars.maxnstar = 1000                        # default: 10000    Maximum number of stars to fit
    # iraf.daopars.maxgroup = 30                          # default: 60       Maximum number of stars to fit per group

    # phot output file
    outphot = 'test.phot'
    photlog = 'phot.log'

    # Delete any existing phot output and log
    for f in (outphot, photlog):
        if os.path.exists(f):
            os.remove(f)

    previous_fwhm = None
    niters = 0
    self.iraf.datapars.datamin = self.datamin

    while previous_fwhm is None or abs(previous_fwhm - fwhm) / fwhm > 0.01 and niters < 5:
        niters += 1
        p_fwhm = fwhm

        # We should use a minimum aperture that is about 1.1 -1.2 the FWHM
        # This is set at the start of the psf loop so the PSF small aperture
        # matches that used in mkapfile after the PSF is built.
        # We do this twice since the first time we build
        # the PSF with the input FWHM
        # and the second time using the FWHM determined from the PSF stars

        # Adjust the large aperture size too, and the sky fitting parameters
        self.iraf.fitskypars.dannulus = DEFAULT_SKY_DANNULUS * apmin
        self.iraf.fitskypars.annulus = apmax + 2

        # set the PSF radius and fitting radius
        self.iraf.daopars.fitrad = apmin
        self.iraf.daopars.psfrad = apmax
        self.iraf.datapars.fwhm = fwhm
        # self.iraf.datamin = datamin

        # time to select upto 25 PSF stars, from the good stars we have left.

        pst2_filename = os.path.splitext(filename)[0] + ".pst.2"
        psg1_filename = os.path.splitext(filename)[0] + ".psg.1"
        os.path.unlink(pst2_filename)
        os.path.unlink(psg1_filename)
        os.path.unliin(self.psf_filename)

        # Estimate the FWHM using the gaussian fit to the PSF
        hselect(t_psf, "PAR1", yes) | scan(xsig)
        hselect(t_psf, "PAR2", yes) | scan(ysig)
        hselect(t_psf, "NPSFSTAR", yes) | scan(npsfstar)
        t_fwhm = sqrt((xsig * xsig + ysig * ysig) / 2.0)

        print "Sigma -> {}".format(t_fwhm)

        if abs((xsig - ysig) / t_fwhm) > 0.2:  # JJK Added a warning message here. March 2002
            touch(t_image // ".psf.WARNING")
            print "ERROR: The PSF is too elliptical for planting. image {}".format(t_image) # we must have 10 or more stars to get a good psf...
            print "{}.psf.WARNING".format(t_image)
        if npsfstar < 9:
            print "Only {} where used to build the image, terminating".format(npsfstar)
            failedflag=1
            goto
            finalproc

        # for MEGACAM the FWHM is about 84% of the
        # value you would get based on a guassian fit to the PSF
        t_fwhm = 0.84 * t_fwhm * sqrt(8 * log(2))
        print "Got FWHM of {}".format(t_fwhm)

        if t_fwhm > 8:
            touch (t_image // ".psf.WARNING")
            print "WARNING: Seeing value large? Measured {}".format(t_fwhm)
            print "{}.psf.WARNING".format(t_image)

        # keep going until change in t_fwhm between
        # loops is less than 0.05 (5%)

    # # Check that the output fwhm is reasonably close to the input one.  # and overwrite the input value with the computed one
    list = "{}.fwhm".format(t_image)
    dum = fscan(list, p_fwhm)
    if (p_fwhm > t_fwhm * 3.0) | | (t_fwhm > p_fwhm * 3.0):
        print "PSF change too  large"
        failedflag = 1
        goto
        finalproc


    kdelete("{}.fwhm".format(t_image))
    print t_fwhm, >> "{}.fwhm".format(t_image)  # compute the apcor using the PSF stars.
    kdelete("{}.nst.1".format(t_image))
    kdelete("{}.nrj.1".format(t_image))
    nstar(t_image, "{}.psg.1".format(t_image),
          psfimage = "{}.psf.fits".format(t_image),
          nstarfile = "{}.nst.1".format(t_image),
          rejfile = "{}.nrj.1".format(t_image))

    kdelete(t_image // ".phot".format(t_image))
    pselect("{}.nst.1".format(t_image), "{}.phot".format(t_image), "CHI < 2.5")  # build a file for use as input into phot
    kdelete("{}.coo.2".format(t_image))
    txdump("{}.nst.1".format(t_image), "XCENTER,YCENTER,MAG,SHARPNESS,ID", "CHI < 2.5",
           headers +, > "{}.coo.2".format(t_image))  # Run PHOT with a range of apertures.
    kdelete("{}.mag.2".format(t_image))
    photpars.apertures = "{}:{}:1".format(apmin, apmax)
    naps = apmax - apmin + 1

    if naps < 3:
        print "inner ({}) and outer ({}) aperatures must differ by more than 2 pixels\n".format(apmin, apmax)
        failedflag = 1
        goto
        finalproc

    phot(t_image, "{}.coo.2".format(t_image), "{}.mag.2".format(t_image)) # fit that range of apertures to a curve of growth.
    kdelete("{}.mkap".format(t_image))

    failedflag = 1
    if err:
        mkapfile("{}.mag.2".format(t_image), naps, "{}.mkap".format(t_image), interactive -, verify -, nparams=t_order)
        goto
        finalproc

    failedflag = 0  # Read in the aperture correction file and convert to CFEPS Pipeline format
    list = "{}.mkap".format(t_image)
    line = "{} {} 30 30 ".format(apmin, apmax) # the last line has the actual correction

    ac = 0
    while fscan(list, line) != EOF:
        ac += 1

    if ac != 3:
        print "Invalid mkapfile output? (check {}.mkap)".format(t_image)
        failedflag = 1
        goto
        finalproc

    # Stick the apcor line in the way the pipeline wants it
    apcor = -30
    aperr = 30
    wc = fscan(line, s1, apcor, aperr)
    kdelete(t_image // ".apcor")
    apcor *= -1
    print "{.1f} {.1f} {6.2f} {6.2f} {}.apcor\n".format(apmin, apmax, apcor, aperr, t_image)


    # spit out some warnings about apcor, if its a bit off.
    if apmin < 2:
        touch ("{}.apcor.WARNING".format(t_image))
        print "WARNING: apmin small: {}".format(apmin)
        print "{}.apcor.WARNING".format(t_image)

    if apmax > 35:
        touch ("{}.apcor.WARNING".format(t_image))
        print "WARNING: apmax huge: {}".format(apmax)
        print "{}.apcor.WARNING".format(t_image)

    if apcor > 1.0:
        touch ("{}.apcor.WARNING".format(t_image))
        print "ERROR: apcor large"
        print "{}.apcor.WARNING".format(t_image)

    if aperr > 0.2:
        touch ("{}.apcor.WARNING".format(t_image))
        print "ERROR: aperr large"
        print "{}.apcor.WARNING".format(t_image)


    finalproc:

    if t_keep_files:
        print "Keeping intermediate files"
    else:
        delete(t_image // ".bright.mag", verify-)
        delete(t_image // ".mag*", verify-)
        delete(t_image // ".pst*", verify-)
        delete(t_image // ".psg*", verify-)
        delete(t_image // ".coo*", verify-)
        delete(t_image // ".nrj*", verify-)
        delete(t_image // ".nst*", verify-)

    if ( failedflag == 0 ):
        if ( imaccess(t_psf) ):
            print("\n" // t_psf // " OK\n")
            touch (t_image // ".jmpmakepsf.OK")
        else:
            print("\n NO PSF BUILT \n")
            failedflag = 1
            pwd | scan(t_base)


    if (failedflag != 0):
        print("\n jmpmakepsf: FAILED in " // t_base // " for " // t_image // "\n", >> t_image // ".jmpmakepsf.FAILED")
        print "failed to build psf"
'''


def meas_psf():
    """
    Calculate psf of stars (and asteroid?)
    """
    # calculate mean psf


def meas_psf_asteroid():
    """
    Calculate psf of asteroid, taking into acount trailing effect
    """


def compare_psf():
    """
    Compare psf of asteroid against mean of stars, check if anomaly in wings
    """


if __name__ == '__main__':
    main()