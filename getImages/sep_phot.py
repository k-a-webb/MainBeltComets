activate_this = '/Users/admin/Desktop/MainBeltComets/bin/activate_this.py'
execfile(activate_this, dict(__file__=activate_this))
import os
import sep
import vos
import getpass
import numpy as np
from astropy.io import fits
from astropy.time import Time
import argparse
from scipy.spatial import cKDTree
import math
import pandas as pd
from astropy.coordinates import SkyCoord
import sys
from get_images import get_image_info
import get_stamps

sys.path.append('User/admin/Desktop/OSSOS/MOP/src/ossos-pipeline/ossos')
from ossos import daophot
from ossos import storage
from ossos import wcs

sys.path.append('User/admin/Desktop/OSSOS/MOP/src/ossos-pipeline/planning')
from planning.plotting.horizons import batch

pd.set_option('display.max_columns', 500)
client = vos.Client()

_VOS_DIR = 'vos:kawebb/postage_stamps/'
_DIR_PATH_BASE = os.path.dirname(os.path.abspath(__file__))
_STAMPS_DIR = '{}/postage_stamps'.format(_DIR_PATH_BASE)
_IMAGE_LISTS = '{}/image_lists'.format(_DIR_PATH_BASE)
_OUTPUT_DIR = '{}/phot_output'.format(_DIR_PATH_BASE)

_INPUT_FILE = 'images.txt'

_OUTPUT_MULT_ID = 'multiple_iden.txt'
_OUTPUT_NO_COND_MET = 'no_cond_met.txt'
_OUTPUT_PHOT_ERR = "phot_param_err.txt"
_OUTPUT_PHOT_MEAS_ERR = 'phot_meas_err.txt'
_OUTPUT_INVOLVED = 'involved.txt'
_OUTPUT_NOT_FOUND = 'no_object_found.txt'
_OUTPUT_ALL = 'output.txt'
_OUTPUT_ERROR = 'unknown_error.txt'
_OUTPUT_VOS_ERR = 'vos_err.txt'

_RADIUS = 0.01  # default radius of cutout
_ERR_ELL_RAD = 30  # minimum radius of nearest neighbour search
_MAG_ERR = 2  # default min range of magnitude from predicted by JPL
_CAT_r_TOL = 5 * 0.184 / 3600  # Tolerance for error in RA DEC position from Stephens catalogue [pixels to degrees]
_CAT_mag_TOL = 1  # Tolerance for error in measured magnitude from Stephens catalogue
_F_ERR = 10.0  # Tolerance for error in calculation focal length
_BUFFER_POLY = 10  # Buffer [pix] around polygon surrounding the object in the catagloue comparison
_A_MAX = 50  # Maximum pixels expected to be measured for semi-major axis
_B_MAX = 40  # Maximum pixels expected to be measured for semi-minor axis

'''
Preforms photometry on .fits files given an input of family name and object name
Identifies object in image from predicted coordinates, elongation, and magnitude
Checks that there is enough background stars and that the object is not involved

Assumes files organised as:
dir_path_base/postage_stamps/familyname_stamps/*.fits
    - where postage stamps are temporarily copied to from VOSpace
dir_path_base/image_lists/*_images.txt
    - list of image exposures, predicted RA and DEC, dates etc.
vos:kawebb/postage_stamps/familyname/*.fits
    - where postage stamps are stored on VOSpace
    - familyname is either 'all' or 'none' depending on family designation
'''


def main():
    parser = argparse.ArgumentParser(
        description='For a given set of fits images: preforms photometry, identifies a specific object, and returns \
                        the orbital elements of the object')
    parser.add_argument("--family", '-f',
                        action="store",
                        default='all',
                        help="Asteroid family name. Usually the asteroid number of the largest member.")
    parser.add_argument("--aperture", '-ap',
                        action='store',
                        default=10.0,
                        help='aperture (degree) of circle for photometry.')
    parser.add_argument("--thresh", '-t',
                        action='store',
                        default=3.,
                        help='threshold value for photometry (sigma above background).')
    parser.add_argument("--object", '-o',
                        action='store',
                        default=None,
                        help='The object to be identified.')
    parser.add_argument("--expnum", '-x',
                        action='store',
                        default=None,
                        help='The exposure that the object is in.')
    parser.add_argument('--type',
                        default='p',
                        choices=['o', 'p', 's'],
                        help="restrict type of image (unprocessed, reduced, calibrated)")

    args = parser.parse_args()

    if args.object is None:
        find_objects_by_phot(args.family, float(args.aperture), float(args.thresh), args.type)
    else:
        assert args.expnum.endswith('p') or args.expnum.endswith('o'), "Must include 'p' or 'o' in expnum"
        image_list_path = '{}/{}_{}'.format(_IMAGE_LISTS, args.family, _INPUT_FILE)
        table = pd.read_table(image_list_path, sep='\t', dtype={'Object': object, 'Image': object}, index_col=False)
        lines = table[(table.Object == args.object) & (table.Image == args.expnum)]
        print lines
        username = raw_input("CADC username: ")
        password = getpass.getpass("CADC password: ")

        postage_stamp_filename = "{}_{}_{:8f}_{:8f}.fits".format(args.object, args.expnum, float(lines.RA.values[0]),
                                                                 float(lines.DEC.values[0]))
        if not storage.exists('{}/all/{}'.format(_VOS_DIR, postage_stamp_filename)):
            print '-- Cutout not found, creating new cutout'
            get_stamps.cutout(username, password, args.family, args.object, args.expnum, lines.RA.values[0],
                              lines.DEC.values[0], _RADIUS)
        assert storage.exists('{}/all/{}'.format(_VOS_DIR, postage_stamp_filename)), 'Error: Cutout not found'
        iterate_thru_images(args.family, args.object, args.expnum, username, password, float(args.aperture),
                            float(args.thresh))


def find_objects_by_phot(family_name, aperture, thresh, imagetype):
    """
    For a given family name and object name, preforms photometry and identifies object in the image.
    If only a family name is given, does the same for all objects in that family
    """

    # get CADC authentification
    username = raw_input("CADC username: ")
    password = getpass.getpass("CADC password: ")

    # establish input/output
    if family_name == 'none':
        vos_dir = '{}/none'.format(_VOS_DIR)
    else:
        vos_dir = '{}/all'.format(_VOS_DIR)
    assert storage.exists(vos_dir)
    image_list_path = '{}/{}_{}'.format(_IMAGE_LISTS, family_name, _INPUT_FILE)  # USING TEST FILE
    assert os.path.exists(image_list_path)

    # Remove any fits files hanging around from failed run
    for fits_file in os.listdir(_STAMPS_DIR):
        if fits_file.endswith('.fits'):
            storage.remove('{}/{}'.format(_STAMPS_DIR, fits_file))

    tkbad_list = []
    with open('catalogue/tkBAD.txt') as infile:
        for line in infile:
            tkbad_list.append(line.split(' ')[0])

    # for each image of each object, make cutout and go to sep_phot
    table = pd.read_table(image_list_path, sep='\t', dtype={'Object': object}, index_col=False)

    for row in range(len(table)):
        print "\n {} --- Searching for asteroid {} in image {} ".format(row, table['Object'][row], table['Image'][row])

        expnum = (table['Image'][row]).strip('{}'.format(imagetype))
        if expnum in tkbad_list:
            print '-- Bad exposure'

        else:
            postage_stamp_filename = "{}_{}_{:8f}_{:8f}.fits".format(table['Object'][row], table['Image'][row],
                                                                     table['RA'][row], table['DEC'][row])

            if not storage.exists('{}/{}'.format(vos_dir, postage_stamp_filename)):
                print '-- Cutout not found, creating new cutout'
                get_stamps.cutout(username, password, family_name, table['Object'][row], table['Image'][row],
                                  table['RA'][row], table['DEC'][row], _RADIUS)
            if not storage.exists('{}/{}'.format(vos_dir, postage_stamp_filename)):
                with open('{}/{}/{}'.format(_OUTPUT_DIR, family_name, _OUTPUT_VOS_ERR), 'a') as outfile:
                    outfile.write('{} {}\n'.format(table['Object'][row], table['Image'][row]))
            else:
                success = False
                attempts = 0
                while (success is False) and (attempts < 3):
                    success = iterate_thru_images(family_name, str(table['Object'][row]), table['Image'][row], username,
                                                  password, aperture, thresh)
                    attempts += 1
                    if attempts == 3:
                        print ' >>>> Last attempt \n'


def iterate_thru_images(family_name, object_name, expnum_p, username, password, ap, th):
    """
    For a given family, object, and exposure number, get orbital information of the object in the image
    """
    try:
        septable, header, data = get_fits_data(object_name, expnum_p, family_name, ap, th)

        print "-- Querying JPL Horizon's ephemeris"
        mag_list_jpl, r_sig = get_mag_rad(object_name, header)
        p_ra, p_dec, ra_dot, dec_dot = get_coords(str(object_name), header)

        table = append_table(septable, header, p_ra, p_dec)
        transients, catalogue = compare_to_catalogue(table)

        objs_met_cond, r_err, p_f, p_mag = iden_good_neighbours(family_name, object_name, transients, header, r_sig,
                                                                p_ra, p_dec, expnum_p, mag_list_jpl, ra_dot, dec_dot)

        if (len(objs_met_cond) != 1) or (len(objs_met_cond) is None):
            return True

        involved = check_involvement(objs_met_cond, catalogue, r_err, header)
        if involved:
            write_to_error_file(object_name, expnum_p, _OUTPUT_INVOLVED, family_name)
            return True

        write_to_file(family_name, object_name, expnum_p, objs_met_cond, header, p_ra, p_dec, p_f,
                      p_mag)
        print '-- Cutting out recentered postage stamp'
        # cut_centered_stamp(familyname, objectname, expnum_p, objs_met_cond, _RADIUS, username, password)
        return objs_met_cond

    except Exception, e:
        print 'Error iterating through images: {}'.format(e)
        if 'NoneType' in e.message:
            return False
        else:
            with open('{}/{}/{}'.format(_OUTPUT_DIR, family_name, _OUTPUT_ERROR), 'a') as outfile:
                outfile.write('{} {}\n'.format(object_name, expnum_p))
            return True
    except AssertionError, e:
        print e
        print 'Error in JPL query, no ephem start'
        with open('{}/{}/{}'.format(_OUTPUT_DIR, family_name, _OUTPUT_ERROR), 'a') as outfile:
            outfile.write('{} {}\n'.format(object_name, expnum_p))
        return True


def get_fits_data(object_name, expnum_p, family_name, ap, th):
    """
    Finds image in VOSpace, determines number of extensions, inputs parameters aperture and threshold into the photometry method,
    returns photometry measurements and header values
    """

    try:
        if family_name == 'none':
            vos_dir = '{}/none'.format(_VOS_DIR)
        else:
            vos_dir = '{}/all'.format(_VOS_DIR)

        assert storage.exists(vos_dir)
        for fits_file in client.listdir(vos_dir):  # images named with convention: object_expnum_RA_DEC.fits

            if fits_file.endswith('.fits'):
                objectname_file = fits_file.split('_')[0]
                expnum_file = fits_file.split('_')[1]

                if (expnum_file == expnum_p) and (objectname_file == object_name):
                    file_path = '{}/{}'.format(_STAMPS_DIR, fits_file)
                    storage.copy('{}/{}'.format(vos_dir, fits_file), file_path)

                    data = fits.getdata(file_path)
                    header = fits.getheader(file_path)
                    objs = sep_phot(data, ap, th)

                    os.unlink(file_path)

                    return objs, header, data

    except TypeError:
        print "WARNING: Image does not exist for {} {}".format(object_name, expnum_p)
        raise
    except Exception, e:
        print 'Error retrieving fits data: {}'.format(e)
        write_to_error_file(object_name, expnum_p, out_filename=_OUTPUT_PHOT_ERR, family_name=family_name)
        raise


def sep_phot(data, ap, th):
    """
    Preforms photometry by SEP, similar to source extractor
    """

    # Measure a spatially variable background of some image data (np array)
    try:
        bkg = sep.Background(data)  # , mask=mask, bw=64, bh=64, fw=3, fh=3) # optional parameters
    except ValueError:
        data = data.byteswap(True).newbyteorder()
        bkg = sep.Background(data)  # , mask=mask, bw=64, bh=64, fw=3, fh=3) # optional parameters

    # Directly subtract the background from the data in place
    bkg.subfrom(data)

    # for the background subtracted data, detect objects in data given some threshold
    thresh = th * bkg.globalrms  # ensure the threshold is high enough wrt background
    objs = sep.extract(data, thresh)
    # calculate the Kron radius for each object, then we perform elliptical aperture photometry within that radius
    kronrad, krflag = sep.kron_radius(data, objs['x'], objs['y'], objs['a'], objs['b'], objs['theta'], ap)
    flux, fluxerr, flag = sep.sum_ellipse(data, objs['x'], objs['y'], objs['a'], objs['b'], objs['theta'],
                                          2.5 * kronrad, subpix=1)
    flag |= krflag  # combine flags into 'flag'

    r_min = 1.75  # minimum diameter = 3.5
    use_circle = kronrad * np.sqrt(objs['a'] * objs['b']) < r_min
    x = objs['x']
    y = objs['y']
    cflux, cfluxerr, cflag = sep.sum_circle(data, x[use_circle], y[use_circle],
                                            r_min, subpix=1)
    flux[use_circle] = cflux
    fluxerr[use_circle] = cfluxerr
    flag[use_circle] = cflag

    return objs


def get_mag_rad(object_name, header):
    """
    Queries the JPL horizon's ephemeris for the variation in magnitude over the time eriod of all images of the object
    and for the radius of the error ellipse
    """

    # search for values for 10 days past observation
    date_start = Time(header['DATE-OBS'], format='iso', scale='utc')
    date_end = Time(date_start.mjd + 10, format='mjd', scale='utc')

    print " Date range in magnitude query: {} -- {}".format(date_start.iso, date_end.iso)

    orbital_elements, ephemerides = batch(object_name, date_start.iso, date_end.iso, step=1, su='d', params=[9, 36])

    mag_list = np.array(ephemerides.icol(2))
    ra_sig = np.mean(np.mean(ephemerides.icol(3)))
    dec_sig = np.mean(np.mean(ephemerides.icol(4)))

    if ra_sig > dec_sig:
        r_sig = ra_sig / 0.184
    else:
        r_sig = dec_sig / 0.184

    if r_sig < _ERR_ELL_RAD:  # 0.003 deg * 3600 "/deg / 0.187 "/pix
        r_sig = _ERR_ELL_RAD

    return mag_list, r_sig


def get_coords(object_name, header):
    """
    Queries the JPL Horizon's ephemeris for rate of change of RA and DEC for a specific day
    """

    time_start = '{} {}'.format(header['DATE-OBS'], header['UTIME'])
    time_end = '{} {}'.format(header['DATEEND'], header['UTCEND'])

    assert type(object_name) is str
    orbital_elements, ephemerides = batch(object_name, time_start, time_end, step=1, su='m', params=[1, 3])

    mid = int(len(ephemerides) / 2)
    ra = ephemerides['R.A._(ICRF/J2000.0)'][mid]
    dec = ephemerides[' DEC_(ICRF/J2000.0)'][mid]
    ra_dot_cos_dec = ephemerides[' dRA*cosD'][mid]
    dec_dot = ephemerides['d(DEC)/dt'][mid]

    c = SkyCoord('{}h{}m{}s'.format(ra.split()[0], ra.split()[1], ra.split()[2]), '{}d{}m{}s'.format(
        dec.split()[0], dec.split()[1], dec.split()[2]), frame='icrs')
    ra_deg = c.ra.degree
    dec_deg = c.dec.degree

    return ra_deg, dec_deg, ra_dot_cos_dec, dec_dot


def append_table(objs, header, p_ra, p_dec):
    """
    Calculates the RA and DEC and MAG and adds columns to the tale from the photometry
    """

    table = pd.DataFrame({'x': objs['x'],
                          'y': objs['y'],
                          'x2': objs['x2'],
                          'y2': objs['y2'],
                          'a': objs['a'],
                          'b': objs['b'],
                          'xmin': objs['xmin'],
                          'xmax': objs['xmax'],
                          'ymin': objs['ymin'],
                          'ymax': objs['ymax'],
                          'theta': objs['theta'],
                          'flux': objs['flux']
                          })

    ra_list = []
    dec_list = []
    pvwcs = wcs.WCS(header)
    for row in range(len(table)):
        ra, dec = pvwcs.xy2sky(table['x'][row], table['y'][row])
        ra_list.append(ra)
        dec_list.append(dec)

    table['ra'] = ra_list
    table['dec'] = dec_list
    table['diff_ra'] = np.subtract(ra_list, p_ra) * 3600
    table['diff_dec'] = np.subtract(dec_list, p_dec) * 3600
    table['mag'] = -2.5 * np.log10(table['flux'].values) + header['PHOTZP']
    table['f'] = np.sqrt(np.square(table['a'].values) - np.square(table['b'].values))
    table['x_mid'] = 0.5 * (table['xmax'].values - table['xmin'].values) + table['xmin'].values
    table['y_mid'] = 0.5 * (table['ymax'].values - table['ymin'].values) + table['ymin'].values

    return table


def compare_to_catalogue(sep_table):
    """
    Compares objects detected by photometry to those in Stephens background object catalogue
    Builds a list of transient objects from objects not in the catalogue
    """

    blocks = ['Hall', 'Lall', 'Oall']  # starting from the second block
    catalogue = pd.read_table('catalogue/Eall.photcat', usecols=[0, 1, 4], header=0, names=['ra', 'dec', 'mag'],
                              sep='      |     |    |   |  ', engine='python')
    for block in blocks:
        temp_table = pd.read_table('catalogue/{}.photcat'.format(block), usecols=[0, 1, 4], header=0,
                                   names=['ra', 'dec', 'mag'], sep='      |     |    |   |  ', engine='python')
        catalogue = pd.concat([catalogue, temp_table])

    catalogue.reset_index(drop=True, inplace=True)

    trans_list = []
    cat_list = []

    for row in range(len(sep_table)):
        # index = catalogue[( abs(catalogue.ra - sep_table['ra'][row]) < 0.0051111 )]
        index = catalogue.query(
            '({} < ra < {}) & ({} < dec < {}) & ({} < mag < {})'.format(sep_table['ra'][row] - _CAT_r_TOL,
                                                                        sep_table['ra'][row] + _CAT_r_TOL,
                                                                        sep_table['dec'][row] - _CAT_r_TOL,
                                                                        sep_table['dec'][row] + _CAT_r_TOL,
                                                                        sep_table['mag'][row] - _CAT_mag_TOL,
                                                                        sep_table['mag'][row] + _CAT_mag_TOL))
        if len(index) == 0:
            trans_list.append(row)

        else:
            cat_list.append(row)

    if len(trans_list) == 0:
        print "WARNING: No transients identified"

    transients = sep_table.loc[trans_list, :]

    return transients, catalogue


def find_neighbours(transients, r_sig, p_x, p_y):
    """
    Parse through sep output to find objects within a radius r_sig of the predicted location of the asteroid
    """
    found = True
    i_list = neighbour_search(transients, r_sig, p_x, p_y)
    if len(i_list) == 0:
        print '-- Expanding radius of nearest neighbour search by 1.5x'
        i_list = neighbour_search(transients, 2 * r_sig, p_x, p_y)
        if len(i_list) == 0:
            print 'WARNING: No nearest neighbours were found within {} ++++++++++++++++++++++'.format(r_sig * 1.5)
            found = False

    return i_list, found


def neighbour_search(transients, r_sig, p_x, p_y):
    """
    Computes the nearest neighbours to predicted coordinates within an RA/DEC uncertainty circle
    """
    tree = cKDTree(zip((np.array(transients['x'])).ravel(), (np.array(transients['y'])).ravel()))

    # parse through table and get RA and DEC closest to predicted coordinates (in pixels)
    coords = np.array([p_x, p_y])
    i_list = tree.query_ball_point(coords, r_sig)

    return i_list


def iden_good_neighbours(family_name, object_name, transients, header, r_sig, p_ra, p_dec,
                         expnum, mag_list_jpl, ra_dot, dec_dot):
    """
    Selects nearest neighbour object from predicted coordinates as object of interest
    In order:
        Compares measured apparent magnitude to predicted, passes if in range of values
        Calculates eccentricity, passes if greater than minimum value that is inputted
    """

    pvwcs = wcs.WCS(header)
    p_x, p_y = pvwcs.sky2xy(p_ra, p_dec)
    print "  Predicted RA, DEC : {:.4f}  {:.4f}".format(p_ra, p_dec)
    print "  Predicted x, y : {:.2f}  {:.2f}".format(p_x, p_y)

    # Set theoretical magnitude as mean of range with uncertainty of _MAG_ERR or magnitude variance
    p_mag = np.mean(mag_list_jpl)
    if _MAG_ERR > np.amax(mag_list_jpl) - np.amin(mag_list_jpl):
        mag_err = _MAG_ERR
    else:
        mag_err = np.amax(mag_list_jpl) - np.amin(mag_list_jpl)

    # calculate theoretical focal length
    p_f = 0.5 * ((ra_dot / 2) ** 2 + (dec_dot / 2) ** 2) ** 0.5 * (header['EXPTIME'] / (3600 * 0.184))
    p_f_err = (_F_ERR / 100) * 0.5 * (abs(ra_dot) + abs(dec_dot)) * (header['EXPTIME'] / (3600 * 0.184))
    assert p_f_err != 0, 'Focal length calculated to be zero'

    if p_f > 10:
        print '>> Asteroid is predicted to have high elongation'
    print '  Theoretical focal length: {:.2f} +/- {:.2f}'.format(p_f, p_f_err)
    print '  Theoretical magnitude: {:.2f} +/- {:.2f}'.format(p_mag, mag_err)

    # parse SEP output for transients in radius r_sig from predicted location of the asteroid
    i_list, found = find_neighbours(transients, r_sig, p_x, p_y)
    if not found:
        write_not_found(family_name, object_name, expnum, p_ra, p_dec, p_x, p_y, p_f, p_mag)
        return

    # build a table of the candidate transients, check if meet elongation and magnitude conditions
    i_table = transients.iloc[i_list, :]  # row, column

    # check that SEP phot values make sense
    '''
    for i in i_list:
        if not (i_table['theta'][i] > - np.pi / 2) or not (i_table['theta'][i] < np.pi / 2) or not \
                (i_table['a'][i] < _A_MAX) or not (i_table['b'][i] < _B_MAX):
            print 'PHOT measurement error: {} {} {}'.format(i_table['theta'][i], i_table['a'][i], i_table['b'][i])
            write_phot_meas_err(family_name, p_ra, p_dec, p_x, p_y, p_f, p_mag, i_table, i)
            return
    '''
    both_cond = i_table.query('({} < f < {}) & ({} < mag < {})'.format(p_f - p_f_err, p_f + p_f_err,
                                                                       p_mag - mag_err, p_mag + mag_err))
    only_f_cond = i_table.query('{} < f < {}'.format(p_f - p_f_err, p_f + p_f_err))
    only_mag_cond = i_table.query('{} < mag < {}'.format(p_mag - mag_err, p_mag + mag_err))

    if len(both_cond) > 0:
        print '>> Both conditions met by:'
        both_cond2 = add_to_object_table(both_cond, 'yes', 'yes')
        print both_cond2
        if len(both_cond2) > 1:
            write_multp_id(object_name, expnum, both_cond2, family_name, p_ra, p_dec, p_x, p_y, p_f, p_mag)
        return both_cond2, p_f_err, p_f, p_mag
    elif len(only_f_cond) > 0:
        print '>> Elongation is consistent, magnitude is not:'
        only_f_cond2 = add_to_object_table(only_f_cond, 'yes', 'no')
        print only_f_cond2
        if len(only_f_cond2) > 1:
            write_multp_id(object_name, expnum, only_f_cond2, family_name, p_ra, p_dec, p_x, p_y, p_f, p_mag)
        return only_f_cond2, p_f_err, p_f, p_mag
    elif len(only_mag_cond) > 0:
        print '>> Magnitude is consistent, elongation is not:'
        only_mag_cond2 = add_to_object_table(only_mag_cond, 'no', 'yes')
        print only_mag_cond2
        if len(only_mag_cond2) > 1:
            write_multp_id(object_name, expnum, only_mag_cond2, family_name, p_ra, p_dec, p_x, p_y, p_f, p_mag)
        return only_mag_cond2, p_f_err, p_f, p_mag
    else:
        print "WARNING: No condition could be satisfied <<<<<<<<<<<<<<<<<<<<<<<<<<<<"
        print '  Nearest neighbours:'
        print i_table
        write_no_cond_met(object_name, expnum, i_table, family_name, p_x, p_y, p_ra, p_dec, p_f, p_mag)
        return


def add_to_object_table(table, first, second):
    """
    Adds 'yes' or 'no' to consistent_f and consistent_mag columns on identified object table
    """

    first_arr = []
    second_arr = []

    for i in range(len(table)):
        first_arr.append(first)
        second_arr.append(second)

    table2 = pd.DataFrame({'consistent_f': first_arr,
                           'consistent_mag': second_arr})

    table.reset_index(drop=True, inplace=True)
    table3 = table.join(table2)

    return table3


def check_involvement(object_data, catalogue, r_err, header):
    """
    Determine whether two objects are involved (ie overlapping psf's)
    """

    pvwcs = wcs.WCS(header)
    size_x = header['NAXIS1']
    size_y = header['NAXIS2']

    object_data.reset_index(drop=True, inplace=True)
    ra = object_data['ra'][0]
    dec = object_data['dec'][0]
    min_x = object_data['xmin'][0] - r_err
    min_y = object_data['ymin'][0] - r_err
    max_x = object_data['xmax'][0] + r_err
    max_y = object_data['ymax'][0] + r_err

    if (max_x > size_x) or (min_x < 0) or (max_y > size_y) or (min_y < 0):
        print ">> Object is on the edge"
        return

    min_ra, min_dec = pvwcs.xy2sky(min_x, min_y)
    max_ra, max_dec = pvwcs.xy2sky(max_x, max_y)

    ra_range = abs(min_ra - max_ra)
    dec_range = abs(min_dec - max_dec)

    cat_list = catalogue.query('({} < ra < {}) & ({} < dec < {})'.format(ra - ra_range, ra + ra_range,
                                                                         dec - dec_range, dec + dec_range))
    if len(cat_list) > 0:
        print ">> Object is involved <<<<<<<<<<<<<<<<<<<"
        print cat_list
        return True
    else:
        print '>> Object is not involved'
        return False


def write_to_file(family_name, object_name, expnum, object_data, header, p_ra, p_dec, p_f, p_mag):
    """
    Prints to outfile
    """

    date_start = Time('{} {}'.format(header['DATE-OBS'], header['UTIME']), format='iso', scale='utc')
    date_end = Time('{} {}'.format(header['DATE-OBS'], header['UTCEND']), format='iso', scale='utc')
    date_mid = 0.5 * (date_end.mjd - date_start.mjd) + date_start.mjd

    pvwcs = wcs.WCS(header)
    p_x, p_y = pvwcs.sky2xy(p_ra, p_dec)

    if object_data is not None:
        with open('{}/{}/{}_{}'.format(_OUTPUT_DIR, family_name, family_name, _OUTPUT_ALL), 'a') as outfile:
            for i in range(len(object_data)):
                outfile.write(
                    '{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n'.format(
                        object_name, expnum, p_ra, p_dec, object_data['ra'][i], object_data['dec'][i],
                        p_x, p_y, object_data['x'][i], object_data['y'][i], object_data['x_mid'][i],
                        object_data['y_mid'][i],
                        object_data['xmin'][i], object_data['xmax'][i], object_data['ymin'][i], object_data['ymax'][i],
                        object_data['a'][i], object_data['b'][i], object_data['theta'][i],
                        p_f, object_data['f'][i], object_data['consistent_f'][i],
                        p_mag, object_data['mag'][i], object_data['flux'][i], object_data['consistent_mag'][i],
                        date_mid))

    else:
        print "WARNING: Could not identify object {} in image".format(object_name, expnum_p)


def write_no_cond_met(object_name, expnum, object_data, family_name, p_x, p_y, p_ra, p_dec, p_f, p_mag):
    """
    prints a list of nearest neighbours to an outfile for human inspection
    """
    object_data.reset_index(drop=True, inplace=True)
    with open('{}/{}/{}'.format(_OUTPUT_DIR, family_name, _OUTPUT_NO_COND_MET), 'a') as outfile:
        try:
            outfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\nnearest neighbours:\n'.format(
                object_name, expnum, p_ra, p_dec, p_x, p_y, p_f, p_mag))
            for i in range(len(object_data)):
                outfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    object_data['ra'][i], object_data['dec'][i], object_data['x'][i], object_data['y'][i],
                    object_data['f'][i], object_data['mag'][i], object_data['diff_ra'][i], object_data['diff_dec'][i]))
            outfile.write('\n')
        except Exception, e:
            print "ERROR: cannot write to outfile {} <<<<<<<<<<<<".format(e)


def write_to_error_file(object_name, expnum, out_filename, family_name):
    """
    Write image information to out file
    """

    with open('{}/{}/{}'.format(_OUTPUT_DIR, family_name, out_filename), 'a') as outfile:
        try:
            outfile.write('{}\t{}\n'.format(object_name, expnum))

        except Exception, e:
            print "ERROR: cannot write to outfile {} <<<<<<<<<<<<".format(e)


def write_multp_id(object_name, expnum, object_data, family_name, p_ra, p_dec, p_x, p_y, p_f, p_mag):
    """
    prints a list of nearest neighbours to an outfile for human inspection
    """
    print '>> More than one object identified, writing to file <<'
    object_data.reset_index(drop=True, inplace=True)
    with open('{}/{}/{}'.format(_OUTPUT_DIR, family_name, _OUTPUT_MULT_ID), 'a') as outfile:
        try:
            outfile.write(
                '{} {} {} {} {} {} {} {}\ncandidates:\n'.format(object_name, expnum, p_ra, p_dec, p_x, p_y, p_f, p_mag))
            for i in range(len(object_data)):
                outfile.write('{} {} {} {} {} {} {} {}\n'.format(
                    object_data['ra'][i], object_data['dec'][i], object_data['x'][i], object_data['y'][i],
                    object_data['f'][i], object_data['mag'][i], object_data['consistent_f'][i],
                    object_data['consistent_mag'][i]))
            outfile.write('\n')
        except Exception, e:
            print "ERROR: cannot write to outfile {} <<<<<<<<<<<<".format(e)


def write_not_found(family_name, object_name, expnum, p_ra, p_dec, p_x, p_y, p_f, p_mag):
    with open('{}/{}/{}'.format(_OUTPUT_DIR, family_name, _OUTPUT_NOT_FOUND), 'a') as infile:
        infile.write('{} {} {} {} {} {} {} {}\n'.format(object_name, expnum, p_ra, p_dec, p_x, p_y, p_f, p_mag))


def write_phot_meas_err(family_name, p_ra, p_dec, p_x, p_y, p_f, p_mag, i_table, row):
    with open('{}/{}/{}'.format(_OUTPUT_DIR, family_name, _OUTPUT_PHOT_MEAS_ERR), 'a') as outfile:
        outfile.write('{} {} {} {} {} {} {} {}\ncandidates:\n{} {} {} {} {} {}'.format(
            i_table['object_name'][row], i_table['expnum'][row], p_ra, p_dec, p_x, p_y, p_f, p_mag,
            i_table['ra'][row], i_table['dec'][row], i_table['x'][row], i_table['y'][row],
            i_table['f'][row], i_table['mag'][row]))


def cut_centered_stamp(family_name, object_name, expnum_p, object_data, r_old, username, password):
    if object_data is not None:
        # for i in range(0, len(object_data)):
        # objectname, expnum_p, r_new, RA, DEC, username, password, familyname
        print '-- Cutting recentered stamp'
        get_stamps.centered_stamp(object_name, expnum_p, r_old, object_data[0][5], object_data[0][6], username,
                                  password, family_name)


if __name__ == '__main__':
    main()
