activate_this = '/Users/admin/Desktop/MainBeltComets/bin/activate_this.py'
execfile(activate_this, dict(__file__=activate_this))
import os
import sep
import vos
import numpy as np
from astropy.io import fits
from astropy.table import vstack
from astropy.time import Time
import argparse
from scipy.spatial import cKDTree
import math
import pandas as pd
from astropy.coordinates import SkyCoord

import sys

sys.path.append('User/admin/Desktop/OSSOS/MOP/src/ossos-pipeline/ossos')
from ossos import daophot
from ossos import storage
from ossos import wcs

sys.path.append('User/admin/Desktop/OSSOS/MOP/src/ossos-pipeline/planning')
from planning.plotting.horizons import batch

from get_images import get_image_info
import get_stamps

client = vos.Client()
_VOS_PATH = 'vos:kawebb/postage_stamps'
_DIR_PATH_BASE = os.path.dirname(os.path.abspath(__file__))
_STAMPS_DIR = '{}/postage_stamps'.format(_DIR_PATH_BASE)
_IMAGE_LISTS = '{}/image_lists'.format(_DIR_PATH_BASE)
_OUTPUT_DIR = '{}/phot_output'.format(_DIR_PATH_BASE)

_RADIUS = 0.01  # default radius of cutout
_ERR_ELL_RAD = 15  # default radius of the error ellipse
_MAG_ERR = 2  # default min range of magnitude from predicted by JPL
_CAT_r_TOL = 5 * 0.184 / 3600  # Tolerance for error in RA DEC position from Stephens catalogue [pixels to degrees]
_CAT_mag_TOL = 1  # Tolerance for error in measured magnitude from Stephens catalogue
_F_ERR = 10.0  # Tolerance for error in calculation focal length
_BUFFER_POLY = 10  # Buffer [pix] around polygon surrounding the object in the catagloue comparison

'''
Preforms photometry on .fits files given an input of family name and object name
Identifies object in image from predicted coordinates, elongation, and magnitude
Checks that there is enough background stars and that the object is not involved

Assumes files organised as:
dir_path_base/asteroid_families/familyname/familyname_stamps/*.fits
    - where postage stamps are temporarily copied to from VOSpace
dir_path_base/asteroid_families/familyname/*_images.txt
    - list of image exposures, predicted RA and DEC, dates etc.
vos:kawebb/postage_stamps/familyname/*.fits
    - where postage stamps are stored on VOSpace
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
    parser.add_argument("--filter",
                        action="store",
                        default='r',
                        dest="filter",
                        choices=['r', 'u'],
                        help="passband, default is r'")
    parser.add_argument('--type',
                        default='p',
                        choices=['o', 'p', 's'],
                        help="restrict type of image (unprocessed, reduced, calibrated)")

    args = parser.parse_args()

    find_objects_by_phot(args.family, args.object, float(args.aperture), float(args.thresh), args.filter, args.type)


def find_objects_by_phot(family_name, object_name, ap, th, filter_type='r', image_type='p'):
    """
    For a given family name and object name, preforms photometry and identifies object in the image.
    If only a family name is given, does the same for all objects in that family
    """

    # get CADC authentification
    username = raw_input("CADC username: ")
    password = getpass.getpass("CADC password: ")

    image_list_path = '{}/{}_images.txt'.format(_IMAGE_LISTS, family_name)

    # Retrieve the predicted coordinates of the object
    if os.path.exists(image_list_path):
        expnum_list = []
        image_list = []
        with open(image_list_path) as infile:
            filestr = infile.read()
            fileline = filestr.split('\n')
            for item in fileline[1:]:
                if len(item.split()) > 0:
                    image_list.append(item.split()[0])
                    expnum_list.append(item.split()[1])
    else:
        image_list, expnum_list, ra_list, dec_list = get_image_info(family_name, filter_type, image_type)

    tkbad_list = []
    with open('catalogue/tkBAD.txt') as infile:
        for line in infile:
            tkbad_list.append(line.split(' ')[0])

    if expnum in tkbad_list:
        print '-- Bad exposure'
    else:
        # If object name is not specified, iterate through all objects in the family
        if object_name is None:
            out_filename = '{}_phot.txt'.format(family_name)
            with open('{}/{}/{}'.format(_OUTPUT_DIR, family_name, out_filename), 'a') as outfile:
                outfile.write("{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n".format(
                    'object', 'expnum', 'ra', 'dec', 'flux', 'mag', 'x', 'y', 'time', 'consistent_f',
                    'consistent_mag', 'diff_ra', 'diff_dec', 'a', 'b', 'theta'))

            for index, image_object in enumerate(image_list):
                print 'Finding asteroid {} in family {} '.format(object_name, family_name)
                iterate_thru_images(family_name, image_object, expnum_list[index], username, password, ap, th)
        else:
            out_filename = '{}_{}_phot.txt'.format(family_name, object_name)
            with open('{}/{}/{}'.format(_OUTPUT_DIR, family_name, out_filename, 'a')) as outfile:
                outfile.write("{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n".format(
                    'object', 'expnum', 'ra', 'dec', 'flux', 'mag', 'x', 'y', 'time', 'consistent_f',
                    'consistent_mag', 'diff_ra', 'diff_dec', 'a', 'b', 'theta'))

            for index, image_object in enumerate(image_list):
                if object_name == image_object:
                    print 'Finding asteroid {} in family {} '.format(object_name, family_name)
                    iterate_thru_images(family_name, object_name, expnum_list[index], username, password, ap, th)


def iterate_thru_images(family_name, object_name, expnum_p, username, password, ap, th):
    """
    For a given family, object, and exposure number, get orbital information of the object in the image
    """
    try:
        print "-- Performing photometry on image {} ".format(expnum_p)
        septable, header, data = get_fits_data(object_name, expnum_p, family_name, ap, th)
        exptime, zeropt, size, pvwcs, start, end, size_x, size_y, size_ccd = get_header_info(header)

        print "-- Querying JPL Horizon's ephemeris"
        mag_list_jpl, r_sig = get_mag_rad(family_name, object_name)
        p_ra, p_dec, ra_dot, dec_dot = get_coords(str(object_name), start, end)

        table = append_table(septable, pvwcs, zeropt, p_ra, p_dec)
        transients, catalogue, catalog_stars = compare_to_catalogue(table)

        good_neighbours, r_err = iden_good_neighbours(object_name, transients, pvwcs, r_sig, p_ra, p_dec,
                                                      expnum_p, mag_list_jpl, ra_dot, dec_dot, exptime, family_name)
        if good_neighbours is None:
            return True

        if len(good_neighbours) == 1:
            involved = check_involvement(good_neighbours, catalogue, r_err, pvwcs, size_x, size_y)
            if involved:
                write_to_error_file2(object_name, expnum_p, out_filename='involved.txt', family_name=family_name)
                return True
            else:
                write_to_file(object_name, expnum_p, good_neighbours, header,
                              out_filename='{}_output.txt'.format(family_name), family_name=family_name)
                print '-- Cutting out recentered postage stamp'
                # cut_centered_stamp(familyname, objectname, expnum_p, good_neighbours, _RADIUS, username, password)
                return good_neighbours

        else:
            print '  More than one object identified <<<'
            return True

    except Exception, e:
        print 'ERROR: {}'.format(e)
        if 'NoneType' in e.message:
            return False
        else:
            return True
    except AssertionError:
        print 'Error in JPL query, no ephem start'
        return True


def get_fits_data(object_name, expnum_p, family_name, ap, th):
    """
    Finds image in VOSpace, determines number of extensions, inputs parameters aperture and threshold into the photometry method,
    returns photometry measurements and header values
    """

    try:
        if family_name == 'none':
            vos_dir = '{}/none'.format(_VOS_PATH)
        else:
            vos_dir = '{}/all'.format(_VOS_PATH)

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
        print 'ERROR: {}'.format(e)
        write_to_error_file2(object_name, expnum_p, out_filename="phot_param_err.txt", family_name=family_name)
        raise


def get_header_info(header):
    """
    Get values from image header
    """

    pvwcs = wcs.WCS(header)
    size_x = header['NAXIS1']
    size_y = header['NAXIS2']
    zeropt = header['PHOTZP']
    exptime = header['EXPTIME']
    size_ccd = header['CCDSIZE']
    start = '{} {}'.format(header['DATE-OBS'], header['UTIME'])
    end = '{} {}'.format(header['DATEEND'], header['UTCEND'])

    if size_x > size_y:
        size = size_x
    else:
        size = size_y

    return exptime, zeropt, size, pvwcs, start, end, size_x, size_y, size_ccd


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


def get_mag_rad(family_name, object_name):
    """
    Queries the JPL horizon's ephemeris for the variation in magnitude over the time eriod of all images of the object
    and for the radius of the error ellipse

    :param family_name:
    :param object_name:
    :return:
    """

    if type(object_name) is not str:
        object_name = str(object_name)

    # from familyname_images.txt get date range of images for objectname
    date_range = []
    with open('{}/{}_images.txt'.format(_IMAGE_LISTS, family_name)) as infile:
        for line in infile.readlines()[1:]:
            if len(line.split()) > 0:
                if object_name == line.split()[0]:
                    date_range.append(float(line.split()[5]))

    date_range_t = Time(date_range, format='mjd', scale='utc')
    assert len(date_range_t.iso) > 0
    time_start = ((date_range_t.iso[0]).split())[0] + ' 00:00:00.0'
    time_end = ((date_range_t.iso[-1]).split())[0] + ' 00:00:00.0'

    if time_start == time_end:
        time_end = add_day(date_range_t)

    print " Date range in query: {} -- {}".format(time_start, time_end)

    # output = batch("Haumea", "2010-12-28 10:00", "2010-12-29 10:00", 1, su='d')
    orbital_elements, ephemerides = batch(object_name, time_start, time_end, step=1, su='d', params=[9, 36])
    # ephemerides = query_jpl(object_name, time_start, time_end, params=[9, 36], step=1)

    mag_list = np.array(ephemerides.icol(2))
    ra_sig = np.mean(np.mean(ephemerides.icol(3)))
    dec_sig = np.mean(np.mean(ephemerides.icol(4)))

    # print '>> RA and DEC 3sigma error: {:.2f} {:.2f}'.format(ra_sig / 0.184, dec_sig / 0.184)

    if ra_sig > dec_sig:
        r_sig = ra_sig / 0.184
    else:
        r_sig = dec_sig / 0.184

    if r_sig < _ERR_ELL_RAD:  # 0.003 deg * 3600 "/deg / 0.187 "/pix
        r_sig = _ERR_ELL_RAD

    return mag_list, r_sig


def get_coords(object_name, time_start, time_end):
    """
    Queries the JPL Horizon's ephemeris for rate of change of RA and DEC for a specific day
    """

    assert type(object_name) is str
    orbital_elements, ephemerides = batch(object_name, time_start, time_end, step=1, su='m', params=[1, 3])
    # ephemerides = query_jpl(object_name, time_start, time_end, params=[1, 3], step='1', su='m')

    mid = int(len(ephemerides) / 2)
    ra = ephemerides['R.A._(ICRF/J2000.0)'][mid]
    dec = ephemerides[' DEC_(ICRF/J2000.0)'][mid]
    ra_dot_cos_dec = ephemerides[' dRA*cosD'][mid]
    dec_dot = ephemerides['d(DEC)/dt'][mid]

    ra_h = ra.split()[0]
    ra_m = ra.split()[1]
    ra_s = ra.split()[2]
    dec_d = dec.split()[0]
    dec_m = dec.split()[1]
    dec_s = dec.split()[2]

    # c = SkyCoord('00h42m30s', '+41d12m00s', frame='icrs')
    c = SkyCoord('{}h{}m{}s'.format(ra_h, ra_m, ra_s), '{}d{}m{}s'.format(dec_d, dec_m, dec_s), frame='icrs')
    ra_deg = c.ra.degree
    dec_deg = c.dec.degree

    return ra_deg, dec_deg, ra_dot_cos_dec, dec_dot


def append_table(objs, pvwcs, zeropt, p_ra, p_dec):
    """
    Calculates the RA and DEC and MAG and adds columns to the tale from the photometry
    """

    table = pd.DataFrame({'x': objs['x'],
                          'y': objs['y'],
                          'x2': objs['x2'],
                          'y2': objs['y2'],
                          'a': objs['a'],
                          'b': objs['b'],
                          # 'cxx': objs['cxx'],
                          # 'cyy': objs['cyy'],
                          'xmin': objs['xmin'],
                          'xmax': objs['xmax'],
                          'ymin': objs['ymin'],
                          'ymax': objs['ymax'],
                          'theta': objs['theta'],
                          'flux': objs['flux']
                          })

    # a2 = 0.25 * math.sqrt((xmax - xmin) ** 2 + (ymax - ymin) ** 2)
    # b2 = math.sqrt((x - math.sqrt(x2)) ** 2 + (y + math.sqrt(y2)) ** 2)
    # f = math.sqrt(a**2 - b**2)

    ra_list = []
    dec_list = []
    mag_sep_list = []
    f_list = []
    diff_ra = []
    diff_dec = []
    a2_list = []
    b2_list = []
    f_list2 = []
    x_mid = []
    y_mid = []

    for row in range(len(table)):
        ra, dec = pvwcs.xy2sky(table['x'][row], table['y'][row])
        ra_list.append(ra)
        dec_list.append(dec)
        mag_sep_list.append(-2.5 * math.log10(table['flux'][row]) + zeropt)
        a2 = 0.25 * math.sqrt(abs((table['xmax'][row] - table['xmin'][row]) ** 2 +
                                  (table['ymax'][row] - table['ymin'][row]) ** 2))
        b2 = math.sqrt(table['x2'][row] + table['y2'][row])
        a2_list.append(a2)
        b2_list.append(b2)
        f_list.append((table['a'][row] ** 2 - table['b'][row] ** 2) ** 0.5)
        f_list2.append(math.sqrt(abs(a2 ** 2 - b2 ** 2)))
        diff_ra.append((p_ra - ra) * 3600)
        diff_dec.append((p_dec - dec) * 3600)
        x_mid.append(0.5 * (table['xmax'][row] - table['xmin'][row]) + table['xmin'][row])
        y_mid.append(0.5 * (table['ymax'][row] - table['ymin'][row]) + table['ymin'][row])

    table['ra'] = ra_list
    table['dec'] = dec_list
    table['mag'] = mag_sep_list
    table['f'] = f_list
    table['diff_ra'] = diff_ra
    table['diff_dec'] = diff_dec
    table['a2'] = a2_list
    table['b2'] = b2_list
    table['f2'] = f_list2
    table['y_mid'] = y_mid
    table['x_mid'] = x_mid
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
    catalog_stars = sep_table.loc[cat_list, :]

    return transients, catalogue, catalog_stars


def find_neighbours(transients, r_sig, p_x, p_y):
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


def iden_good_neighbours(object_name, transients, pvwcs, r_sig, p_ra, p_dec,
                         expnum, mag_list_jpl, ra_dot, dec_dot, exptime, family_name):
    """
    Selects nearest neighbour object from predicted coordinates as object of interest
    In order:
        Compares measured apparent magnitude to predicted, passes if in range of values
        Calculates eccentricity, passes if greater than minimum value that is inputted
    """

    p_x, p_y = pvwcs.sky2xy(p_ra, p_dec)
    print "  Predicted RA, DEC : {:.4f}  {:.4f}".format(p_ra, p_dec)
    print "  Predicted x, y : {:.2f}  {:.2f}".format(p_x, p_y)

    # Set theoretical magnitude as mean of range with uncertainty of _MAG_ERR or magnitude variance
    mean = np.mean(mag_list_jpl)
    maxmag = np.amax(mag_list_jpl)
    minmag = np.amin(mag_list_jpl)
    if _MAG_ERR > maxmag - minmag:
        magrange = _MAG_ERR
    else:
        magrange = maxmag - minmag

    # calculate theoretical focal length
    # print '>> Error allowance is set to {} percent'.format(_F_ERR)
    # print '>> RA_dot: {:.2f}, DEC_dot: {:.2f}, exptime: {:.1f}'.format(ra_dot, dec_dot, exptime)
    f_pix = 0.5 * ((ra_dot / 2) ** 2 + (dec_dot / 2) ** 2) ** 0.5 * (exptime / (3600 * 0.184))
    f_pix_err = (_F_ERR / 100) * 0.5 * (abs(ra_dot) + abs(dec_dot)) * (exptime / (3600 * 0.184))
    assert f_pix_err != 0

    if f_pix > 10:
        print '>> Asteroid is predicted to have high elongation'

    print '  Theoretical focal length: {:.2f} +/- {:.2f}'.format(f_pix, f_pix_err)
    print '  Theoretical magnitude: {:.2f} +/- {:.2f}'.format(mean, magrange)

    print '-- Identifying object from nearest neighbours within {} pixels'.format(r_sig)
    i_list, found = find_neighbours(transients, r_sig, p_x, p_y)
    if not found:
        with open('{}/{}/no_object_found.txt'.format(_OUTPUT_DIR, family_name), 'a') as infile:
            infile.write('{}\t{}\t{}\t{}\t{}\n'.format(object_name, expnum, p_ra, p_dec, f_pix))
        raise Exception

    i_table = transients.iloc[i_list, :]  # row, column

    both_cond = i_table.query('({} < f2 < {}) & ({} < mag < {})'.format(f_pix - f_pix_err, f_pix + f_pix_err,
                                                                        mean - magrange, mean + magrange))

    only_f_cond = i_table.query('{} < f < {}'.format(f_pix - f_pix_err, f_pix + f_pix_err))
    only_mag_cond = i_table.query('{} < mag < {}'.format(mean - magrange, mean + magrange))

    if len(both_cond) > 0:
        print '>> Both conditions met by:'
        both_cond2 = add_to_object_table(both_cond, 'yes', 'yes')
        print both_cond2
        write_all_to_file(object_name, expnum, both_cond2, p_x, p_y, f_pix,
                          out_filename='{}_all_output.txt'.format(family_name), family_name=family_name)
        if len(both_cond) > 1:
            print '>> More than one object identified, writing to file <<'
            write_to_error_file(object_name, expnum, both_cond, out_filename='multiple_iden.txt',
                                family_name=family_name)
        return both_cond2, f_pix_err
    elif len(only_f_cond) > 0:
        print '>> Elongation is consistent, magnitude is not:'
        only_f_cond2 = add_to_object_table(only_f_cond, 'yes', 'no')
        print only_f_cond2
        if len(only_f_cond) > 1:
            print '>> More than one object identified, writing to file <<'
            write_to_error_file(object_name, expnum, only_f_cond2, out_filename='multiple_iden.txt',
                                family_name=family_name)
        return only_f_cond2, f_pix_err
    elif len(only_mag_cond) > 0:
        print '>> Magnitude is consistent, elongation is not:'
        only_mag_cond2 = add_to_object_table(only_mag_cond, 'no', 'yes')
        print only_mag_cond2
        if len(only_mag_cond) > 1:
            print '>> More than one object identified, writing to file <<'
            write_to_error_file(object_name, expnum, only_mag_cond2, out_filename='multiple_iden.txt',
                                family_name=family_name)
        return only_mag_cond2, f_pix_err
    else:
        print "WARNING: No condition could be satisfied <<<<<<<<<<<<<<<<<<<<<<<<<<<<"
        print '  Nearest neighbours:'
        print i_table
        write_to_error_file(object_name, expnum, i_table, out_filename='no_cond_satif.txt', family_name=family_name)
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


def check_involvement(object_data, catalogue, r_err, pvwcs, size_x, size_y):
    """
    Determine whether two objects are involved (ie overlapping psf's)
    """
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


def write_to_file(object_name, expnum_p, object_data, header, out_filename, family_name):
    """
    Prints to outfile
    """

    date = header['DATE-OBS']
    start_mjd = Time('{}T{}'.format(header['DATE-OBS'], header['UTIME']), format='isot').mjd
    end_mjd = Time('{}T{}'.format(header['DATE-OBS'], header['UTCEND']), format='isot').mjd
    time = (end_mjd - start_mjd) / 2
    date_time = '{}-{}'.format(date, time)

    if object_data is not None:
        with open('{}/{}/{}'.format(_OUTPUT_DIR, family_name, out_filename), 'a') as outfile:
            for i in range(0, len(object_data)):
                outfile.write('{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n'.format(
                    object_name, expnum_p, object_data['ra'][i], object_data['dec'][i], object_data['flux'][i],
                    object_data['mag'][i], object_data['x'][i], object_data['y'][i], date_time,
                    object_data['consistent_f'][i], object_data['consistent_mag'][i], object_data['diff_ra'][i],
                    object_data['diff_dec'][i], object_data['a'][i], object_data['b'][i], object_data['theta'][i]))

    else:
        print "WARNING: Could not identify object {} in image".format(object_name, expnum_p)


def write_all_to_file(object_name, expnum_p, object_data, p_x, p_y, p_f, out_filename, family_name):
    """
    Prints to outfile
    """

    if object_data is not None:
        print object_data
        with open('{}/{}/{}'.format(_OUTPUT_DIR, family_name, out_filename), 'a') as outfile:
            for i in range(0, len(object_data)):
                outfile.write('{} {}\n{} {} {} {} {} {} {} {}\n{} {} {} {}\n{} {} {} {} {}\n\n'.format(
                    object_name, expnum_p, p_x, p_y, object_data['x'][i], object_data['y'][i],
                    object_data['x_mid'][i], object_data['y_mid'][i],
                    object_data['x'][i] - object_data['x_mid'][i], object_data['y'][i] - object_data['y_mid'][i],
                    object_data['a'][i], object_data['b'][i], object_data['a2'][i], object_data['b2'][i],
                    p_f, object_data['f'][i], object_data['f2'][i],
                    object_data['consistent_f'][i], object_data['consistent_mag'][i]))

    else:
        print "WARNING: Could not print all output"


def write_to_error_file(object_name, expnum, object_data, out_filename, family_name):
    """
    prints a list of nearest neighbours to an outfile for human inspection
    """
    object_data.reset_index(drop=True, inplace=True)
    with open('{}/{}/{}'.format(_OUTPUT_DIR, family_name, out_filename), 'a') as outfile:
        try:
            for i in range(len(object_data)):
                outfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    object_name, expnum, object_data['ra'][i], object_data['dec'][i], object_data['flux'][i],
                    object_data['mag'][i], object_data['x'][i], object_data['y'][i],
                    object_data['diff_ra'][i], object_data['diff_dec'][i]))
            outfile.write('\n')
        except Exception, e:
            print "ERROR: cannot write to outfile {} <<<<<<<<<<<<".format(e)


def write_to_error_file2(object_name, expnum, out_filename, family_name):
    """
    Write image information to out file
    """

    with open('{}/{}/{}'.format(_OUTPUT_DIR, family_name, out_filename), 'a') as outfile:
        try:
            outfile.write('{}\t{}\n'.format(object_name, expnum))

        except Exception, e:
            print "ERROR: cannot write to outfile {} <<<<<<<<<<<<".format(e)


def cut_centered_stamp(family_name, object_name, expnum_p, object_data, r_old, username, password):
    if object_data is not None:
        # for i in range(0, len(object_data)):
        # objectname, expnum_p, r_new, RA, DEC, username, password, familyname
        print '-- Cutting recentered stamp'
        get_stamps.centered_stamp(object_name, expnum_p, r_old, object_data[0][5], object_data[0][6], username,
                                  password, family_name)


def add_day(date_range_t):
    """
    Return given date plus one day (or there abouts)
    """

    print "WARNING: only searching for one day"
    time_end_date = (((date_range_t.iso[-1]).split())[0]).split('-')
    day_add_one = int(time_end_date[2]) + 1
    if day_add_one < 10:
        day = '0{}'.format(day_add_one)
    else:
        day = day_add_one
    month = int(time_end_date[1])
    year = int(time_end_date[0])
    if day > 27:
        day = 1
        if month == 12:
            month = 1
            year += 1
        else:
            month += 1

    time_end = '{}-{}-{} 00:00:00.0'.format(year, month, day)
    return time_end


if __name__ == '__main__':
    main()
