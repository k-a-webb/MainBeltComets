__author__ = 'kw'

activate_this = '/Users/admin/Desktop/MainBeltComets/bin/activate_this.py'
execfile(activate_this, dict(__file__=activate_this))
import pandas as pd
import numpy as np
import vos
import argparse
import os
import sep
from astropy.io import fits
from scipy.ndimage.interpolation import rotate
import math
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import linregress
from scipy.stats import chisquare
from shapely.geometry import Polygon
import sys

pd.set_option('display.max_columns', 500)
client = vos.Client()
sys.path.append('User/admin/Desktop/OSSOS/MOP/src/ossos-pipeline/ossos')
from ossos import storage
from ossos import wcs

_VOS_DIR = 'vos:kawebb/postage_stamps'
_DIR_PATH_BASE = os.path.dirname(os.path.abspath(__file__))
_STAMPS_DIR = '{}/postage_stamps'.format(_DIR_PATH_BASE)
_OSSOS_PATH_BASE = 'vos:OSSOS/dbimages'
_PHOT_DIR = '{}/phot_output'.format(_DIR_PATH_BASE)

_CCD = 'EXTNAME'
_ZMAG = 'PHOTZP'
_TARGET = 'OBJECT'

_BUFFER1 = 2.5  # aperture of asteroid is this * the fwhm.
_BUFFER2 = 30.  # the size of the cutout of the asteroid before rotation.

_OUTPUT_NO_MKPSF = '{}/psf_output/no_image_psf.txt'.format(_DIR_PATH_BASE)
_INPUT_FILE = 'all_output_test.txt'
_OUTPUT_DETECTION = 'above_3_sig.txt'

'''
headers in {}_all_output:
object expnum p_x p_y x y x_mid y_mid x-xmid y-ymid a b a2 b2 pf f f2 consist_f consist_mag
'''


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
                        default=None,
                        help='The object to be identified.')
    parser.add_argument("--expnum", '-x',
                        action='store',
                        default='1735257p',
                        help='The expsure in which the object is being identified.')

    args = parser.parse_args()

    # detect_mbc(args.family, args.object, args.expnum)

    table = pd.read_table('{}/{}/{}_{}'.format(_PHOT_DIR, args.family, args.family, _INPUT_FILE), sep=' ',
                          dtype={'object': object},
                          index_col=False)

    for i in range(6, len(table)):
        detect_mbc('434', table['object'][i], table['expnum'][i], i)
        print '\n'


def detect_mbc(family_name, object_name, expnum, i):
    """
    Compare psf of asteroid with mean of stars to detect possible activity
    """

    # read in asteroid identification values from the photometry output
    phot_output_file = '{}/{}/{}_{}'.format(_PHOT_DIR, family_name, family_name, _INPUT_FILE)
    phot_output_table = pd.read_table(phot_output_file, sep=' ', dtype={'object': object}, index_col=False)
    asteroid_id = phot_output_table.query('object == "{}" & expnum == "{}"'.format(object_name, expnum))
    print asteroid_id
    if not len(asteroid_id) == 1:
        print 'Multiple object identified'
        return

    # read in postage stamp header and data, do photometry to measure background (needed for saturation check)
    header, data, fits_file = fits_data(object_name, expnum, family_name)

    # make sure that a mean star psf has been created form the OSSOS pipeline
    if not storage.get_status(expnum.strip('p'), header[_CCD].split('d')[1], 'mkpsf'):
        expnum_ccd = '{} {}'.format(expnum.strip('p'), header[_CCD].split('d')[1])
        print 'PSF does not exist'
        with open(_OUTPUT_NO_MKPSF, 'r') as outfile:
            contents = []
            for line in outfile:
                contents.append(line.strip('\n'))
        if expnum_ccd not in contents:
            with open(_OUTPUT_NO_MKPSF, 'a') as outfile:
                outfile.write('{}\n'.format(expnum_ccd))
        return

    '''
    target = header[_TARGET]
    if 'BH' not in target:
        print '-- NOT in H Block'
        return
    '''

    # reject any object too bright that will definetly be saturated
    mag = asteroid_id['mag'].values[0]
    if mag < 18.5:
        print '>> Object is too bright for accurate photometry'
        return

    # get fwhm from OSSOS VOSpace file
    fwhm = storage.get_fwhm(expnum.strip('p'), header[_CCD].split('d')[1])

    ast_data = get_asteroid_data(asteroid_id, data, i, fwhm)
    star_data = get_star_data(expnum, header, asteroid_id, fits_file, data)

    print '-- Comparing PSFs'
    detection = compare_psf(star_data, ast_data, fwhm)
    if detection:
        write_to_file(asteroid_id, family_name)


def fits_data(object_name, expnum, family_name):
    """
    Creates local copy of fits file from VOSpace
    """

    if family_name == 'none':
        vos_dir = '{}/none'.format(_VOS_DIR)
    else:
        vos_dir = '{}/all'.format(_VOS_DIR)

    assert storage.exists(vos_dir)
    for fits_file in client.listdir(vos_dir):  # images named with convention: object_expnum_RA_DEC.fits

        if fits_file.endswith('.fits'):
            objectname_file = fits_file.split('_')[0]
            expnum_file = fits_file.split('_')[1]
            if (expnum_file == expnum) and (objectname_file == object_name):
                storage.copy('{}/{}'.format(vos_dir, fits_file), '{}/{}'.format(_STAMPS_DIR, fits_file))
                file_path = '{}/{}'.format(_STAMPS_DIR, fits_file)

                with fits.open(file_path) as hdulist:
                    data = hdulist[0].data
                    header = hdulist[0].header

                os.unlink(file_path)
                return header, data, fits_file


def get_star_data(expnum, header, object_data, fits_file, exp_data):
    """
    From ossos psf fitted image, calculate line profile
    """
    # calculate mean psf
    uri = storage.get_uri(expnum.strip('p'), header[_CCD].split('d')[1])
    ossos_psf = '{}.psf.fits'.format(uri.strip('.fits'))
    local_psf = '{}{}.psf.fits'.format(expnum, header[_CCD].split('d')[1])
    local_file_path = '{}/{}'.format(_STAMPS_DIR, local_psf)
    storage.copy(ossos_psf, local_file_path)

    # run seepsf on the mean psf image
    from pyraf import iraf

    stamp_uri = '{}/all/{}'.format(_VOS_DIR, fits_file)
    header = storage.get_header(stamp_uri)
    pvwcs = wcs.WCS(header)
    # x, y = pvwcs.sky2xy(object_data['ra'].values, object_data['dec'].values)
    x = object_data['x_mid'].values
    y = object_data['y_mid'].values

    iraf.set(uparm="./")
    iraf.digiphot(_doprint=0)
    iraf.apphot(_doprint=0)
    iraf.daophot(_doprint=0)
    iraf.seepsf(local_file_path, local_psf, xpsf=x, ypsf=y, magnitude=object_data['mag'].values)

    with fits.open(local_psf) as hdulist:
        data = hdulist[0].data

    os.unlink(local_file_path)
    os.unlink(local_psf)

    th = math.degrees(object_data['theta'].values)
    data_rot = rotate(data, th)
    data_rot = np.ma.masked_where(data_rot == 0, data_rot)

    # add in background value of exposure instead of subtracting it from the asteroid PSF
    data2 = np.ones(exp_data.shape)
    np.copyto(data2, exp_data)
    try:
        bkg = sep.Background(data2).globalback
    except ValueError:
        data3 = data2.byteswap(True).newbyteorder()
        bkg = sep.Background(data3).globalback
    data_rot += bkg

    return np.ma.mean(data_rot, axis=1)


def get_asteroid_data(object_data, data, i, fwhm):
    """
    Calculate psf of asteroid, taking into acount trailing effect
    Cut a square around the polygon of the object, remove saturated rows, rotate so that the ellipse is parallel to the
    horizontal axis, then mask shape of ellipse around the object
    """

    # make sure there's only one object, otherwise pd.DataFrame.values won't work
    assert len(object_data) == 1

    # define parameters of the square to cut out from the polygon around the elongated object
    x_min = int(object_data['xmin'].values[0] - _BUFFER2)
    x_max = int(object_data['xmax'].values[0] + _BUFFER2)
    y_min = int(object_data['ymin'].values[0] - _BUFFER2)
    y_max = int(object_data['ymax'].values[0] + _BUFFER2)

    # ensure that you are selecting data that is on the CCD
    data_y_max, data_x_max = data.shape
    if x_min < 0:
        x_min = 0
    if y_min < 0:
        y_min = 0
    if x_max > data_x_max:
        x_max = data_x_max
    if y_max > data_y_max:
        y_max = data_y_max

    data_obj = np.ones((len(range(y_min, y_max)), len(range(x_min, x_max))))
    np.copyto(data_obj, data[y_min:y_max, x_min:x_max])

    # rotate the data about the angle of elongation, semi major axis is along x as reading out rows (not columns)
    assert (math.degrees(object_data['theta'].values[0]) > -90) and (math.degrees(object_data['theta'].values[0]) < 90)
    data_rot = rotate(data_obj, math.degrees(object_data['theta'].values[0]))

    # cut out a circular aperture around the object
    a = 0.5 * math.sqrt((object_data['xmax'].values[0] - object_data['xmin'].values[0]) ** 2 +
                        (object_data['ymax'].values[0] - object_data['ymin'].values[0]) ** 2)
    ell_buffer = _BUFFER1 * fwhm
    cy, cx = np.divide(data_rot.shape, 2)

    mask = np.ones(data_rot.shape, dtype=bool)
    mask2 = np.zeros(data_rot.shape)
    x_range = range(int(cx - a - ell_buffer + 1), int(cx + a + ell_buffer - 1))
    for x in x_range:
        y1 = int(cy - math.sqrt((a + ell_buffer) ** 2 - (x - cx) ** 2))
        y2 = int(cy + math.sqrt((a + ell_buffer) ** 2 - (x - cx) ** 2))
        for y in range(y1, y2):
            mask[y][x] = False
            mask2[y][x] = 1

    data_cutout = np.ma.array(data_rot, mask=mask)
    '''
    # instead of ellipse, try cutting out rectangle
    x_min2 = cx2 - (a + ell_buffer)
    x_max2 = cx2 + (a + ell_buffer)
    y_min2 = cy2 - (b + ell_buffer)
    y_max2 = cy2 + (a + ell_buffer)
    data_cutout = (data_rot.values)[y_min2:y_max2, x_min2:x_max2]
    '''
    hdu = fits.PrimaryHDU()
    hdu.data = np.multiply(data_rot, mask2)
    hdu.writeto('cutout_{}.fits'.format(i), clobber=True)

    data_mean = np.ma.mean(data_cutout, axis=1)

    # take the mean of all the values in each ROW
    return data_mean[np.nonzero(np.ma.fix_invalid(data_mean, fill_value=0))[0]]


def get_bkg(data):
    """
    Subtract the background from the data
    """

    data2 = np.ones(data.shape)
    np.copyto(data2, data)
    try:
        bkg = sep.Background(data2)  # , mask=mask, bw=64, bh=64, fw=3, fh=3) # optional parameters

    except ValueError:
        data3 = data2.byteswap(True).newbyteorder()
        bkg = sep.Background(data3)  # , mask=mask, bw=64, bh=64, fw=3, fh=3) # optional parameters

    bkg.subfrom(data2)

    return data2


def gauss(x, *p):
    amp, mu, sigma, b = p
    return amp * np.exp(-(x - mu) ** 2 / (2. * sigma ** 2)) + b


def compare_psf(data_str, data_ast, fwhm):
    """
    Compare psf of asteroid against mean of stars, check if anomaly in wings
    >> changed data to not normalized nor baseline subtracted
    """

    x_ast = range(len(data_ast))
    x_str = range(len(data_str))
    p_str = [np.amax(data_str), float(len(x_str)) / 2, fwhm, 0.]
    p_ast = [np.amax(data_ast), float(len(x_ast)) / 2., fwhm, 0.]

    # fit a gaussian to the asteroid psf and then the star psf
    fitp_ast, fitco_ast = curve_fit(gauss, x_ast, data_ast, p_ast)
    fitp_str, fitco_str = curve_fit(gauss, x_str, data_str, p_str)
    gauss_ast = fitp_ast[0] * np.exp(-(x_ast - fitp_ast[1]) ** 2 / (2. * fitp_ast[2] ** 2)) + fitp_ast[3]
    '''
    with sns.axes_style('ticks'):
        plt.plot(x_ast, data_ast, label='Ast PSF', ls='--')
        plt.plot(x_ast, gauss_ast, label='Gaussian fit', ls=':')
        plt.legend()
        plt.show()
    '''
    # shift the center of the star psf to the center of the asteroid psf
    gauss_shift = fitp_str[1] - fitp_ast[1]
    x_str_shift = np.subtract(x_str, gauss_shift)

    # interpolate the psf values of the star at the x values that there is data for the asteroid psf
    y_str = lambda x: np.interp(x, x_str_shift, data_str)
    data_str_at_astpts = y_str(x_ast)

    print ">> Ratio test ast/str"
    data_str_sig = np.absolute(np.array(data_str_at_astpts)) ** 0.5
    data_ast_sig = np.absolute(np.array(data_ast)) ** 0.5
    r = np.divide(data_ast, data_str_at_astpts)
    r_sig = r * (np.divide(data_ast_sig, data_ast) + np.divide(data_str_sig, data_str_at_astpts))
    r_sig2 = r * np.sqrt(
        np.square(np.divide(data_ast_sig, data_ast)) + np.square(np.divide(data_str_sig, data_str_at_astpts)))

    print r.compressed()
    print np.ma.mean(r), np.ma.mean(np.sort(r)[2:-3])
    print '>> (r - r_mean) / r_sig'
    ratio = (abs(r) - np.ma.mean(np.sort(r)[2:-3])) / r_sig2
    print ratio
    print np.ma.mean(ratio)

    # normalize the star PSF to the height of the asteroid PSF
    ast_baseline = np.mean(np.sort(data_ast)[:8])
    str_baseline = np.mean(np.sort(data_str)[:8])
    data_str_norm = np.multiply(data_str_at_astpts,
                                np.amax(data_ast - ast_baseline) / np.amax(data_str_at_astpts - str_baseline))
    data_str_like_ast = data_str_norm + ast_baseline - np.mean(np.sort(data_str_norm)[:8])

    # slope, intercept, r_value, p_value, std_err = linregress(data_ast, data_str_at_astpts)
    # data_str_normed = slope * data_str_at_astpts + intercept
    '''
    with sns.axes_style('ticks'):
        plt.plot(x_ast, data_ast, label='Ast PSF', ls='--')
        plt.plot(x_str_shift, data_str, label='Star PSF', ls=':')
        plt.legend()
        plt.show()
    '''
    with sns.axes_style('ticks'):
        plt.plot(x_ast, data_ast, label='Ast PSF', ls='--')
        plt.plot(x_ast, data_str_like_ast, label='Interp star PSF', ls='-.')
        plt.legend()
        plt.show()

    off_mean_pts = []
    for i in ratio:
        if abs(i) > 3:
            off_mean_pts.append(i)
    if len(off_mean_pts) > 1:
        print '>> Detect possible comae <<'
        return True
    else:
        print '>> No detection'
        return False


def write_to_file(asteroid_id, family_name):
    with open('{}/{}/{}'.format(_PHOT_DIR, family_name, _OUTPUT_DETECTION), 'a') as outfile:
        outfile.write('{} {}'.format(asteroid_id['object'].values[0], asteroid_id['expnum'].values[0]))


if __name__ == '__main__':
    main()