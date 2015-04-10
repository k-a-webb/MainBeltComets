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

sys.path.append('User/admin/Desktop/OSSOS/MOP/src/ossos-pipeline/ossos')
from ossos import storage
from ossos import wcs

_VOS_DIR = 'vos:kawebb/postage_stamps'
_DIR_PATH_BASE = os.path.dirname(os.path.abspath(__file__))
_STAMPS_DIR = '{}/postage_stamps'.format(_DIR_PATH_BASE)
_OSSOS_PATH_BASE = 'vos:OSSOS/dbimages'
_PHOT_DIR = '{}/phot_output'.format(_DIR_PATH_BASE)
_PSF_DIR = '{}/psf_output'.format(_DIR_PATH_BASE)

_CCD = 'EXTNAME'
_ZMAG = 'PHOTZP'
_TARGET = 'OBJECT'

_APERTURE = 10.0
_THRESHOLD = 3.0
_BUFFER1 = 2.5  # used in the shape of the ellipse
_BUFFER2 = 40.  # used in the intial cutout of the polygon
_BUFFER3 = 5  # masking center of star psf
_SATURATION_LEVEL = 3  # sigma? above background
_SATURATION_THRESHOLD = 3  # maximum allowed pixels above saturation level

_OUTPUT_NO_MKPSF = '{}/no_image_psf.txt'.format(_PSF_DIR)

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

    # detect_mbc(args.family, args.object, args.expnum)

    table = pd.read_table('{}/434/434_output_test.txt'.format(_PHOT_DIR), sep=' ', dtype={'object': object})

    for i in range(5, 6):  # len(input)):

        detect_mbc('434', table['object'][i], table['expnum'][i], i)
        print '\n'


def detect_mbc(family_name, object_name, expnum, i):
    """
    Compare psf of asteroid with mean of stars to detect possible comae
    """

    # read in asteroid identification values from the photometry output
    phot_output_file = '{}/{}/{}_output_test.txt'.format(_PHOT_DIR, family_name, family_name)
    phot_output_table = pd.read_table(phot_output_file, sep=' ', dtype={'object': object})
    asteroid_id = phot_output_table.query('object == "{}" & expnum == "{}"'.format(object_name, expnum))

    print 'buffers: ', _BUFFER1, _BUFFER2
    print asteroid_id

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
    mag = asteroid_id['mag'].values
    if mag < 18.5:
        print '>> Object is too bright for accurate photometry'
        return

    # get fwhm from OSSOS VOSpace file
    fwhm = storage.get_fwhm(expnum.strip('p'), header[_CCD].split('d')[1])

    ast_data = get_asteroid_data(asteroid_id, data, i, fwhm)
    star_data = get_star_data(expnum, header, asteroid_id, fits_file, data)

    print '-- Comparing PSFs'
    compare_psf(star_data, ast_data, fwhm)


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
    x, y = pvwcs.sky2xy(object_data['ra'].values, object_data['dec'].values)

    iraf.set(uparm="./")
    iraf.digiphot(_doprint=0)
    iraf.apphot(_doprint=0)
    iraf.daophot(_doprint=0)
    iraf.seepsf(local_file_path, local_psf, xpsf=x, ypsf=y)

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
    a = object_data['a'].values[0]
    b = object_data['b'].values[0]
    th = object_data['theta'].values
    cx = object_data['x'].values
    cy = object_data['y'].values
    a_x = a * math.cos(th)
    a_y = a * math.sin(th)
    b_x = b * math.sin(th)
    b_y = b * math.cos(th)
    p1 = (cx + b_x - a_x, cy - a_y + b_y)
    p2 = (cx + a_x + b_x, cy + a_y + b_y)
    p3 = (cx + a_x - b_x, cy + a_y - b_y)
    p4 = (cx - a_x - b_x, cy - a_y - b_y)
    polygon = Polygon([p1, p2, p3, p4]).buffer(_BUFFER2)
    x_min, y_min, x_max, y_max = polygon.bounds

    data_obj = np.ones((len(range(int(y_min), int(y_max))), len(range(int(x_min), int(x_max)))))
    np.copyto(data_obj, data[y_min:y_max, x_min:x_max])

    # rotate the data about the angle of elongation, semi major axis is along x as reading out rows (not columns)
    data_rot = pd.DataFrame(data=rotate(data_obj, math.degrees(th)))

    # subtract background
    # data_rot = get_bkg(data_rot)

    rows, cols = data_rot.shape
    cx2 = int(cols / 2)
    cy2 = int(rows / 2)

    ell_buffer = _BUFFER1 * fwhm
    print 'buffer, fwhm: ', ell_buffer, fwhm

    # cut out a circular aperture around the object
    mask = np.ones(data_rot.shape, dtype=bool)
    x_range = range(cx2 - int(a + ell_buffer), cx2 + int(a + ell_buffer))
    for x in x_range:
        y1 = int(cy2 - ((a + ell_buffer) ** 2 - (x - cx2) ** 2) ** 0.5)
        y2 = int(cy2 + ((a + ell_buffer) ** 2 - (x - cx2) ** 2) ** 0.5)
        for y in range(y1, y2):
            mask[y][x] = False

    data_cutout = np.ma.array(data_rot, mask=mask)

    '''
    # instead of ellipse, try cutting out rectangle
    x_min2 = cx2 - (a + ell_buffer)
    x_max2 = cx2 + (a + ell_buffer)
    y_min2 = cy2 - (b + ell_buffer)
    y_max2 = cy2 + (a + ell_buffer)
    data_cutout = (data_rot.values)[y_min2:y_max2, x_min2:x_max2]
    '''
    '''
    hdu = fits.PrimaryHDU()
    hdu.data = np.ma.fix_invalid(data_cutout, fill_value=0)
    hdu.writeto('cutout_{}.fits'.format(i), clobber=True)
    '''
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
    r_sig2 = r * np.sqrt(np.square(np.divide(data_ast_sig, data_ast)) +
                         np.square(np.divide(data_str_sig, data_str_at_astpts)))
    print r.compressed()
    print np.ma.mean(r)
    print np.ma.mean(np.sort(r)[2:-3])
    print np.ma.median(r)
    print '>> (r - r_mean) / r_sig'
    print (abs(r) - np.ma.mean(r)) / r_sig
    print np.ma.mean((abs(r) - np.ma.mean(r)) / r_sig)
    print (abs(r) - np.ma.mean(np.sort(r)[2:-3])) / r_sig
    print np.ma.mean((abs(r) - np.ma.mean(np.sort(r)[2:-3])) / r_sig)

    # normalize the star PSF to the height of the asteroid PSF

    ast_baseline = np.mean(np.sort(data_ast)[:8])
    str_baseline = np.mean(np.sort(data_str)[:8])
    data_str_norm = np.multiply(data_str_at_astpts,
                                np.amax(data_ast - ast_baseline) / np.amax(data_str_at_astpts - str_baseline))
    data_str_like_ast = data_str_norm + ast_baseline - np.mean(np.sort(data_str_norm)[:8])

    slope, intercept, r_value, p_value, std_err = linregress(data_ast, data_str_at_astpts)
    data_str_normed = slope * data_str_at_astpts + intercept

    with sns.axes_style('ticks'):
        plt.plot(x_ast, data_ast, label='Ast PSF', ls='--')
        plt.plot(x_str_shift, data_str, label='Star PSF', ls=':')
        plt.legend()
        plt.show()

    with sns.axes_style('ticks'):
        plt.plot(x_ast, data_ast, label='Ast PSF', ls='--')
        plt.plot(x_ast, data_str_like_ast, label='Interp star PSF', ls='-.')
        plt.legend()
        plt.show()


if __name__ == '__main__':
    main()