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
from scipy.optimize import leastsq
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
_BUFFER2 = 20.  # used in the intial cutout of the polygon
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

    print '-- Getting asteroid data'
    ast_data, x_ast = get_asteroid_data(asteroid_id, data, i, fwhm)

    print '-- Getting star data'
    star_data, x_star = get_star_data(expnum, header, asteroid_id, fits_file)

    print '-- Comparing PSFs'
    compare_psf(star_data, x_star, [np.argmax(star_data), float(len(x_star)) / 2, 8., 0.], ast_data, x_ast,
                [np.argmax(ast_data), float(len(x_ast)) / 2., 5., 0.])

    compare_psf_jj(star_data, ast_data)


def compare_psf_jj(data_str, data_ast):
    """
    Calculate the difference between the two psfs
    """

    data_str_sig = np.array(data_str) ** 0.5
    data_ast_sig = np.array(data_ast) ** 0.5

    i = np.argmax(data_str) - np.argmax(
        data_ast)  # THIS IS WRONG, COMPARING TWO DIFFERENT INDEX SYSTEMS, not the ame size bins
    r = np.divide(data_ast[:i], data_str[:-i])

    r_sig = r * (np.divide(data_ast_sig[i:], data_ast[i:]) + np.divide(data_str_sig[:-i], data_str[:-i]))

    r_mean = np.mean(r)
    print (r - r_mean) / r_sig


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

                    '''
                    # Measure a spatially variable background of some image data (np array)
                    try:
                        bkg = sep.Background(data)  # , mask=mask, bw=64, bh=64, fw=3, fh=3) # optional parameters
                    except ValueError:
                        data = data.byteswap(True).newbyteorder()
                        bkg = sep.Background(data)  # , mask=mask, bw=64, bh=64, fw=3, fh=3) # optional parameters

                    # Directly subtract the background from the data in place
                    bkg.subfrom(data)
                    '''

                    return header, data, fits_file

                    # os.unlink(file_path)


def get_star_data(expnum, header, object_data, fits_file):
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
        header = hdulist[0].header
        print header['NAXIS2']

    os.unlink(local_file_path)

    # data_masked = remove_saturated_rows(data, 117)
    data_masked = data

    th = math.degrees(object_data['theta'].values)
    data_rot = rotate(data_masked, th)
    # np.savetxt('star_rot.txt', data_rot, fmt='%.0e')

    # sum all the values in each ROW
    totaled = np.ma.sum(data_rot, axis=1)

    return totaled, range(len(totaled))


def get_asteroid_data(object_data, data, i, fwhm):
    """
    Calculate psf of asteroid, taking into acount trailing effect
    """

    data_cutout, x = cutout_data(object_data, data, fwhm)

    hdu = fits.PrimaryHDU()
    hdu.data = data_cutout
    hdu.writeto('cutout_{}.fits'.format(i), clobber=True)

    # sum all the values in each ROW
    totaled = np.ma.sum(data_cutout, axis=1)

    return totaled, range(len(totaled))


def cutout_data(object_data, data, fwhm):
    """
    Cut a square around the polygon of the object, remove saturated rows, rotate so that the ellipse is parallel to the
    horizontal axis, then mask shape of ellipse around the object
    """

    # make sure there's only one object, otherwise pd.DataFrame.values won't work
    assert len(object_data) == 1

    # subtract background
    data_sub = get_bkg(data)

    # define parameters of the square to cut out from the polygon around the elongated object
    a = object_data['a'].values[0]
    b = object_data['b'].values[0]
    th = math.degrees(object_data['theta'].values)
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
    np.copyto(data_obj, data_sub[y_min:y_max, x_min:x_max])

    # rotate the data about the angle of elongation, semi major axis is along x as reading out rows (not columsn)
    data_rot = pd.DataFrame(data=rotate(data_obj, th))

    rows, cols = data_rot.shape
    cx2 = int(cols / 2)
    cy2 = int(rows / 2)

    ell_buffer = _BUFFER1 * fwhm
    print 'buffer, fwhm: ', ell_buffer, fwhm

    # instead of ellipse, try cutting out rectangle
    x_min2 = cx2 - (a + ell_buffer)
    x_max2 = cx2 + (a + ell_buffer)
    y_min2 = cy2 - (b + ell_buffer)
    y_max2 = cy2 + (a + ell_buffer)
    data_cutout = (data_rot.values)[y_min2:y_max2, x_min2:x_max2]

    return data_cutout, y_max2 - y_min2


def remove_saturated_rows(data, sat_level):
    """
    Check for saturation
    """
    print '-- Checking for saturated pixels'

    data_pd = pd.DataFrame(data=data)

    satur_row = []
    mask = []
    for row in data_pd.index:
        data_row = data_pd.iloc[row, :]
        sat = data_row.where(data_row > sat_level)
        sat_np = sat.values
        sat_masked = np.ma.masked_array(sat_np, np.isnan(sat_np))

        if np.ma.sum(sat_masked) > 5 * sat_level:
            satur_row.append(row)
            mask_row = np.zeros(len(data_row))
        else:
            mask_row = np.ones(len(data_row))

        mask.append(mask_row)

    mask_arr = np.asarray(mask)
    print 'saturated rows: ', satur_row
    data_masked = np.multiply(data, mask_arr)

    # data_pd.drop(data_pd.index[satur_row])
    # print data_pd

    return data_masked


def get_bkg(data):
    data2 = np.ones(data.shape)
    np.copyto(data2, data)
    try:
        bkg = sep.Background(data2)  # , mask=mask, bw=64, bh=64, fw=3, fh=3) # optional parameters

    except ValueError:
        data3 = data2.byteswap(True).newbyteorder()
        bkg = sep.Background(data3)       # , mask=mask, bw=64, bh=64, fw=3, fh=3) # optional parameters

    bkg.subfrom(data2)

    return data2


def gauss(x, *p):
    amp, mu, sigma, b = p
    return amp * np.exp(-(x - mu) ** 2 / (2. * sigma ** 2)) + b


def gauss2(x, *p):
    mu, sigma1, sigma2 = p
    return np.exp(-(x - mu) ** 2 / (2. * sigma1 ** 2)) + np.exp(-(x - mu) ** 2 / (2. * sigma2 ** 2))


def compare_psf(data_str, x_str, p_str, data_ast, x_ast, p_ast):
    """
    Compare psf of asteroid against mean of stars, check if anomaly in wings
    >> changed data to not normalized nor baseline subtracted
    """

    # fit a gaussian to the asteroid psf, use parameter 'b' to zero the baseline, refit \
    # use parameter 'amp' to normalize, refit
    fitp_ast, fitco_ast = curve_fit(gauss, x_ast, data_ast, p_ast)
    data_ast_norm = np.divide(data_ast, fitp_ast[0])
    fitp_ast3, fitco_ast3 = curve_fit(gauss, x_ast, data_ast_norm, p_ast)
    perr_ast3 = np.sqrt(np.diag(fitco_ast3))
    gauss_ast = fitp_ast3[0] * np.exp(-(x_ast - fitp_ast3[1]) ** 2 / (2. * fitp_ast3[2] ** 2)) + fitp_ast3[3]

    # fit a gaussian to the star psf, use parameter 'b' to zero the baseline, refit \
    # use parameter 'amp' to normalize, refit
    fitp_str, fitco_str = curve_fit(gauss, x_str, data_str, p_str)
    data_str_norm = np.divide(data_str, fitp_str[0])
    fitp_str3, fitco_str3 = curve_fit(gauss, x_str, data_str_norm, p_str)
    perr_str3 = np.sqrt(np.diag(fitco_str3))
    gauss_str = fitp_str3[0] * np.exp(-(x_str - fitp_str3[1]) ** 2 / (2. * fitp_str3[2] ** 2)) + fitp_str3[3]

    # fit the asteroid pst to a gaussian with the same mu (centerpoint) as the star
    gauss_ast2 = fitp_ast3[0] * np.exp(-(x_ast - fitp_str3[1]) ** 2 / (2. * fitp_ast3[2] ** 2)) + fitp_ast3[3]

    print fitp_str, fitp_str3  # , perr_str3
    print fitp_ast, fitp_ast3  # , perr_ast3
    '''
    with sns.axes_style('ticks'):
        plt.plot(x_str, gauss_str, label='Gauss star', ls=':')
        plt.plot(x_str, data_str_norm, label='PSF star', ls='-')
        plt.legend()
        plt.show()

    with sns.axes_style('ticks'):
        plt.plot(x_ast, gauss_ast, label='Gauss ast', ls=':')
        plt.plot(x_ast, data_ast_norm, label='PSF ast', ls='-')
        plt.legend()
        plt.show()

    with sns.axes_style('ticks'):
        print 'Plotting gausian fits of asteroid and star psfs'
        plt.plot(x_ast, gauss_ast2, label='Gauss ast', ls=':')
        plt.plot(x_str, gauss_str, label='Gauss star', ls='-.')
        plt.legend()
        plt.show()
    '''
    # calculate difference in peak values
    print np.argmax(data_str_norm), np.argmax(data_ast_norm)
    i = (x_str[np.argmax(data_str_norm)] - x_ast[np.argmax(data_ast_norm)])
    x_ast_shifted = np.add(x_ast, i)
    '''
    with sns.axes_style('ticks'):
        plt.plot(x_ast_shifted, data_ast_norm, label='PSF ast shifted', ls='--')
        plt.plot(x_str, data_str_norm, label='PSF star', ls='-.')
        plt.legend()
        plt.show()
    '''
    i2 = fitp_str3[1] - fitp_ast3[1]
    x_ast_shifted2 = np.add(x_ast, i2)
    '''
    with sns.axes_style('ticks'):
        plt.plot(x_ast_shifted2, data_ast_norm, label='PSF ast shifted', ls='--')
        plt.plot(x_str, data_str_norm, label='PSF star', ls='-.')
        plt.legend()
        plt.show()
    '''
    data_ast_longer = []
    for i in data_str_norm:
        try:
            data_ast_longer.append(data_ast[i])
        except:
            data_ast_longer.append(0)

    guess = float(fitp_str3[1]) - float(fitp_ast3[1])
    print guess, type(guess)
    fitp_shift, fitco_shirt = curve_fit(shift, data_ast_longer, x_str, p0=guess)
    print fitp_shift
    x_ast_shifted2 = np.add(x_ast, fitp_shift[0])

    with sns.axes_style('ticks'):
        plt.plot(x_ast_shifted2, data_ast_norm, label='PSF ast shifted', ls='--')
        plt.plot(x_str, data_str_norm, label='PSF star', ls='-.')
        plt.legend()
        plt.show()

def shift(x, *p):
    shifty = p
    return x + shifty

if __name__ == '__main__':
    main()