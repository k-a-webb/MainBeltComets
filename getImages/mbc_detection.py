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
import copy

from scipy.ndimage.interpolation import rotate

import math
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
from shapely.geometry import Polygon

import sys

sys.path.append('User/admin/Desktop/OSSOS/MOP/src/ossos-pipeline/ossos')
from ossos import daophot
from ossos import storage
import sep_phot

_VOS_DIR = 'vos:kawebb/postage_stamps'
_DIR_PATH_BASE = os.path.dirname(os.path.abspath(__file__))
_STAMPS_DIR = '{}/postage_stamps'.format(_DIR_PATH_BASE)
_OSSOS_PATH_BASE = 'vos:OSSOS/dbimages'
_PHOT_DIR = '{}/phot_output'.format(_DIR_PATH_BASE)

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

    for i in range(5, 15):  # len(input)):

        detect_mbc('434', table['object'][i], table['expnum'][i], i)
        print '\n'


def detect_mbc(family_name, object_name, expnum, i):
    """
    Compare psf of asteroid with mean of stars to detect possible comae
    """
    phot_output_file = '{}/{}/{}_output_test.txt'.format(_PHOT_DIR, family_name, family_name)
    phot_output_table = pd.read_table(phot_output_file, sep=' ', dtype={'object': object})
    asteroid_id = phot_output_table.query('object == "{}" & expnum == "{}"'.format(object_name, expnum))
    print 'buffers: ', _BUFFER1, _BUFFER2
    print asteroid_id

    header, data = fits_data(object_name, expnum, family_name)
    bkg = get_bkg(data)

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

    ossos_path = '{}/{}/{}'.format(_OSSOS_PATH_BASE, expnum.strip('p'), header[_CCD])
    file_fwhm = '{}{}.fwhm'.format(expnum, header[_CCD].split('d')[1])
    storage.copy('{}/{}'.format(ossos_path, file_fwhm), '{}/{}'.format(_STAMPS_DIR, file_fwhm))
    with open('{}/{}'.format(_STAMPS_DIR, file_fwhm)) as infile:
        fwhm = float(infile.readline())
    os.unlink('{}/{}'.format(_STAMPS_DIR, file_fwhm))

    print '-- Getting asteroid data'
    ast_data, x_ast = get_asteroid_data(asteroid_id, data, i, fwhm, bkg)

    print '-- Getting star data'
    star_data, x_star = get_star_data(expnum, header, asteroid_id)

    print '-- Mean star fit'
    meas_psf(star_data, x_star, amp=100, sigma=5., alpha=1.)
    print '-- Asteroid fit'
    meas_psf(ast_data, x_ast, amp=100, sigma=12., alpha=10.)


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

                    return header, data

                    # os.unlink(file_path)


def get_star_data(expnum, header, object_data):
    """
    From ossos psf fitted image, calculate line profile
    """
    # calculate mean psf
    ossos_path = '{}/{}/{}'.format(_OSSOS_PATH_BASE, expnum.strip('p'), header[_CCD])
    file_psf = '{}{}.psf.fits'.format(expnum, header[_CCD].split('d')[1])
    storage.copy('{}/{}'.format(ossos_path, file_psf), '{}/{}'.format(_STAMPS_DIR, file_psf))
    local_file_path = '{}/{}'.format(_STAMPS_DIR, file_psf)
    local_psf = '{}{}.psf'.format(expnum, header[_CCD].split('d')[1])

    from pyraf import iraf

    iraf.set(uparm="./")
    iraf.digiphot(_doprint=0)
    iraf.apphot(_doprint=0)
    iraf.daophot(_doprint=0)
    iraf.seepsf(local_file_path, local_psf)

    with fits.open(local_file_path) as hdulist:
        data = hdulist[0].data

    os.unlink(local_psf)

    # data_masked = remove_saturated_rows(data, 117)
    data_masked = data

    th = math.degrees(object_data['theta'].values)
    data_rot = rotate(data_masked, th)
    np.savetxt('star_rot.txt', data_rot, fmt='%.0e')

    # sum all the values in each ROW
    totaled = []
    for row in range(len(data_rot)):
        sum_col = np.ma.sum(data_rot[row])
        if sum_col != 0.0:
            totaled.append(sum_col)

    # shift the baseline to zero
    zeroed = []
    for item in totaled:
        zeroed.append(item - np.amin(totaled))

    # normalize the data, perhaps not necessary?
    normed = []
    for item in zeroed:
        normed.append(item / np.amax(zeroed) * 100)

    midpt = int(len(normed) / 2)

    # remove core from psf
    x = []
    y = []
    for i in range(len(normed)):
        if (i < midpt - _BUFFER3) or (i > midpt + _BUFFER3):
            x.append(i)
            y.append(normed[i])

    return y, x


def get_asteroid_data(object_data, data, i, fwhm, bkg):
    """
    Calculate psf of asteroid, taking into acount trailing effect
    """

    data_cutout = cutout_data(object_data, data, fwhm, bkg)

    hdu = fits.PrimaryHDU()
    hdu.data = data_cutout
    hdu.writeto('cutout_-th_{}_fwhm.fits'.format(i), clobber=True)

    # sum all the values in each ROW
    totaled = []
    for row in range(len(data_cutout)):
        sum_col = np.ma.sum(data_cutout[row])
        if sum_col != 0.0:
            totaled.append(sum_col)

    # shift the baseline to zero
    zeroed = []
    for item in totaled:
        zeroed.append(item - np.amin(totaled))

    # normalize the data, perhaps not necessary?
    normed = []
    for item in zeroed:
        normed.append(item / np.amax(zeroed) * 100)

    x = np.array(range(len(normed)))

    return normed, x


def cutout_data(object_data, data, fwhm, bkg):
    """
    Cut a square around the polygon of the object, remove saturated rows, rotate so that the ellipse is parallel to the
    horizontal axis, then mask shape of ellipse around the object
    """

    # make sure there's only one object, otherwise pd.DataFrame.values won't work
    assert len(object_data) == 1

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
    np.copyto(data_obj, data[y_min:y_max, x_min:x_max])

    print 'saturation level, bkg: ', bkg * _SATURATION_LEVEL, bkg
    data_masked = remove_saturated_rows(data_obj, bkg * _SATURATION_LEVEL)

    # rotate the data about the angle of elongation, semi major axis is along x as reading out rows (not columsn)
    data_rot = pd.DataFrame(data=rotate(data_masked, th))

    rows, cols = data_rot.shape
    cx2 = int(cols / 2)
    cy2 = int(rows / 2)

    ell_buffer = _BUFFER1 * fwhm
    print 'buffer, fwhm: ', ell_buffer, fwhm

    '''
    # make a mask the shape of the ellipse THIS IS NOT RECCOMMENDED, NEED CORNER BACKGROUND FLUX
    y_range = range(cy2 - int(b + ell_buffer), cy2 + int(b + ell_buffer))
    mask = np.zeros((rows, cols))
    for y in y_range:
        x1 = int(cy2 - float(a + ell_buffer + 1) * (1. - (y - cy2) ** 2 / (b + ell_buffer + 1) ** 2) ** 0.5)
        x2 = int(cy2 + float(a + ell_buffer + 1) * (1. - (y - cy2) ** 2 / (b + ell_buffer + 1) ** 2) ** 0.5)
        for x in range(x1, x2):
            mask[y][x] = 1

    data_cutout = (np.multiply(data_rot, mask)).values
    '''

    # instead of ellipse, try cutting out rectangle
    x_min2 = cx2 - (a + ell_buffer)
    x_max2 = cx2 + (a + ell_buffer)
    y_min2 = cy2 - (b + ell_buffer)
    y_max2 = cy2 + (a + ell_buffer)
    data_cutout = (data_rot.values)[y_min2:y_max2, x_min2:x_max2]

    # np.savetxt('mask.txt', mask, fmt='%.0e')
    np.savetxt('data_cutout.txt', data_cutout)  # , fmt='%.0e')

    return data_cutout


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

        if np.ma.sum(sat_masked) > 3 * sat_level:
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
        bkg = (sep.Background(data2)).globalback  # , mask=mask, bw=64, bh=64, fw=3, fh=3) # optional parameters

    except ValueError:
        data3 = data2.byteswap(True).newbyteorder()
        bkg = (sep.Background(data3)).globalback  # , mask=mask, bw=64, bh=64, fw=3, fh=3) # optional parameters

    return bkg


def meas_psf(data_obj, x, amp, sigma, alpha):
    """
    Apply a fit to the point spread function and return parameters of fit
    """

    # calculate best parameters for the gaussian fit
    # x = data_obj.index
    mid_pt = np.argmax(data_obj)

    print x
    print len(x), len(data_obj)

    try:

        fitparams, fitcovariances = curve_fit(gauss, x, data_obj, p0=[amp, mid_pt, sigma])
        perr = np.sqrt(np.diag(fitcovariances))

        fitparams2, fitcovariances2 = curve_fit(moffat, x, data_obj, p0=[amp, mid_pt, sigma, alpha])
        perr2 = np.sqrt(np.diag(fitcovariances2))

        # fitparams3, fitcovariances3 = curve_fit(lorentz, x, data_obj, p0=[amp, mid_pt, 1.])
        # perr3 = np.sqrt(np.diag(fitcovariances3))

        # fitparams4, fitcovariances4 = curve_fit(penny, x, data_obj, p0=[1.5, mid_pt, 7., 10.])
        # perr4 = np.sqrt(np.diag(fitcovariances4))

        print fitparams, perr
        print fitparams2, perr2
        # print fitparams3, perr3
        # print fitparams4, perr4

        x2 = range(x[0], x[-1])

        gauss_fit = fitparams[0] * np.exp(-(x2 - fitparams[1]) ** 2 / (2. * fitparams[2] ** 2))

        moffat_fit = fitparams2[0] * (1 + (x2 - fitparams2[1]) ** 2 / fitparams2[2] ** 2) ** (
            -fitparams2[3])

        # lorentz_fit = (fitparams3[0] * np.ones(len(x2))) / (
        # 1 + (x2 / (fitparams3[1] * np.ones(len(x2)))) ** 2 + x2 * (fitparams3[2] * np.ones(len(x2))))

        # penny_fit = fitparams4[0] * ((1. - fitparams4[2]) / (1. + (x / fitparams4[1]) ** 2 + fitparams4[2] * np.exp(-0.63 * x ** 2 / fitparams4[1] ** 2 + x * fitparams4[3])))

        with sns.axes_style('ticks'):
            plt.plot(x, data_obj, label='Object psf')
            plt.plot(x2, gauss_fit, label='Gauss fit', ls='--')
            plt.plot(x2, moffat_fit, label='Moffat fit', ls=':')
            # plt.plot(x2, lorentz_fit, label='Lorentz fit', ls='-.')
            # plt.plot(x, penny_fit, label='Penny fit', ls='-.')
            plt.legend()
            plt.show()

    except Exception, e:
        print 'ERROR: ', e
        with sns.axes_style('ticks'):
            plt.plot(x, data_obj, label='Object psf')
            plt.legend()
            plt.show()


def gauss(x, *p):
    amp, mu, sigma = p
    return amp * np.exp(-(x - mu) ** 2 / (2. * sigma ** 2))


def moffat(x, *p):
    amp, mu, sigma, alpha = p
    return amp * (1. + (x - mu) ** 2 / (2. * sigma ** 2)) ** (-alpha)


def lorentz(x, *p):
    amp, p1, p2 = p
    return amp / (1. + (x ** 2 / p1 ** 2) + (x * p2))


def penny(x, *p):
    amp, p1, p2, p3 = p
    return amp * ((1. - p2) / (1. + (x / p1) ** 2) + p2 * np.exp(-0.693 * ((x / p1) ** 2 + x * p3)))


def compare_psf():
    """
    Compare psf of asteroid against mean of stars, check if anomaly in wings
    """


if __name__ == '__main__':
    main()