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
from astropy.io import fits
from scipy.ndimage.interpolation import rotate
import math
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from shapely.geometry import Polygon
from collections import Counter

count = Counter()

import sep_phot
import ossos_scripts.wcs as wcs

import sys

sys.path.append('User/admin/Desktop/OSSOS/MOP/src/ossos-pipeline/ossos')
from ossos import daophot
from ossos import storage

sys.path.append('/Users/admin/Desktop/PythonPhot/PythonPhot')
from PythonPhot import getpsf, aper

_VOS_DIR = 'vos:kawebb/postage_stamps/all'
_DIR_PATH_BASE = os.path.dirname(os.path.abspath(__file__))
_DIR_PATH = '{}/asteroid_families'.format(_DIR_PATH_BASE)
_STAMPS_DIR = '{}/434/434_stamps'.format(_DIR_PATH)
_OSSOS_PATH_BASE = 'vos:OSSOS/dbimages'

_CCD = 'EXTNAME'
_ZMAG = 'PHOTZP'
_TARGET = 'OBJECT'

_APERTURE = 10.0
_THRESHOLD = 3.0
_BUFFER1 = 15.
_BUFFER2 = 10.
_SATURATION_LEVEL = 970
_SATURATION_THRESHOLD = 0.05  # percent fraction

client = vos.Client()


def main():
    parser = argparse.ArgumentParser(
        description='For a given set of fits images: preforms photometry, identifies a specific object, and returns \
                        the orbital elements of the object')
    parser.add_argument("--family", '-f',
                        action="store",
                        default='all',
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

    table = pd.read_table('{}/434/phot_output/434_output_test.txt'.format(_DIR_PATH), sep=' ', dtype={'object': object})

    for i in range(0, 1):  # len(input)):

        detect_mbc('434', table['object'][i], table['expnum'][i])


def detect_mbc(family_name, object_name, expnum):
    """
    Compare psf of asteroid with mean of stars to detect possible comae
    """
    phot_output_file = '{}/{}/phot_output/{}_output_test.txt'.format(_DIR_PATH, family_name, family_name)
    phot_output_table = pd.read_table(phot_output_file, sep=' ', dtype={'object': object})
    asteroid_id = phot_output_table.query('object == "{}" & expnum == "{}"'.format(object_name, expnum))
    print _BUFFER1, _BUFFER2
    print asteroid_id

    header, data = fits_data(object_name, expnum)

    '''
    target = header[_TARGET]
    if 'BH' not in target:
        print '-- NOT in H Block'
        return
    '''
    # get_star_data(expnum, header)
    ast_data = get_asteroid_data(asteroid_id, data, header)

    meas_psf(ast_data)


def fits_data(object_name, expnum):
    """
    Creates local copy of fits file from VOSpace
    """

    assert storage.exists(_VOS_DIR)
    for fits_file in client.listdir(_VOS_DIR):  # images named with convention: object_expnum_RA_DEC.fits

        if fits_file.endswith('.fits'):
            objectname_file = fits_file.split('_')[0]
            expnum_file = fits_file.split('_')[1]
            if (expnum_file == expnum) and (objectname_file == object_name):
                storage.copy('{}/{}'.format(_VOS_DIR, fits_file), '{}/{}'.format(_STAMPS_DIR, fits_file))
                file_path = '{}/{}'.format(_STAMPS_DIR, fits_file)

                with fits.open(file_path) as hdulist:
                    data = hdulist[0].data
                    header = hdulist[0].header

                    # Measure a spatially variable background of some image data (np array)
                    try:
                        bkg = sep.Background(data)  # , mask=mask, bw=64, bh=64, fw=3, fh=3) # optional parameters
                    except ValueError:
                        data = data.byteswap(True).newbyteorder()
                        bkg = sep.Background(data)  # , mask=mask, bw=64, bh=64, fw=3, fh=3) # optional parameters

                    # Directly subtract the background from the data in place
                    # bkg.subfrom(data)

                    return header, data

                    # os.unlink(file_path)


def get_star_data(expnum, header):
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

    data = fits.getdata('{}/{}.fits'.format(_STAMPS_DIR, local_psf))

    return data


def get_asteroid_data(object_data, data, header):
    """
    Calculate psf of asteroid, taking into acount trailing effect
    """

    # reject any object too bright that will definetly be saturated
    mag = object_data['mag'].values
    if mag < 15.5:
        print '>> Object is too bright for accurate photometry'

    data_cutout = cutout_data(object_data, data, header)

    # sum all the values in each column
    totaled = []
    for row in range(len(data_cutout)):

        sum_col = np.sum(data_cutout[row])
        if sum_col != 0.0:
            totaled.append(sum_col)

    # shift the baseline to zero
    zeroed = []
    for item in totaled:
        zeroed.append(item - np.amin(totaled))

    # normalize the data, perhaps not necessary?
    normed = []
    for item in zeroed:
        normed.append(item / np.amax(zeroed))

    return normed


def cutout_data(object_data, data, header):
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
    polygon = Polygon([p1, p2, p3, p4]).buffer(_BUFFER1)
    x_min, y_min, x_max, y_max = polygon.bounds
    data_obj = data[y_min:y_max, x_min:x_max]
    (pd.DataFrame(data=data_obj)).to_csv('rawdata.txt')

    # data_masked = build_mask(object_data, header)

    # rotate the data about the angle of elongation, and cut into a square again
    data_rot = pd.DataFrame(data=rotate(data_obj, th))
    print data_rot.shape
    rows, cols = data_rot.shape
    cx2 = int(cols / 2)
    cy2 = int(rows / 2)

    y_range = range(cy2 - int(b + _BUFFER1), cy2 + int(b + _BUFFER1))
    in_ell = []

    mask = np.zeros((rows, cols))
    for y in y_range:
        x1 = int(float(a + _BUFFER1 + 1) * (1. - (y - cy2) ** 2 / float(b + _BUFFER1 + 1) ** 2) ** 0.5 + cx2)
        x2 = int(-1 * float(a + _BUFFER1 + 1) * (1. - (y - cy2) ** 2 / float(b + _BUFFER1 + 1) ** 2) ** 0.5 + cx2)
        for x in range(x2, x1):
            mask[y - 1][x - 1] = 1

    data_cutout = np.ma.masked_array(data_rot, mask)
    np.savetxt('data_cutout.txt', data_cutout, fmt='%.1e')
    print data_cutout

    #hdu = fits.PrimaryHDU()
    #hdu.data = data_masked
    #hdu.writeto('cutout.fits', clobber=True)

    '''
    mid_col = int(len(data_rot.columns) / 2)
    mid_row = int(len(data_rot) / 2)
    data_rot_sq = data_rot.iloc[mid_row - int(a + _BUFFER2):mid_row + int(a + _BUFFER2),
                  mid_col - int(b + _BUFFER2):mid_col + int(b + _BUFFER2)]

    # data_rot_sq.to_csv('data.txt')
    '''

    return data_cutout


def build_mask(obj_data, header):
    """
    Check for saturation
    """
    # ccd_satur = header['SATURATE']
    satur_col = []
    for row in obj_data:
        data_row = obj_data[row]
        sat = data_row.where(data_row > _SATURATION_LEVEL)
        # sat_np = sat.values
        sat_arr = sat_np[~np.isnan(sat)]

        if (len(sat_arr) / len(obj_data)) > _SATURATION_THRESHOLD:
            satur_col.append(row)

    print satur_col

    return obj_data


def meas_psf(obj_data):
    """
    Apply a fit to the point spread function and return parameters of fit
    """

    # calculate best parameters for the gaussian fit
    x = range(len(obj_data))
    mid_pt = np.argmax(obj_data)
    fitparams, fitcovariances = curve_fit(gauss, x, obj_data, p0=[1.5, mid_pt, 10.])
    fitparams2, fitcovariances2 = curve_fit(moffat, x, obj_data, p0=[1.5, mid_pt, 10., 1.])
    fitparams3, fitcovariances3 = curve_fit(lorentz, x, obj_data, p0=[1.5, mid_pt, 1.])
    # fitparams4, fitcovariances4 = curve_fit(penny, x, obj_data, p0=[1.5, mid_pt, 7., 10.])
    print fitparams
    print fitparams2
    print fitparams3
    # print fitparams4
    gauss_fit = fitparams[0] * np.exp(-(x - fitparams[1]) ** 2 / (2. * fitparams[2] ** 2))  # + fitparams[3]
    moffat_fit = fitparams2[0] * (1 + (x - fitparams2[1]) ** 2 / fitparams2[2] ** 2) ** (
        -fitparams2[3])  # + fitparams2[4]
    lorentz_fit = (fitparams3[0] * np.ones(len(x))) / (
        1 + (x / (fitparams3[1] * np.ones(len(x)))) ** 2 + x * (fitparams3[2] * np.ones(len(x))))
    # penny_fit = fitparams4[0] * ((1. - fitparams4[2]) / (1. + (x / fitparams4[1]) ** 2 + fitparams4[2] * np.exp(-0.63 * x ** 2 / fitparams4[1] ** 2 + x * fitparams4[3])))

    with sns.axes_style('ticks'):
        plt.plot(x, obj_data, label='Test data')
        plt.plot(x, gauss_fit, label='Gauss fit', ls='--')
        plt.plot(x, moffat_fit, label='Moffat fit', ls=':')
        plt.plot(x, lorentz_fit, label='Lorentz fit', ls='-.')
        # plt.plot(x, penny_fit, label='Penny fit', ls='-.')
        plt.legend()
        plt.show()


def gauss(x, *p):
    amp, mu, sigma = p
    return amp * np.exp(-(x - mu) ** 2 / (2. * sigma ** 2))


def moffat(x, *p):
    amp, mu, sigma, alpha = p
    return amp * (1. + (x - mu) ** 2 / sigma ** 2) ** (-alpha)


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