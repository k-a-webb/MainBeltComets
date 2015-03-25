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
    # meas_psf(expnum, header)

    meas_psf_asteroid(asteroid_id, data)


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
                    print file_path
                    data = hdulist[0].data
                    header = hdulist[0].header

                    return header, data

                    # os.unlink(file_path)


def meas_psf(expnum, header):
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


def meas_psf_asteroid(object_data, data):
    """
    Calculate psf of asteroid, taking into acount trailing effect
    """

    # make sure there's only one object, otherwise pd.DataFrame.values won't work
    assert len(object_data) == 1

    # define parameters of the square to cut out from the polygon around the elongated object
    a = object_data['a'].values
    b = object_data['b'].values
    th = object_data['theta'].values
    x = object_data['x'].values
    y = object_data['y'].values
    mag = object_data['mag'].values
    a_x = a * math.cos(th)
    a_y = a * math.sin(th)
    b_x = b * math.sin(th)
    b_y = b * math.cos(th)
    p1 = (x + b_x - a_x, y - a_y + b_y)
    p2 = (x + a_x + b_x, y + a_y + b_y)
    p3 = (x + a_x - b_x, y + a_y - b_y)
    p4 = (x - a_x - b_x, y - a_y - b_y)
    polygon = Polygon([p1, p2, p3, p4]).buffer(_BUFFER1)
    x_min, y_min, x_max, y_max = polygon.bounds

    # grab the fits data in a square around the polygon
    object_data = data[y_min:y_max, x_min:x_max]

    # rotate the data about the angle of elongation, and cut into a square again
    rot_od = pd.DataFrame(data=rotate(object_data, th))
    mid_col = int(len(rot_od.columns) / 2)
    mid_row = int(len(rot_od) / 2)
    rot_sq0 = rot_od.iloc[mid_row - int(a + _BUFFER2):mid_row + int(a + _BUFFER2), mid_col - int(b + _BUFFER2):mid_col + int(b + _BUFFER2)]
    rot_sq = rot_sq0.values

    # sum all the values in each column
    summed_cols = []
    for row in range(len(rot_sq)):
        sum_col = np.sum(rot_sq[row])
        if sum_col != 0.0:
            summed_cols.append(sum_col)

    zerod = []
    for item in summed_cols:
        zerod.append(item - np.amin(summed_cols))

    mid_pt = np.argmax(summed_cols)

    '''
    # shift the values so the zeropoint is 0
    shifted = []
    for item in summed_cols:
        shifted.append(item - np.amin(summed_cols))

    # normalize data to 1
    norm_cols = []
    for idx, item in enumerate(shifted):
        norm_cols.append((item / np.amax(shifted)))

    mid_pt = np.argmax(norm_cols)

    # build mask for saturated values
    mask = []

    for i in range(len(norm_cols)):
        if (i < int(mid_pt - b * 2)) or (i > int(mid_pt + b * 2)):
            mask.append(0)
        else:
            mask.append(1)
    norm_cols_masked = np.ma.masked_array(norm_cols, mask)
    '''

    if mag < 19:
        mask = []
        for i in range(len(summed_cols)):
            if (i < int(mid_pt - b * 2)) or (i > int(mid_pt + b * 2)):
                mask.append(0)
            else:
                mask.append(1)
        summed_cols = np.ma.masked_array(summed_cols, mask)

    # calculate best parameters for the gaussian fit
    p0 = [1.2, mid_pt, 1., 0.]
    p02 = [10000., mid_pt, 5., 1., 20000.]
    x = range(len(summed_cols))
    #fitparams, fitcovariances = curve_fit(gauss, x, summed_cols, p0=p0)
    fitparams2, fitcovariances2 = curve_fit(moffat, x, summed_cols, p0=p02)
    #print fitparams
    print fitparams2
    #gauss_fit = fitparams[0] * np.exp(-(x - fitparams[1]) ** 2 / (2. * fitparams[2] ** 2)) + fitparams[3]
    moffat_fit = fitparams2[0] * (1 + (x - fitparams2[1]) ** 2 / fitparams2[2] ** 2) ** (-fitparams2[3]) + fitparams2[4]

    with sns.axes_style('ticks'):
        plt.plot(x, summed_cols, label='Test data')
        #plt.plot(x, gauss_fit, label='Gauss fit')
        plt.plot(x, moffat_fit, label='Moffat fit')
        plt.legend()
        plt.show()


def gauss(x, *p):
    amp, mu, sigma, b = p
    return amp * np.exp(-(x - mu) ** 2 / (2. * sigma ** 2)) + b


def moffat(x, *p):
    amp, mu, sigma, alpha, b = p
    return amp * (1 + (x - mu) ** 2 / sigma ** 2) ** (-alpha) + b


def compare_psf():
    """
    Compare psf of asteroid against mean of stars, check if anomaly in wings
    """


if __name__ == '__main__':
    main()