__author__ = 'kw'

import pandas as pd
import numpy as np
import vos
import argparse

import sep_phot
import ossos_scripts.wcs as wcs

_VOS_PATH = 'vos:kawebb/postage_stamps/all'
_DIR_PATH_BASE = os.path.dirname(os.path.abspath(__file__))
_DIR_PATH = '{}/asteroid_families'.format(_DIR_PATH_BASE)
_APERTURE = 10.0
_THRESHOLD = 3.0
_MAG_LIM = 20.0 # minimum brightness of catalog star to compare psf with



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
                        default=None,
                        help='The expsure in which the object is being identified.')

    args = parser.parse_args()


    iden_catalog_stars(args.family, args.object, args.expnum)

def iden_catalog_stars(family_name, object_name, expnum):
    """
    Compare photometry output with Stephen's catalogue to identify stars to use for PSF comparison
    """
    phot_output_table = '{}/{}/phot_output/{}_output_test'.format(_DIR_PATH, family_name, family_name)
    phot_output = pd.read_table(phot_output_table, sep='\t', dtype={'object': object})
    asteroid_id = phot_output.query('object == "{}" & expnum == "{}"'.format(object_name, expnum))

    septable, header = sep_phot.get_fits_data(object_name, expnum, _APERTURE, _THRESHOLD)
    pvwcs = wcs.WCS(header)
    table = sep_phot.append_table(septable, pvwcs, header['PHOTZP'], asteroid_id['ra'][0], asteroid_id['dec'][0])
    transients, catalogue, num_cat_objs, catalog_stars = sep_phot.compare_to_catalogue(table)

    bright_stars = catalog_stars.query('mag < {}'.format(_MAG_LIM))


    # calculate mean psf


def meas_psf():
    """
    Calculate psf of stars (and asteroid?)
    """


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