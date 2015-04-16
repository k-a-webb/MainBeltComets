
activate_this = '/Users/admin/Desktop/MainBeltComets/bin/activate_this.py'
execfile(activate_this, dict(__file__=activate_this))
import argparse
import os
import getpass
import pandas as pd

import sys
sys.path.append('User/admin/Desktop/OSSOS/MOP/src/ossos-pipeline/ossos')
from ossos import storage

from get_images import get_image_info
from get_stamps import cutout
from sep_phot import iterate_thru_images

_DIR_PATH_BASE = os.path.dirname(os.path.abspath(__file__))
_VOS_DIR = 'vos:kawebb/postage_stamps/'
_IMAGE_LISTS = '{}/image_lists'.format(_DIR_PATH_BASE)
_STAMPS_DIR = '{}/postage_stamps'.format(_DIR_PATH_BASE)
_PHOT_DIR = '{}/phot_output'.format(_DIR_PATH_BASE)


def main():
    """
    Input asteroid family name and an asteroid number and get out photometry values
    """

    parser = argparse.ArgumentParser(
        description='For an object in an asteroid family, parses AstDys for a list of members, \
                        parses MPC for images of those objects from CFHT/MegaCam in specific filter and exposure time,\
                        cuts out postage stamps images of given radius (should eventually be uncertainty ellipse), \
                         preforms photometry on a specified object given an aperture size and threshold, \
                        and then selects the object in the image from the predicted coordinates, magnitude, and eventually shape')

    parser.add_argument("--filter", '-fil',
                        action="store",
                        default='r',
                        dest="filter",
                        choices=['r', 'u'],
                        help="passband: default is r'")
    parser.add_argument("--family", '-f',
                        action="store",
                        default='all',
                        help='list of objects to query')
    parser.add_argument('--type',
                        default='p',
                        choices=['o', 'p', 's'],
                        help="restrict type of image (unprocessed, reduced, calibrated)")
    parser.add_argument("--radius", '-r',
                        action='store',
                        default=0.01,
                        help='Radius (degree) of circle of cutout postage stamp.')
    parser.add_argument("--aperture", '-ap',
                        action='store',
                        default=10.0,
                        help='aperture (degree) of circle for photometry.')
    parser.add_argument("--thresh", '-th',
                        action='store',
                        default=3.,
                        help='threshold value.')

    args = parser.parse_args()

    do_all_things(args.family, args.filter, args.type, float(args.radius), float(args.aperture),
                  float(args.thresh))


def do_all_things(familyname, filtertype, imagetype, radius, aperture, thresh):
    """
    For every postage stamps of an object in a family, determine orbital parameters of the object
    """

    # get CADC authentification
    username = raw_input("CADC username: ")
    password = getpass.getpass("CADC password: ")

    # establish input/output
    if familyname == 'none':
        vos_dir = '{}/none'.format(_VOS_DIR)
    else:
        vos_dir = '{}/all'.format(_VOS_DIR)
    assert storage.exists(vos_dir)
    image_list_path = '{}/{}_images.txt'.format(_IMAGE_LISTS, familyname)  # USING TEST FILE

    # Remove any fits files hanging around from failed run
    for fits_file in os.listdir(_STAMPS_DIR):
        if fits_file.endswith('.fits'):
            storage.remove('{}/{}'.format(_STAMPS_DIR, fits_file))

    tkbad_list = []
    with open('catalogue/tkBAD.txt') as infile:
        for line in infile:
            tkbad_list.append(line.split(' ')[0])

    # for each image of each object, make cutout and go to sep_phot
    if os.path.exists(image_list_path):

        try:
            table = pd.read_table(image_list_path, usecols=[0, 1, 3, 4], header=0,
                                  names=['Object', 'Image', 'RA', 'DEC'],
                                  sep=' ', dtype={'Object': object, 'Image': object})
        except pd.parser.CParserError:
            table = pd.read_table(image_list_path, usecols=[0, 1, 3, 4], header=0,
                                  names=['Object', 'Image', 'RA', 'DEC'],
                                  sep='\t', dtype={'Object': object, 'Image': object})

        for row in range(len(table)):
            print '\n{} --- Searching for {} in exposure {} -----'.format(row, table['Object'][row],
                                                                          table['Image'][row])
            expnum = (table['Image'][row]).strip('{}'.format(imagetype))
            if expnum in tkbad_list:
                print '-- Bad exposure'

            else:
                postage_stamp_filename = "{}_{}_{:8f}_{:8f}.fits".format(table['Object'][row], table['Image'][row],
                                                                         table['RA'][row], table['DEC'][row])

                if not storage.exists('{}/{}'.format(vos_dir, postage_stamp_filename)):
                    print '-- Cutout not found, creating new cutout'
                    cutout(username, password, familyname, table['Object'][row], table['Image'][row], table['RA'][row],
                           table['DEC'][row], radius)

                assert storage.exists('{}/{}'.format(vos_dir, postage_stamp_filename))
                success = False
                attempts = 0
                while (success is False) and (attempts < 3):
                    success = iterate_thru_images(familyname, str(table['Object'][row]), table['Image'][row], username,
                                                  password, aperture, thresh)
                    attempts += 1
                    if attempts == 3:
                        print ' >>>> Last attempt \n'
                        with open('{}/{}/vos_error.txt'.format(_PHOT_DIR, familyname), 'a') as outfile:
                            outfile.write('{} {}'.format(table['Object'][row], table['Image'][row]))

    else:
        print 'image list path does not exist'

if __name__ == '__main__':
    main()                        
