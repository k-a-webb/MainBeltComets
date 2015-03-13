import argparse
import os
import getpass
import pandas as pd

from find_family import find_family_members
from find_family import get_all_families_list
from get_images import get_image_info
from get_stamps import get_stamps, cutout
from sep_phot import iterate_thru_images
from ossos_scripts import storage

_LOCALPATH = 'asteroid_families/all/all_stamps'


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
                        default=0.02,
                        help='Radius (degree) of circle of cutout postage stamp.')
    parser.add_argument("--aperture", '-ap',
                        action='store',
                        default=10.0,
                        help='aperture (degree) of circle for photometry.')
    parser.add_argument("--thresh", '-th',
                        action='store',
                        default=5.0,
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
    vos_dir = 'vos:kawebb/postage_stamps/{}'.format(familyname)
    assert storage.exists(vos_dir)
    image_list_path = 'asteroid_families/{}/{}_images_test.txt'.format(familyname, familyname)  # USING TEST FILE
    print "WARNING: USING A TEST FILE ***************************************************************"

    '''
    # initiate output file
    out_filename = '{}_r{}_t{}_output.txt'.format(familyname, aperture, thresh)
    with open('asteroid_families/{}/{}_stamps/{}'.format(familyname, familyname, out_filename), 'w') as outfile:
        outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format('Object', "Image", 'RA', 'DEC', 'flux', 'mag', 'x', 'y'))
    '''

    # Remove any fits files hanging around from failed run
    for fits_file in os.listdir(_LOCALPATH):
        if fits_file.endswith('.fits'):
            storage.remove('{}/{}'.format(_LOCALPATH, fits_file))

    tkbad_list = []
    with open('asteroid_families/tkBAD.txt') as infile:
        for line in infile:
            tkbad_list.append(line.split(' ')[0])

    # for each image of each object, make cutout and go to sep_phot
    if os.path.exists(image_list_path):
        table = pd.read_table(image_list_path, usecols=[0, 1, 3, 4], header=0, names=['Object', 'Image', 'RA', 'DEC'],
                              sep=' ', dtype={'Object': object, 'Image': object})

        for row in range(11, 18):  # len(table)):
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
                                                  password, aperture, thresh, filtertype, imagetype)
                    attempts += 1

                    if attempts == 3:
                        print ' >>>> Last attempt'
                    print '\n'

    else:
        go_the_long_way(familyname, filtertype, imagetype)


def go_the_long_way(familyname, filtertype, imagetype):
    image_list, expnum_list, ra_list, dec_list = get_image_info(familyname, filtertype, imagetype)

    for index, objectname in enumerate(image_list):

        print '\n----- Searching for {} {} -----'.format(objectname, expnum_list[index])

        vos_dir = 'vos:kawebb/postage_stamps/{}'.format(familyname)
        postage_stamp_filename = "{}_{}_{:8f}_{:8f}.fits".format(objectname, expnum_list[index], ra_list[index],
                                                                 dec_list[index])
        if storage.exists('{}/{}'.format(vos_dir, postage_stamp_filename)):
            print "-- Stamp already exists"
        else:
            cutout(objectname, expnum_list[index], ra_list[index], dec_list[index], radius, username, password,
                   familyname)

        success = False
        attempts = 0
        while (success is False) and (attempts < 3):
            success = iterate_thru_images(familyname, str(table['Object'][row]), table['Image'][row], username,
                                          password, aperture, thresh, filtertype, imagetype)
            attempts += 1
            if attempts == 3:
                print ' >>>> Last attempt'
            print '\n'


if __name__ == '__main__':
    main()                        
