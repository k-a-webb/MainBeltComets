import argparse
from cStringIO import StringIO
import getpass
from astropy.io import fits
import requests
import os
import vos
from astropy.time import Time
import urllib2 as url
import numpy as np
from astropy.table import Table, Column
import pandas as pd

import sys
sys.path.append('User/admin/Desktop/OSSOS/MOP/src/ossos-pipeline/ossos')
from ossos import storage
from ossos import coding
from ossos import mpc
from ossos import util

_TARGET = "TARGET"
_DIR_PATH_BASE = os.path.dirname(os.path.abspath(__file__))
_IMAGE_LISTS = '{}/image_lists'.format(_DIR_PATH_BASE)
_VOS_PATH = 'vos:kawebb/postage_stamps'
_STAMPS_DIR = '{}/postage_stamps'.format(_DIR_PATH_BASE)

BASEURL = "http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/vospace/auth/synctrans"

"""
Retrieval of cutouts of the FITS images associated with the CFHT/MegaCam detections.
Takes a table (get_images.py output) as input
An Example URL for cutouts from OSSOS (not CFHT/MegaCam)
http://www.canfar.phys.uvic.ca/vospace/auth/synctrans?
TARGET=vos://cadc.nrc.ca~vospace/OSSOS/dbimages/1625356/1625356p.fits&
DIRECTION=pullFromVoSpace&
PROTOCOL=ivo://ivoa.net/vospace/core%23httpget&
view=cutout&
cutout=CIRCLE+ICRS+242.1318+-12.4747+0.05
"""


def main():
    parser = argparse.ArgumentParser(
        description='Parse an familyname.images.txt input file (get_images.py output) and create links in the postage stamp directory '
                    'that allow retrieval of cutouts of the FITS images associated with the CHT/MegaCam detections. '
                    'Cutouts are defined on the WCS ra/dec of the object position.')
    parser.add_argument("--family", '-f',
                        action="store",
                        default=None,
                        help="The input .txt files of astrometry/photometry measurements.")
    parser.add_argument("--radius", '-r',
                        action='store',
                        default=0.02,
                        help='Radius (degree) of circle of cutout postage stamp.')

    args = parser.parse_args()

    # CADC PERMISSIONS
    username = raw_input("CADC username: ")
    password = getpass.getpass("CADC password: ")

    get_stamps(str(args.family), username, password, args.radius)


def get_stamps(familyname, username, password, radius):
    """
    Cutout postage stamps of objects from VOSpace OSSOS mosaic, upload to VOSpace again
    """

    if familyname == 'none':
        vos_dir = '{}/none'.format(_VOS_PATH)
    else:
        vos_dir = '{}/all'.format(_VOS_PATH)

    print "----- Cutting postage stamps of objects in family {} -----".format(familyname)

    image_list_path = '{}/{}_images_new.txt'.format(_IMAGE_LISTS, familyname)
    try:
        table = pd.read_table(image_list_path, usecols=[0, 1, 3, 4], header=0, names=['object', 'expnum', 'ra', 'dec'],
                              sep=' ', dtype={'Object': object, 'Image': object})
    except pd.parser.CParserError:
        table = pd.read_table(image_list_path, usecols=[0, 1, 3, 4], header=0, names=['object', 'expnum', 'ra', 'dec'],
                              sep='\t', dtype={'Object': object, 'Image': object})

    for row in range(4669,len(table)):
        postage_stamp_filename = "{}_{}_{:8f}_{:8f}.fits".format(table['object'][row],
                                                                 table['expnum'][row],
                                                                 table['ra'][row],
                                                                 table['dec'][row])
        if storage.exists('{}/{}'.format(vos_dir, postage_stamp_filename)):
            print "  Stamp already exists"
        else:
            cutout(username, password, familyname, table['object'][row], table['expnum'][row], table['ra'][row],
                   table['dec'][row], radius)


def get_one_stamp(object_name, expnum, radius, username, password, familyname):
    """
    Pull stamp from VOSpace for a specific object in an exposure
    """

    print "-- Cutting postage stamps of {} {}".format(object_name, expnum)

    if familyname == 'none':
        vos_dir = '{}/none'.format(_VOS_PATH)
    else:
        vos_dir = '{}/all'.format(_VOS_PATH)

    try:
        image_list = '{}/{}_images.txt'.format(_IMAGE_LISTS, familyname)
    except:
        image_list, expnum_list, ra_list, dec_list = get_image_info(familyname, filtertype='r', imagetype='p')

    with open(image_list) as infile:
        for line in infile.readlines()[1:]:  # skip header info
            if len(line.split()) > 0:
                object_name_file = line.split()[0]
                expnum_file = line.split()[1]
                ra = float(line.split()[3])
                dec = float(line.split()[4])

                if (expnum == expnum_file) & (object_name == object_name_file):

                    postage_stamp_filename = "{}_{}_{:8f}_{:8f}.fits".format(object_name, expnum, ra, dec)
                    storage.remove('{}/{}'.format(vos_dir, postage_stamp_filename))

                    file_path = '{}/{}'.format(_STAMPS_DIR, postage_stamp_filename)
                    if os.path.exists(file_path):
                        os.remove(file_path)

                    cutout(username, password, familyname, object_name, expnum, ra, dec, radius)
                    return

    '''
    image_list_path = 'asteroid_families/{}/{}_images_test.txt'.format(familyname, familyname)
    table = pd.read_table(image_list_path, usecols=[0, 1, 3, 4], header=0, names=['object', 'expnum', 'ra', 'dec'], sep=' ', dtype={'object':object, 'expnum':object})
    
    row = table.query('(object == "{}") & (expnum == "{}")'.format(object_name, expnum))
    
    postage_stamp_filename = "{}_{}_{:8f}_{:8f}.fits".format(table['object'][row], table['expnum'][row], table['ra'][row], table['dec'][row])
    storage.remove('vos:kawebb/postage_stamps/{}/{}'.format(familyname, postage_stamp_filename))          
    cutout(username, password, familyname, table['object'][row], table['expnum'][row], table['ra'][row], table['dec'][row], radius)              
    '''


def centered_stamp(username, password, familyname, object_name, expnum, ra, dec, radius):

    init_dirs(familyname)

    postage_stamp_filename = "{}_{}_{:8f}_{:8f}_centered.fits".format(object_name, expnum, ra, dec)

    if storage.exists('{}/{}'.format(vos_dir, postage_stamp_filename)):
        print "  Stamp already exists"
    else:
        cutout(username, password, familyname, object_name, image, ra, dec, radius)
    return


def cutout(username, password, family_name, object_name, image, ra, dec, radius, test=False):

    """
    Test for image known to work   
    image = '1667879p'
    ra = 21.1236333333
    dec = 11.8697277778
    """

    if family_name == 'none':
        vos_dir = '{}/none'.format(_VOS_PATH)
    else:
        vos_dir = '{}/all'.format(_VOS_PATH)

    this_cutout = "CIRCLE ICRS {} {} {}".format(ra, dec, 2 * 0.18 / 3600.0)
    print "  cut out: {}".format(this_cutout)

    expnum = image.split('p')[0]  # only want calibrated images
    target = storage.vospace.fixURI(storage.get_uri(expnum))
    direction = "pullFromVoSpace"
    protocol = "ivo://ivoa.net/vospace/core#httpget"
    view = "cutout"
    params = {_TARGET: target,
              "PROTOCOL": protocol,
              "DIRECTION": direction,
              "cutout": this_cutout,
              "view": view}
    r = requests.get(BASEURL, params=params, auth=(username, password))

    try:
        r.raise_for_status()
        small_fobj = fits.open(StringIO(r.content))
        extname = small_fobj[0].header.get('EXTNAME', None)
        del small_fobj
    except requests.HTTPError, e:
        print 'Connection Failed, {}'.format(e)
        write_to_file(object_name, image)
        return

    this_cutout2 = "CIRCLE ICRS {} {} {}".format(ra, dec, radius)
    print "  cut out: {}".format(this_cutout2)

    expnum = image.split('p')[0]  # only want calibrated images
    target = storage.vospace.fixURI(storage.get_uri(expnum))
    direction = "pullFromVoSpace"
    protocol = "ivo://ivoa.net/vospace/core#httpget"
    view = "cutout"
    params = {"TARGET": target,
              "PROTOCOL": protocol,
              "DIRECTION": direction,
              "cutout": this_cutout2,
              "view": view}

    r = requests.get(BASEURL, params=params, auth=(username, password))

    try:
        # print "Getting extension: {}".format(extname)
        r.raise_for_status()
        full_fobj = fits.open(StringIO(r.content))
        if extname is not None:
            cutout_fobj = full_fobj[extname]
        else:
            cutout_fobj = full_fobj[0]
    except requests.HTTPError, e:
        print 'Connection Failed, {}'.format(e)
        write_to_file(object_name, image)
        return

    # print "Got this much data: {}".format(cutout_fobj.data.shape)
    cutout_fobj = fits.PrimaryHDU(data=cutout_fobj.data, header=cutout_fobj.header)

    postage_stamp_filename = "{}_{}_{:8f}_{:8f}.fits".format(object_name, image, float(ra), float(dec))
    cutout_fobj.writeto("{}/{}".format(_STAMPS_DIR, postage_stamp_filename), clobber=True)
    del cutout_fobj
    if test:
        return
    try:
        storage.copy('{}/{}'.format(_STAMPS_DIR, postage_stamp_filename),
                     '{}/{}'.format(vos_dir, postage_stamp_filename))
        os.unlink('{}/{}'.format(_STAMPS_DIR, postage_stamp_filename))
    except Exception, e:
        print e


def write_to_file(object_name, expnum):

    with open('{}/cutout_error.txt'.format(_IMAGE_LISTS), 'a') as outfile:
        outfile.write('{} {}\n'.format(object_name, expnum))

if __name__ == '__main__':
    main()