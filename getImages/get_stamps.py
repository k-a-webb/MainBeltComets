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
sys.path.append('/Users/admin/Desktop/MainBeltComets/getImages/ossos_scripts/')

from ossos_scripts import storage
from ossos_scripts import coding
from ossos_scripts import mpc
from ossos_scripts import util

_TARGET = "TARGET"

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
                        default=0.01,
                        help='Radius (degree) of circle of cutout postage stamp.')
    parser.add_argument("--suffix", '-s',
                        action='store',
                        default=None,
                        help='Suffix of mba without family designation')
    args = parser.parse_args()
    
    # CADC PERMISSIONS
    username = raw_input("CADC username: ")
    password = getpass.getpass("CADC password: ")
    
    get_stamps(args.family, username, password, args.radius, args.suffix)
    
def get_stamps(familyname, username, password, radius=0.01, suffix=None):
    
    print "----- Cutting postage stamps of objects in family {}  from CFHT/MegaCam images -----".format(familyname)	                  
    
    dir_path = os.path.dirname(os.path.abspath(__file__))
    dir_path_base = '{}/asteroid_families'.format(dir_path)
    family_dir = os.path.join(dir_path_base, familyname)
    if os.path.isdir(family_dir) == False:
        print "Invalid family name or directory does not exist"
        
    vos_dir = 'vos:kawebb/postage_stamps/{}'.format(familyname)
    if not storage.exists(vos_dir, force=True):
        storage.mkdir(vos_dir)
    assert storage.exists(vos_dir, force=True)    

    image_list_path = 'asteroid_families/{}/{}_images.txt'.format(familyname, familyname) # USING TEST FILE
    table = pd.read_table(image_list_path, usecols=[0, 1, 3, 4], header=0, names=['object', 'expnum', 'ra', 'dec'], sep=' ', dtype={'Object':object})
    
    for row in range(len(table)):
        postage_stamp_filename = "{}_{}_{:8f}_{:8f}.fits".format(table['object'][row], table['expnum'][row], table['ra'][row], table['dec'][row])
        if storage.exists('{}/{}'.format(vos_dir, postage_stamp_filename)) == True:
            print "  Stamp already exists"
        else:                
            cutout(username, password, familyname, table['object'][row], table['expnum'][row], table['ra'][row], table['dec'][row], radius=0.02)
        
                
def get_one_stamp(objectname, expnum, radius, username, password, familyname):
    
    print "-- Cutting postage stamps of {} {}".format(objectname, expnum)	                  
    
    dir_path = os.path.dirname(os.path.abspath(__file__))
    dir_path_base = '{}/asteroid_families'.format(dir_path)
    family_dir = os.path.join(dir_path_base, familyname)
    if os.path.isdir(family_dir) == False:
        print "Invalid family name or directory does not exist"
        
    vos_dir = 'vos:kawebb/postage_stamps/{}'.format(familyname)
    if not storage.exists(vos_dir, force=True):
        storage.mkdir(vos_dir)
    
    image_list = '{}/{}_images.txt'.format(family_dir, familyname)
    with open(image_list) as infile: 
        for line in infile.readlines()[1:]: # skip header info
            assert len(line.split()) > 0
            objectname = line.split()[0]
            expnum_file = line.split()[1]
            ra = float(line.split()[3])   
            dec = float(line.split()[4])
                    
            if expnum == expnum_file:
            
                postage_stamp_filename = "{}_{}_{:8f}_{:8f}.fits".format(objectname, expnum, ra, dec)
                storage.remove('vos:kawebb/postage_stamps/{}/{}'.format(familyname, postage_stamp_filename))
                file_path = '{}/{}_stamps/{}'.format(family_dir, familyname, postage_stamp_filename)
                if os.path.exists(file_path):
                    os.remove(file_path)
                cutout(username, password, familyname, objectname, expnum, ra, dec, radius)
                return    
    
    '''
    image_list_path = 'asteroid_families/{}/{}_images_test.txt'.format(familyname, familyname)
    table = pd.read_table(image_list_path, usecols=[0, 1, 3, 4], header=0, names=['object', 'expnum', 'ra', 'dec'], sep=' ', dtype={'object':object, 'expnum':object})
    
    row = table.query('(object == "{}") & (expnum == "{}")'.format(objectname, expnum))
    
    postage_stamp_filename = "{}_{}_{:8f}_{:8f}.fits".format(table['object'][row], table['expnum'][row], table['ra'][row], table['dec'][row])
    storage.remove('vos:kawebb/postage_stamps/{}/{}'.format(familyname, postage_stamp_filename))          
    cutout(username, password, familyname, table['object'][row], table['expnum'][row], table['ra'][row], table['dec'][row], radius)              
	'''
def centered_stamp(username, password, familyname, objectname, expnum, ra, dec, radius):
    
    vos_dir = 'vos:kawebb/postage_stamps/{}'.format(familyname)
    postage_stamp_filename = "{}_{}_{:8f}_{:8f}_centered.fits".format(objectname, expnum, ra, dec)
    if storage.exists('{}/{}'.format(vos_dir, postage_stamp_filename)) == True:
        print "  Stamp already exists"
    else:
        print type(objectname, expnum, ra, dec, radius, username, password, familyname)
        cutout(username, password, familyname, objectname, image, ra, dec, radius)
    return


def cutout(username, password, family_name, object_name, image, ra, dec, radius, test=False):
    """
    Test for image known to work   
    image = '1667879p'
    ra = 21.1236333333
    dec = 11.8697277778
    """

    vos_dir = 'vos:kawebb/postage_stamps/{}'.format(family_name)
    output_dir = 'asteroid_families/{}/{}_stamps'.format(family_name, family_name)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    this_cutout = "CIRCLE ICRS {} {} {}".format(ra, dec, 2*0.18/3600.0)
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
        return

    this_cutout = "CIRCLE ICRS {} {} {}".format(ra, dec, radius)
    print "  cut out: {}".format(this_cutout)

    expnum = image.split('p')[0]  # only want calibrated images
    target = storage.vospace.fixURI(storage.get_uri(expnum))
    direction = "pullFromVoSpace"
    protocol = "ivo://ivoa.net/vospace/core#httpget"
    view = "cutout"
    params = {"TARGET": target,
              "PROTOCOL": protocol,
              "DIRECTION": direction,
              "cutout": this_cutout,
              "view": view}

    r = requests.get(BASEURL, params=params, auth=(username, password))

    try:
        r.raise_for_status()
        full_fobj = fits.open(StringIO(r.content))
        if extname is not None:
            cutout_fobj = full_fobj[extname]
        else:
            cutout_fobj = full_fobj[0]
    except requests.HTTPError, e:
        print 'Connection Failed, {}'.format(e)
        return

    cutout_fobj = fits.PrimaryHDU(data=cutout_fobj.data, header=cutout_fobj.header)

    postage_stamp_filename = "{}_{}_{:8f}_{:8f}.fits".format(object_name, image, float(ra), float(dec))
    cutout_fobj.writeto("{}/{}".format(output_dir, postage_stamp_filename))
    del cutout_fobj
    if test:
        return
    storage.copy('{}/{}'.format(output_dir, postage_stamp_filename),
                 '{}/{}'.format(vos_dir, postage_stamp_filename))
    os.unlink('{}/{}'.format(output_dir, postage_stamp_filename))
    	
if __name__ == '__main__':
    main()	
		
		
		
		
		