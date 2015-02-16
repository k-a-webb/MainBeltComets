"""Retrieval of cutouts of the FITS images associated with the CFHT/MegaCam detections.
Takes a table (getmbc.py output) as input

An Example URL for cutouts from OSSOS (not CFHT/MegaCam)
http://www.canfar.phys.uvic.ca/vospace/auth/synctrans?
TARGET=vos://cadc.nrc.ca~vospace/OSSOS/dbimages/1625356/1625356p.fits&
DIRECTION=pullFromVoSpace&
PROTOCOL=ivo://ivoa.net/vospace/core%23httpget&
view=cutout&
cutout=CIRCLE+ICRS+242.1318+-12.4747+0.05
"""

import argparse
import getpass
import requests
import os
import vos

import sys
sys.path.append('/Users/admin/Desktop/MainBeltComets/getImages/ossos_scripts/')

from ossos_scripts import storage
from ossos_scripts import coding
from ossos_scripts import mpc
from ossos_scripts import util

import numpy as np
from astropy.table import Table, Column

BASEURL = "http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/vospace/auth/synctrans"


    
def get_stamps(familyname, radius):
    
    print "----- Cutting postage stamps of objects in family {}  from CFHT/MegaCam images -----".format(familyname)	

    # CADC PERMISSIONS
    username = raw_input("CADC username: ")
    password = getpass.getpass("CADC password: ")

    # PARSE THROUGH INPUT FILE FOR IMAGE INFORMATION
        # format into lines, parse for image, RA and DEC
        # CUT OUT IMAGE
	        # PASS: object, directory, radius, CADC permissions                    
    
    dir_path_base = '/Users/admin/Desktop/MainBeltComets/getImages/asteroid_families'
    family_dir = os.path.join(dir_path_base, familyname)
    if os.path.isdir(family_dir) == False:
        print "Invalid family name or directory does not exist"
    
    image_list = '{}/{}_images.txt'.format(family_dir, familyname)
    
    with open(image_list) as infile: 
        for line in infile.readlines()[1:]: # skip header info
            assert len(line.split()) > 0
            objectname = line.split()[0]
            expnum = line.split()[1]
            RA = float(line.split()[3])   
            DEC = float(line.split()[4])
            
            vos_dir = 'vos:kawebb/postage_stamps/{}'.format(familyname)
            if not storage.exists(vos_dir, force=True):
                storage.mkdir(vos_dir)
            #assert storage.exists(vos_dir, force=True)

            cutout(objectname, expnum, RA, DEC, radius, username, password, familyname, vos_dir)
	
def cutout(objectname, image, RA, DEC, radius, username, password, familyname, vos_dir):
    # CUT OUT (image, RA, DEC, radius, CADC permissions)
	# for each attribute in mpc_observations:
		# storage.vospace.fixURI and 
			# storage.get_uri (Build the uri for an OSSOS image stored in the dbimages containerNode)
		# define parameters for request: target, protocol, direction, cutout, view
		# request image : url, parameters, cadc permissions
		# assign to a file
    
    ''' Test for image known to work   
    image = '1667879p'
    RA = 21.1236333333
    DEC = 11.8697277778
    '''
    
    output_dir = 'asteroid_families/{}/{}_stamps'.format(familyname, familyname)
    if os.path.isdir(output_dir) == False:
            os.makedirs(output_dir)
        
    this_cutout = "CIRCLE ICRS {} {} {}".format(RA, DEC, radius)                                 
    print "cut out: ", this_cutout

    expnum = image.split('p')[0] # only want calibrated images    
    target = storage.vospace.fixURI(storage.get_uri(expnum)) 
    direction = "pullFromVoSpace"
    protocol = "ivo://ivoa.net/vospace/core#httpget"
    view = "cutout"
    params = {"TARGET": target,
                  "PROTOCOL": protocol,
                  "DIRECTION": direction,
                  "cutout": this_cutout,
                  "view": view}
    
    try:
        r = requests.get(BASEURL, params=params, auth=(username, password))
        
        postage_stamp_filename = "{}_{}_{:8f}_{:8f}.fits".format(objectname, image, RA, DEC)
    
        with open('{}/{}'.format(output_dir, postage_stamp_filename), 'w') as tmp_file:
            object_dir = 'asteroid_families/{}/{}_stamps/{}'.format(familyname, familyname, postage_stamp_filename)
            assert os.path.exists(object_dir)
            tmp_file.write(r.content)
            storage.copy(object_dir, '{}/{}'.format(vos_dir, postage_stamp_filename))
        os.unlink(object_dir)  # easier not to have them hanging around
        
    except requests.exceptions.ConnectionError as e:
        print "These aren't the domains we're looking for."
    
    #r = requests.get(BASEURL, params=params, auth=(username, password))
    #r.raise_for_status()  # confirm the connection worked as hoped
    
    
def main():
    
    # PARSE INFORMATION INPUTTED FROM THE COMMAND LINE
	# VERSION - OSSOS DATA RELEASE VERSION THE STAMPS ARE TO BE ASSIGNED TO
	# INPUT FILE
	# RADIUS - SIZE OF CIRCULAR CUTOUT
    
    # INPUT LIST OF IMAGES IN COMMAND LINE
    # IDENTIFY PARAMETERS FOR QUERY OF SSOIS FROM INPUT

    parser = argparse.ArgumentParser(
        description='Parse an obects.txt input file (find_family.py output) and create links in the postage stamp directory '
                    'that allow retrieval of cutouts of the FITS images associated with the CHT/MegaCam detections. '
                    'Cutouts are defined on the WCS RA/DEC of the object position.')
    parser.add_argument("--family", '-f',
                        action="store",
                        default="lixImages.txt",
                        help="The input .txt files of astrometry/photometry measurements.")
    parser.add_argument("--radius", '-r',
                        action='store',
                        default=0.005,
                        help='Radius (degree) of circle of cutout postage stamp.')
    
    args = parser.parse_args()
    
    get_stamps(args.family, args.radius)
		
if __name__ == '__main__':
    main()	
		
		
		
		
		