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
                    default=5.0,
                    help='threshold value.')
    parser.add_argument("--object", '-obj',
                    action='store',
                    default='54286',
                    help='the object to preform photometry on')
                            
    args = parser.parse_args()

    do_all_things(args.family, args.object, args.filter, args.type, args.radius, args.aperture, args.thresh)

def do_all_things(familyname, objectname=None, filtertype='r', imagetype='p', radius=0.01, aperture=10.0, thresh=5.0):
   
    username = raw_input("CADC username: ")
    password = getpass.getpass("CADC password: ")
    
    family_list_path = 'asteroid_families/{}/{}_family.txt'.format(familyname, familyname)
    
    if  os.path.exists(family_list_path):
        print "----- List of objects in family {} exists already -----".format(familyname)
        with open(family_list_path) as infile:
            filestr = infile.read()
            all_object_list = filestr.split('\n')
    else:    
        all_object_list = find_family_members(familyname)
    
    
    out_filename = '{}_r{}_t{}_output.txt'.format(familyname, aperture, thresh)
    with open('asteroid_families/{}/{}_stamps/{}'.format(familyname, familyname, out_filename), 'w') as outfile:
        outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format('Object', "Image", 'flux', 'mag', 'RA', 'DEC', 'ecc', 'index'))
    

    image_list_path = 'asteroid_families/{}/{}_images_test.txt'.format(familyname, familyname) # USING TEST FILE
    print "WARNING: USING A TEST FILE ***************************************************************" 
    if  os.path.exists(image_list_path):
        table = pd.read_table(image_list_path, usecols=[0, 1, 3, 4], header=0, names=['Object', 'Image', 'RA', 'DEC'], sep=' ', dtype={'Object':object})
        for row in range(0,4):#len(table)):
            print '\n----- Searching for {} {} -----'.format(table['Object'][row], table['Image'][row])
            vos_dir = 'vos:kawebb/postage_stamps/{}'.format(familyname)
            postage_stamp_filename = "{}_{}_{:8f}_{:8f}.fits".format(table['Object'][row], table['Image'][row], table['RA'][row], table['DEC'][row])
            if storage.exists('{}/{}'.format(vos_dir, postage_stamp_filename)) == True:
                print "-- Stamp already exists"
            else:
                cutout(table['Object'][row], table['Image'][row], table['RA'][row], table['DEC'][row], radius, username, password, familyname)
                
            object_data = iterate_thru_images(familyname, str(table['Object'][row]), table['Image'][row], username, password, aperture, thresh, filtertype, imagetype)
    else:  
        go_the_long_way(familyname, filtertype, imagetype)

def go_the_long_way(familyname, filtertype, imagetype):
        
    image_list, expnum_list, ra_list, dec_list = get_image_info(familyname, filtertype, imagetype) 
    
    for index, objectname in enumerate(image_list):
        
        print '\n----- Searching for {} {} -----'.format(objectname, expnum_list[index])
        
        vos_dir = 'vos:kawebb/postage_stamps/{}'.format(familyname)
        postage_stamp_filename = "{}_{}_{:8f}_{:8f}.fits".format(objectname, expnum_list[index], ra_list[index], dec_list[index])
        if storage.exists('{}/{}'.format(vos_dir, postage_stamp_filename)) == True:
            print "-- Stamp already exists"
        else:
            cutout(objectname, expnum_list[index], ra_list[index], dec_list[index], radius, username, password, familyname)
                
        object_data = iterate_thru_images(familyname, objectname, expnum_list[index], username, password, aperture, thresh, filtertype, imagetype)  
                        
if __name__ == '__main__':
    main()                        
    