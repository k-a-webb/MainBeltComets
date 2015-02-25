import argparse
import os
import getpass

from find_family import find_family_members
from find_family import get_all_families_list
from get_images import get_image_info
from get_stamps import get_stamps
from sep_phot import select_object

from ossos_scripts import storage
import ossos_scripts.coding
import ossos_scripts.mpc
import ossos_scripts.util

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

    parser.add_argument("--filter",
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
                    default=0.005,
                    help='Radius (degree) of circle of cutout postage stamp.')
    parser.add_argument("--aperture", '-a',
                    action='store',
                    default=10.0,
                    help='aperture (degree) of circle for photometry.')
    parser.add_argument("--thresh", '-t',
                    action='store',
                    default=5.0,
                    help='threshold value.')
    parser.add_argument("--object", '-o',
                    action='store',
                    default='54286',
                    help='the object to preform photometry on')
                            
    args = parser.parse_args()
    
    # CADC PERMISSIONS
    username = raw_input("CADC username: ")
    password = getpass.getpass("CADC password: ")
    
    do_all_things(username, password, args.family, args.object, args.filter, args.type, args.radius, args.aperture, args.thresh)
    
    '''
    family_file = 'asteroid_families/{}/{}_family.txt'.format(args.family, args.family)
    if os.path.exists(family_file):
        with open(family_file) as infile:
            for line in infile:
                do_all_things

    
        for familyname in families_list[32:]:
            do_all_things(username, password, familyname, args.object, args.filter, args.type, args.radius, args.aperture, args.thresh)
    
    else:
        family_list = find_family_members(args.family)
        do_all_things(username, password, args.family, args.object, args.filter, args.type, args.radius, args.aperture, args.thresh)
    '''
    
def do_all_things(username, password, familyname, objectname=None, filtertype='r', imagetype='p', radius=0.005, aperture=10.0, thresh=5.0):
   
    family_list_path = 'asteroid_families/{}/{}_family.txt'.format(familyname, familyname) # TEST FILE
    #print "WARNING: USING A TEST FILE ******************************"
    if  os.path.exists(family_list_path):
        print "----- List of objects in family {} exists already -----".format(familyname)
        with open(family_list_path) as infile:
            filestr = infile.read()
            all_object_list = filestr.split('\n')
    else:    
        all_object_list = find_family_members(familyname)
    
    image_list = []
    image_list_path = 'asteroid_families/{}/{}_images.txt'.format(familyname, familyname)  
    if  os.path.exists(image_list_path):
        print "----- List of images in family {} exists already -----".format(familyname)
        with open(image_list_path) as infile:
            filestr = infile.read()
            fileline = filestr.split('\n')
            for item in fileline:
                if len(item.split()) > 0:
                    image_list.append(item.split()[0])
    else:  
        image_list = get_image_info(familyname, filtertype, imagetype) 
                
    object_list = []
    for item in all_object_list:
        if item in image_list:
            object_list.append(item)

    '''
    stamps_dir = 'asteroid_families/{}/{}_stamps'.format(familyname, familyname)
    if os.path.exists(stamps_dir):
        print '----- Stamps already exist in VOSpace -----'
    else:
        get_stamps(familyname, radius, username, password)
    '''
    
    for objectname in object_list:
        select_object(familyname, objectname, aperture, thresh)    
              

                        
if __name__ == '__main__':
    main()                        
    