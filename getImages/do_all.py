import argparse
import os

from find_family import find_family_members
from find_family import get_all_families_list
from get_images import get_image_info
from get_stamps import get_stamps
from sep_phot import find_objects_by_phot

from ossos_scripts import storage
import ossos_scripts.coding
import ossos_scripts.mpc
import ossos_scripts.util



def do_all_things(familyname, objectname=None, filtertype='r', imagetype='p', radius=0.005, aperture=10.0, thresh=5.0):
   
    family_list_path = 'asteroid_families/{}/{}_family.txt'.format(familyname, familyname)
    if  os.path.exists(family_list_path):
        print "----- List of objects in family {} exists already -----".format(familyname)
        with open(family_list_path) as infile:
            filestr = infile.read()
            asteroid_list = filestr.split('\n')
    else:    
        asteroid_list = find_family_members(familyname)

    '''image_list_path = 'asteroid_families/{}/{}_images.txt'.format(familyname, familyname)  
    if  os.path.exists(image_list_path):
        print "----- List of images in family {} exists already -----".format(familyname)
        with open(image_list_path) as infile:
            filestr = infile.read()
            image_list = filestr.split('\n')
    else:'''    
    image_list = get_image_info(familyname, filtertype, imagetype)       
    
    if objectname == None:
        objectname = image_list[0]

    get_stamps(familyname, radius)
                
    #find_objects_by_phot(familyname, objectname, aperture, thresh)    
              
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
                    default=None,
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
    
    allfamily_list_path = 'asteroid_families/all_families.txt'
    if args.family == None: 
        if os.path.exists(allfamily_list_path):
            print "----- List of all family names exists already -----"
            with open(allfamily_list_path) as infile:
                filestr = infile.read()
                families_list = filestr.split('\n')
        else:    
            families_list = get_all_families_list()
    
        for familyname in families_list:
            do_all_things(familyname, args.object, args.filter, args.type, args.radius, args.aperture, args.thresh)
    
    else:
        do_all_things(args.family, args.object, args.filter, args.type, args.radius, args.aperture, args.thresh)
                        
if __name__ == '__main__':
    main()                        
    