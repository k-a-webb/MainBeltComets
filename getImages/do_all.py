import argparse
import os
import getpass

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
                    default=0.01,
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

    do_all_things(args.family, args.object, args.filter, args.type, args.radius, args.aperture, args.thresh)
    
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
    
def do_all_things(familyname, objectname=None, filtertype='r', imagetype='p', radius=0.01, aperture=10.0, thresh=3.5):
   
    family_list_path = 'asteroid_families/{}/{}_family.txt'.format(familyname, familyname)
    
    if  os.path.exists(family_list_path):
        print "----- List of objects in family {} exists already -----".format(familyname)
        with open(family_list_path) as infile:
            filestr = infile.read()
            all_object_list = filestr.split('\n')
    else:    
        all_object_list = find_family_members(familyname)
    
    image_list_path = 'asteroid_families/{}/{}_images_test.txt'.format(familyname, familyname) # USING TEST FILE
    print "WARNING: USING A TEST FILE ***************************************************************" 
    if  os.path.exists(image_list_path):
        expnum_list = []
        image_list = []
        ra_list = []
        dec_list = []
        print "----- List of images in family {} exists already -----".format(familyname)
        with open(image_list_path) as infile:
            filestr = infile.read()
            fileline = filestr.split('\n')
            for item in fileline[1:]:
                if len(item.split()) > 0:
                    image_list.append(item.split()[0])
                    expnum_list.append(item.split()[1])
                    ra_list.append(float(item.split()[3]))
                    dec_list.append(float(item.split()[4]))
    else:  
        image_list, expnum_list, ra_list, dec_list = get_image_info(familyname, filtertype, imagetype) 

    username = raw_input("CADC username: ")
    password = getpass.getpass("CADC password: ")
    
    out_filename = '{}_r{}_t{}_output.txt'.format(familyname, aperture, thresh)
    with open('asteroid_families/{}/{}_stamps/{}'.format(familyname, familyname, out_filename), 'w') as outfile:
        outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format('Object', "Image", 'flux', 'mag', 'RA', 'DEC', 'ecc', 'index'))
    
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
    