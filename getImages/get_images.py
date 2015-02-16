import datetime
import os
from astropy.io import ascii
from astropy.time import Time
import requests
import argparse
import requests

from ossos_scripts.ssos import Query

def get_image_info(familyname, filtertype='r', imagetype='p'):
    '''
    Query the ssois ephemeris for images of objects in a given family. Then parse through for desired image type, 
    filter, exposure time, and telescopy instrument
    '''
    
    dir_path_base = '/Users/admin/Desktop/MainBeltComets/getImages/asteroid_families'
    family_dir = os.path.join(dir_path_base, familyname)
    if os.path.isdir(family_dir) == False:
        print "Invalid family name or directory does not exist"
    
    family_list = '{}/{}_family.txt'.format(family_dir, familyname)
    
    with open(family_list) as infile: 
        filestr = infile.read()
    object_list = filestr.split('\n') # array of objects to query

    # TO DO: confirm that the input is the proper format to search for the appropriate ephemeris
    
    # From the given input, identify the desired filter and rename appropriately
    if filtertype.lower().__contains__('r'):
        filtertype = 'r.MP9601'  # this is the old (standard) r filter for MegaCam
    if filtertype.lower().__contains__('u'):
        filtertype = 'u.MP9301'

    # Define time period of image search, basically while MegaCam in operation
    search_start_date=Time('2013-01-01', scale='utc')   # epoch1=2013+01+01
    search_end_date=Time('2017-01-01', scale='utc')     # epoch2=2017+1+1
        
    # Setup output, label columns
    with open('{}/{}_images.txt'.format(family_dir, familyname), 'w') as outfile:
        outfile.write("{:>10s} {:>10s} {:>10s} {:>16s} {:>16s} {:>16s} {:>12s}\n".format(
            "Object", "Image", "Exp_time", "RA", "DEC", "time", "filter"))

    # Query a specific ephemeris depending on object naming convention. CADC : search = bynameCADC , MPC : search = bynameMPC
    # Specify search parameters
        # ephemeris, name, date range, resolve image extension, resolve to x,y, positional uncertainty?
        # http:// www3.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/cadcbin/ssos/ssosclf.pl?lang=en; object=elst-pizarro; search=bynameMPC; epoch1=2013+01+01; epoch2=2015+1+16; eellipse=; eunits=arcseconds; extres=yes; xyres=yes; format=tsv
        
    print "-------------------- \n Searching for images of objects in family {} from CFHT/Megacam from the MPC ephemeris".format(familyname)
    print " with filter {} and exposure time of 287, 387, 500 seconds (OSSOS data) \n--------------------".format(filtertype)
        
    image_list = []
    for object_name in object_list[:len(object_list)-1]:
        query = Query(object_name, search_start_date=search_start_date, search_end_date=search_end_date)

        # GET/REVIEW THE DATA RETURNED FROM THE SEARCH
        # Parse the data for Image, Instrument, Filter and create table for each object
        # Download files with appropriate names? ie object+SSOISfilename ??
        obs_in_filter = parse_ssois_return(query.get(), object_name, imagetype, camera_filter=filtertype)

        # output the data is previously initiated output file
        if len(obs_in_filter) > 0:
            with open('{}/{}_images.txt'.format(family_dir, familyname), 'a') as outfile:
                for line in obs_in_filter:
                    image_list.append(object_name)
                    try:
                        outfile.write("{} {} {} {} {} {} {}\n".format(object_name,
                            line['Image'], line['Exptime'], line['Object_RA'], line['Object_Dec'],
                            Time(line['MJD'], format='mjd', scale='utc'), line['Filter']))
                    except:
                        print "cannot write to outfile"
    return image_list
                    
def parse_ssois_return(ssois_return, object_name, imagetype, camera_filter='r.MP9601', telescope_instrument='CFHT/MegaCam'):
    '''
    Parse through objects in ssois query and filter out images of desired filter, type, exposure time, and instrument
    '''
    
    assert camera_filter in ['r.MP9601', 'u.MP9301']
    ret_table = []
    good_table = 0
    bad_table = 0

    table_reader = ascii.get_reader(Reader=ascii.Basic)
    table_reader.inconsistent_handler = _skip_missing_data
    table_reader.header.splitter.delimiter = '\t'
    table_reader.data.splitter.delimiter = '\t'
    table = table_reader.read(ssois_return)

    for row in table:
        # check if a dbimages object exists
        mjd = row['MJD']

        # Excludes the OSSOS wallpaper.
        # note: 'Telescope_Insturment' is a typo in SSOIS's return format
        if (row['Telescope_Insturment'] == telescope_instrument) and (row['Filter'] == camera_filter) \
                and (row['Image'].endswith('{}'.format(imagetype))) and not row['Image_target'].startswith('WP'):
            if (int(row['Exptime']) == 287) or (int(row['Exptime']) == 387 ) or (int(row['Exptime']) == 500 ):
                good_table += 1
                ret_table.append(row)
            else:
                bad_table += 1
        if (row['Telescope_Insturment'] == telescope_instrument) and (row['Filter'] == camera_filter) \
            and not (int(row['Exptime']) > 286):
            bad_table += 1
            
    print "Searching for object %s " % object_name
    
    if (good_table > 0):
        print " %d usable images found" % good_table
    
    if (bad_table > 0):        
        print " %d unusable images found " % bad_table
    
    if (good_table == 0) and (bad_table == 0):
        print " no images found "

    return ret_table

def _skip_missing_data(str_vals, ncols):
    """
    add a extra column if one is missing, else return None.
    """
    if len(str_vals) == ncols - 1:
        str_vals.append('None')
        return str_vals
    else:
        #raise ValueError("not enough columns in table")
        print '  Not enough columns in data table'
    
def main():
    '''
    Input asteroid family, filter type, and image type to query SSOIS
    '''
    
    parser = argparse.ArgumentParser(description='Run SSOIS and return the available images in a particular filter.')

    parser.add_argument("--filter",
                    action="store",
                    default='r',
                    dest="filter",
                    choices=['r', 'u'],
                    help="passband: default is r'")
    parser.add_argument("--family", '-f',
                        action="store",
                        default="testfamily/testfamily_family.txt",
                        help='list of objects to query')
    parser.add_argument('--type',
                        default='p',
                        choices=['o', 'p', 's'], 
                        help="restrict type of image (unprocessed, reduced, calibrated)")

    args = parser.parse_args()
    
    get_image_info(args.family, args.filter, args.type)


if __name__ == '__main__':
    main()
