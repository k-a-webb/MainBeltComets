import datetime
import os
from astropy.io import ascii
from astropy.time import Time
import requests
import argparse
from ossos_scripts.ssos import Query
import time
import pandas as pd

dir_path_base = '/Users/admin/Desktop/MainBeltComets/getImages/asteroid_families'

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
    parser.add_argument('--suffix',
                        default=None,
                        help="Suffix of object list file")
                        
    args = parser.parse_args()
            
    get_image_info(args.family, args.suffix, args.filter, args.type)

def get_image_info(familyname, suffix=None, filtertype='r', imagetype='p'):
    '''
    Query the ssois ephemeris for images of objects in a given family. Then parse through for desired image type, 
    filter, exposure time, and telescope instrument
    '''
    
    family_dir = os.path.join(dir_path_base, familyname)
    if os.path.isdir(family_dir) == False:
        print "Invalid family name or directory does not exist"
    
    if suffix == None:
        family_list = '{}/{}_family.txt'.format(family_dir, familyname)
        output = '{}_images.txt'.format(familyname)
    else:
        family_list = '{}/none_family_{}.txt'.format(family_dir, suffix)
        output = '{}_images_{}'.format(familyname, suffix)
    
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

    # Query a specific ephemeris depending on object naming convention. CADC : search = bynameCADC , MPC : search = bynameMPC
    # Specify search parameters
        # ephemeris, name, date range, resolve image extension, resolve to x,y, positional uncertainty?
        # http:// www3.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/cadcbin/ssos/ssosclf.pl?lang=en; object=elst-pizarro; search=bynameMPC; epoch1=2013+01+01; epoch2=2015+1+16; eellipse=; eunits=arcseconds; extres=yes; xyres=yes; format=tsv
        
    print "-------------------- \n Searching for images of objects in family {} from CFHT/Megacam from the MPC ephemeris".format(familyname)
    print " with filter {} and exposure time of 287, 387, 500 seconds (OSSOS data) \n--------------------".format(filtertype)
        
    image_list = []
    for object_name in object_list[:50]:#len(object_list)-1]:
        query = Query(object_name, search_start_date=search_start_date, search_end_date=search_end_date)

        # GET/REVIEW THE DATA RETURNED FROM THE SEARCH
        # Parse the data for Image, Instrument, Filter and create table for each object
        # Download files with appropriate names? ie object+SSOISfilename ??
        
        try:
            objects = parse_ssois_return(query.get(), object_name, imagetype, camera_filter=filtertype)
        except IOError:
            time.sleep(20)
            print "Sleeping 20 seconds"
            objects = parse_ssois_return(query.get(), object_name, imagetype, camera_filter=filtertype)
        x=0
        index=[]
        if len(objects) > 0:
            for line in objects:
                image_list.append(object_name)
                index.append(x)
                x+=1
                table_arrays = {'object': object_name, 'image': line['Image'], 'Exptime': line['Exptime'], 'RA': line['Object_RA'], \
                                        'DEC': line['Object_Dec'], 'date': Time(line['MJD'], format='mjd', scale='utc'), 'filter': line['Filter']}
                                        
    table = pd.DataFrame(data=table_arrays, index=index)
    table.to_csv('{}/{}'.format(family_dir, output), sep='\t', encoding='utf-8', index=False) 
               
    return image_list
                    
def parse_ssois_return(ssois_return, object_name, imagetype, camera_filter='r.MP9601', telescope_instrument='CFHT/MegaCam'):
    '''
    Parse through objects in ssois query and filter out images of desired filter, type, exposure time, and instrument
    '''
    
    assert camera_filter in ['r.MP9601', 'u.MP9301']
    
    ret_table = []
    good_table = 0

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
                and (row['Image'].endswith('{}'.format(imagetype))) and not str(row['Image_target']).startswith('WP'):
                
            if (int(row['Exptime']) == 287) or (int(row['Exptime']) == 387 ) or (int(row['Exptime']) == 500 ):
                good_table += 1
                ret_table.append(row)
            
    print "Searching for object %s " % object_name
    
    if (good_table > 0):
        print " %d images found" % good_table

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
    

if __name__ == '__main__':
    main()
