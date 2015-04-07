import datetime
import os
from astropy.io import ascii
from astropy.time import Time
import requests
import argparse

import sys

sys.path.append('User/admin/Desktop/OSSOS/MOP/src/ossos-pipeline/ossos')
# from ossos.ssos import Query WHY DOESNT THIS WORK?
from ossos_scripts.ssos import Query

import time
import pandas as pd
import find_family

_DIR_PATH_BASE = os.path.dirname(os.path.abspath(__file__))
_FAMILY_LISTS = '{}/family_lists'.format(_DIR_PATH_BASE)
_OUTPUT_DIR = '{}/image_lists'.format(_DIR_PATH_BASE)


def main():
    """
    Input asteroid family, filter type, and image type to query SSOIS
    """

    parser = argparse.ArgumentParser(description='Run SSOIS and return the available images in a particular filter.')

    parser.add_argument("--filter",
                        action="store",
                        default='r',
                        dest="filter",
                        choices=['r', 'u'],
                        help="passband: default is r'")
    parser.add_argument("--family", '-f',
                        action="store",
                        default="all",
                        help='list of objects to query')
    parser.add_argument('--type',
                        default='p',
                        choices=['o', 'p', 's'],
                        help="restrict type of image (unprocessed, reduced, calibrated)")

    args = parser.parse_args()

    get_image_info(str(args.family), args.filter, args.type)


def get_image_info(familyname, filtertype='r', imagetype='p'):
    """
    Query the ssois ephemeris for images of objects in a given family. Then parse through for desired image type, 
    filter, exposure time, and telescope instrument

    CADC : search = bynameCADC , MPC : search = bynameMPC
    ephemeris, name, date range, resolve image extension, resolve to x,y, positional uncertainty?
    http:// www3.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/cadcbin/ssos/ssosclf.pl?lang=en; object=elst-pizarro; 
        search=bynameMPC; epoch1=2013+01+01; epoch2=2015+1+16; eellipse=; eunits=arcseconds; extres=yes; xyres=yes; format=tsv
    """

    # establish input/output
    # family_list = '{}/{}_family.txt'.format(_FAMILY_LISTS, familyname)
    family_list = 'image_lists/all_images_new_missing.txt'
    output = '{}/{}_images_new_found.txt'.format(_OUTPUT_DIR, familyname)

    if os.path.exists(family_list):
        with open(family_list) as infile:
            filestr = infile.read()
        object_list = filestr.split('\n')  # array of objects to query
    else:
        if familyname == 'all':
            object_list = find_family.get_all_families_list()
        else:
            object_list = find_family.find_family_members(familyname)

    # From the given input, identify the desired filter and rename appropriately
    if filtertype.lower().__contains__('r'):
        filtertype = 'r.MP9601'  # this is the old (standard) r filter for MegaCam
    if filtertype.lower().__contains__('u'):
        filtertype = 'u.MP9301'

    # Define time period of image search, basically while MegaCam in operation
    search_start_date = Time('2013-01-01', scale='utc')  # epoch1=2013+01+01
    search_end_date = Time('2017-01-01', scale='utc')  # epoch2=2017+1+1

    # Setup output, label columns
    with open(output, 'w') as outfile:
        outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            "Object", "Image", "Exp_time", "RA", "DEC", "time", "filter"))

    print "----- \n Searching for images of objects in family {} from CFHT/Megacam from the MPC ephemeris".format(
        familyname)
    print " with filter {} and exposure time of 287, 387, 500 seconds (OSSOS data) \n-----".format(filtertype)

    image_list = []
    expnum_list = []
    ra_list = []
    dec_list = []
    for object_name in object_list[0:len(object_list) - 1]: # skip header lines
        query = Query(object_name, search_start_date=search_start_date, search_end_date=search_end_date)

        try:
            objects = parse_ssois_return(query.get(), object_name, imagetype, camera_filter=filtertype)
        except IOError:
            time.sleep(20)
            print "Sleeping 20 seconds"
            try:
                objects = parse_ssois_return(query.get(), object_name, imagetype, camera_filter=filtertype)
            except IOError:
                time.sleep(40)
                print "Sleeping 40 seconds"
                objects = parse_ssois_return(query.get(), object_name, imagetype, camera_filter=filtertype)

        for line in objects:
            with open(output, 'a') as outfile:
                image_list.append(object_name)
                expnum_list.append(line['Exptime'])
                ra_list.append(line['Object_RA'])
                dec_list.append(line['Object_Dec'])
                try:
                    outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(object_name,
                                                                        line['Image'], line['Exptime'],
                                                                        line['Object_RA'], line['Object_Dec'],
                                                                        Time(line['MJD'], format='mjd', scale='utc'),
                                                                        line['Filter']))
                except Exception, e:
                    print "Error writing to outfile: {}".format(e)

    return image_list, expnum_list, ra_list, dec_list


def parse_ssois_return(ssois_return, object_name, imagetype, camera_filter='r.MP9601',
                       telescope_instrument='CFHT/MegaCam'):
    """
    Parse through objects in ssois query and filter out images of desired filter, type, exposure time, and instrument
    """

    assert camera_filter in ['r.MP9601', 'u.MP9301']

    ret_table = []
    good_table = 0

    table_reader = ascii.get_reader(Reader=ascii.Basic)
    table_reader.inconsistent_handler = _skip_missing_data
    table_reader.header.splitter.delimiter = '\t'
    table_reader.data.splitter.delimiter = '\t'
    table = table_reader.read(ssois_return)

    for row in table:
        # Excludes the OSSOS wallpaper.
        # note: 'Telescope_Insturment' is a typo in SSOIS's return format
        exp_times = [287, 387, 500]
        if (row['Telescope_Insturment'] == telescope_instrument) and \
                (row['Filter'] == camera_filter) and \
                        int(row['Exptime']) in exp_times and \
                (row['Image'].endswith('{}'.format(imagetype))) and not \
                str(row['Image_target']).startswith('WP'):
            good_table += 1
            ret_table.append(row)

    print "Searching for object %s " % object_name

    if good_table > 0:
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
        # raise ValueError("not enough columns in table")
        print '  Not enough columns in data table'


if __name__ == '__main__':
    main()
