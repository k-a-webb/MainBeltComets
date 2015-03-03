import os
import io
import sep
import re
import vos
import tempfile
import urllib2 as url
import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.io import ascii
from astropy.time import Time
import argparse
from scipy.spatial import cKDTree
import math
#from pyraf import iraf
import pandas as pd

client = vos.Client()

import sys
sys.path.append('/Users/admin/Desktop/MainBeltComets/getImages/ossos_scripts/')

from ossos_scripts import storage
import ossos_scripts.wcs as wcs
from ossos_scripts.storage import get_astheader, exists
from get_images import get_image_info
from find_family import find_family_members
import get_stamps

''' 
Preforms photometry on .fits files given an input of family name and object name
Identifies object in image from predicted coordinates, magnitude (and eventually shape)
Assumes files organised as:
dir_path_base/familyname/familyname_objectname/*.fits       - images to do photometry on
  or dir_path_base/familyname/familyname_stamps/*.fits
dir_path_base/familyname/*_images.txt                       - list of image exposures, predicted RA and DEC, dates etc.
'''

'''
NOTES, to do:
    - add in psf.py to get PSF
    - reformat list appending
'''

def main(): 
    
    parser = argparse.ArgumentParser(
                        description='For a given set of fits images, \
                         preforms photometry on a specified object given an aperture size and threshold, \
                        and then selects the object in the image from the predicted coordinates, magnitude, and eventually shape')
    parser.add_argument("--family", '-f',
                        action="store",
                        default='all',
                        help="Asteroid family name. Usually the asteroid number of the largest member.")
    parser.add_argument("--aperture", '-ap',
                        action='store',
                        default=10.0,
                        help='aperture (degree) of circle for photometry.')
    parser.add_argument("--thresh", '-t',
                        action='store',
                        default=3.5,
                        help='threshold value for photometry (sigma above background).')
    parser.add_argument("--object", '-o',
                        action='store',
                        default=None,
                        help='The object to preform photometry on.')
    parser.add_argument("--filter",
                        action="store",
                        default='r',
                        dest="filter",
                        choices=['r', 'u'],
                        help="passband: default is r'")
    parser.add_argument('--type',
                        default='p',
                        choices=['o', 'p', 's'], 
                        help="restrict type of image (unprocessed, reduced, calibrated)")
    parser.add_argument('--elim', '-el',
                        default='0.6',
                        help="limit on degree of ellipticity")
                            
    args = parser.parse_args()
    
    find_(args.family, args.object, float(args.aperture), float(args.thresh), args.filter, args.type)
    
def find_objects_by_phot(familyname, objectname=None, ap=10.0, th=3.5, filtertype='r', imagetype='p'):

    out_filename = '{}_r{}_t{}_output.txt'.format(familyname, ap, th)
    with open('asteroid_families/{}/{}_stamps/{}'.format(familyname, familyname, out_filename), 'w') as outfile:
        outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format('Object', "Image", 'flux', 'mag', 'RA', 'DEC', 'ecc', 'index'))
    
    # initiate directories
    init_dirs(familyname, objectname)
    
    # From the given input, identify the desired filter and rename appropriately
    if filtertype.lower().__contains__('r'):
        filtertype = 'r.MP9601'  # this is the old (standard) r filter for MegaCam
    if filtertype.lower().__contains__('u'):
        filtertype = 'u.MP9301'
    
    if  os.path.exists('asteroid_families/{}/{}_images.txt'.format(familyname, familyname)):
        expnum_list = []
        image_list = []
        with open(image_list_path) as infile:
            filestr = infile.read()
            fileline = filestr.split('\n')
            for item in fileline[1:]:
                if len(item.split()) > 0:
                    image_list.append(item.split()[0])
                    expnum_list.append(item.split()[1])
    else:  
        image_list, expnum_list, ra_list, dec_list = get_image_info(familyname, filtertype, imagetype)
        
    if objectname == None:
        for index, imageobject in enumerate(image_list):
            print 'Finding asteroid {} in family {} '.format(objectname, familyname)
            iterate_thru_images(familyname, imageobject, expnum_list[index], ap, th, filtertype, imagetype, elim)    
    else:  
        for index, imageobject in enumerate(image_list):
            if objectname == imageobject:
                print 'Finding asteroid {} in family {} '.format(objectname, familyname)
                iterate_thru_images(familyname, objectname, expnum_list[index], ap, th, filtertype, imagetype, elim)
        

def iterate_thru_images(familyname, objectname, expnum_p, username, password, ap=10.0, th=5.0, filtertype='r', imagetype='p'):
            
    # initiate directories
    init_dirs(familyname, objectname) 
    
    try:
        print "-- Querying JPL Horizon's ephemeris"
        e_jpl = get_ell(familyname, objectname, expnum_p)
        mag_list_jpl, r_err = get_mag_radius(familyname, objectname)
    except Exception, e:
        print 'WARNING: No images exist, {}'.format(e)
        raise
      
    try:
        print "-- Preforming photometry on image {} ".format(expnum_p)
        table, exptime, zeropt, size, pvwcs, stamp_found = find_all_objects(familyname, objectname, expnum_p, username, password, ap, th, filtertype, imagetype)
    except Exception, e:
        print "ERROR: Error while doing photometry, {}".format(e)
        return
        
    if stamp_found == False:
        print "WARNING: no stamps exist"
        get_stamps.get_one_stamp(objectname, expnum_p, 0.02, username, password, familyname)
        return
    
    r_new, r_old, enough = check_num_stars(table, size, objectname, expnum_p, username, password, familyname)
    if enough == False:
        return
                    
    try:
        i_list = find_neighbours(table, objectname, expnum_p, pvwcs, r_err)
        if len(i_list) == 0:
            print '-- Expanding radius of nearest neighbour search by 1.5x'
            i_list = find_neighbours(table, objectname, expnum_p, pvwcs, 1.5*r_err)
            if len(i_list) == 0:
                print 'WARNING: No nearest neighbours were found within {} ++++++++++++++++++++++'.format(r_err*1.5)
                ascii.write(table, 'asteroid_families/temp_phot_files/{}_phot.txt'.format(expnum_p))
                return
                
        good_neighbours, mean = iden_good_neighbours(expnum_p, i_list, table, zeropt, mag_list_jpl, e_jpl, pvwcs)
        print good_neighbours

        print_output(familyname, objectname, expnum_p, good_neighbours, ap, th)
        
        '''if len(good_neighbours) == 1:
            print '-- Cutting out recentered postage stamp'
            recut_stamp(familyname, objectname, expnum_p, good_neighbours, r_old, username, password)
        '''
        
    except TypeError, e:
        return
    except ValueError:
        print 'WARNING: Only one object in image'
        get_stamps.get_one_stamp(objectname, expnum_p, r_new, username, password, familyname)
        return            
                                    
def find_all_objects(familyname, objectname, expnum_p, username, password, ap, th, filtertype, imagetype):    
    stamp_found = False      
    for file in client.listdir(vos_dir): # images named with convention: object_expnum_RA_DEC.fits
        if file.endswith('.fits') == True:
            objectname_file = file.split('_')[0]
            expnum_file = file.split('_')[1]
            if (expnum_file == expnum_p) and (objectname_file == objectname):
                stamp_found = True
                file_path = '{}/{}'.format(stamps_dir, file)
                storage.copy('{}/{}'.format(vos_dir, file), file_path)
                try:
                    with fits.open('{}/{}'.format(stamps_dir, file)) as hdulist: 
                        print hdulist.info()
                        
                        try:
                            if hdulist[0].data is None:
                                data1 = fits.getdata(file_path, 1)
                                data2 = fits.getdata(file_path, 2)
                                header = fits.getheader(file_path, 1)
                                header2 = fits.getheader(file_path, 2)
                                size = header['NAXIS1']+header2['NAXIS1']
                            
                                table1 = sep_phot(data1, ap, th)
                                table2 = sep_phot(data2, ap, th)
                                table = vstack([table1, table2])
                                #ascii.write(table, os.path.join(stamps_dir, '{}_phot.txt'.format(expnum_p))) 
                            
                            else:
                                data = fits.getdata(file_path)
                                header = fits.getheader(file_path)
                                size = header['NAXIS1']
                                table = sep_phot(data, ap, th)
                                #ascii.write(table, os.path.join(stamps_dir, '{}_phot.txt'.format(expnum_p)))
                                
                        except LookupError, e: # maybe not correct error type?
                            print "ERROR: {} xxxxxxxxxx".format(e)
                            #get_stamps.get_one_stamp(objectname, expnum_p, 0.025, username, password, familyname)
                            raise     
                        except ValueError, e:
                            print 'ERROR: {} **********'.format(e)
                            get_stamps.get_one_stamp(objectname, expnum_p, 0.02, username, password, familyname)
                            raise
                
                    os.unlink('{}/{}'.format(stamps_dir, file))    
                    
                    pvwcs = wcs.WCS(header)
                    zeropt = header['PHOTZP']
                    exptime = header['EXPTIME']
                    
                    return table, exptime, zeropt, size, pvwcs, stamp_found
                    
                except IOError, e:
                    print 'ERROR: {} []][][][][][][][]'.format(e)
                    #get_stamps.get_one_stamp(objectname, expnum_p, 0.02, username, password, familyname)
                    raise
        
def recut_stamp(familyname, objectname, expnum_p, object_data, r_old, username, password):
        
    if object_data is not None:
        #for i in range(0, len(object_data)):
            # objectname, expnum_p, r_new, RA, DEC, username, password, familyname
        print '-- Cutting recentered stamp'
        get_stamps.centered_stamp(objectname, expnum_p, r_old, object_data[0][5], object_data[0][6], username, password, familyname)
                            
def query_jpl(objectname, time_start, time_end, params, step=1, su='d'):
    
    # Construct the query url
    s = "'"
    for p in params:
        s += "{},".format(p)

    # form URL pieces that Horizon needs for its processing instructions
    urlArr = ["http://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&COMMAND=",
                  '',
                  "&MAKE_EPHEM='YES'&TABLE_TYPE='OBSERVER'&START_TIME=",
                  '',
                  "&STOP_TIME=",
                  '',
                  "&STEP_SIZE=",
                  '',
                  "&QUANTITIES=" + s,
                  "&CSV_FORMAT='YES'"]
          
    # change the object name, start and end times, and time step into proper url-formatting
    url_style_output = []
    for obj in [objectname, time_start, time_end]:
        os = obj.split()
        if len(os) > 1:
            ob = "'" + os[0] + '%20' + os[1] + "'"
        else:
            ob =  "'" + objectname + "'"
        url_style_output.append(ob)
    step = "'" + str(step) + "%20" + su + "'"
 
    # URL components
    urlArr[1] = url_style_output[0]  # formatted object name
    urlArr[3] = url_style_output[1]  # start time
    urlArr[5] = url_style_output[2]  # end time
    urlArr[7] = step  # timestep   
    urlStr = "".join(urlArr)  # create the url to pass to Horizons
   
    # Query Horizons; if it's busy, wait and try again in a minute
    done = 0
    while not done:
        urlHan = url.urlopen(urlStr)
        urlData = urlHan.readlines()
        urlHan.close()
        if len(urlData[0].split()) > 1:
            if "BUSY:" <> urlData[0].split()[1]:
                done = 1
            else:
                print urlData[0],
                print "Sleeping 60 s and trying again"
                time.sleep(60)
        else:
            done = 1
    
    EPHEM_CSV_START_MARKER = '$$SOE'
    EPHEM_CSV_END_MARKER = '$$EOE'
    ephemCSV_start = None
    ephemCSV_end = None
    for i, dataLine in enumerate(urlData):
        if dataLine.strip() == EPHEM_CSV_START_MARKER:
            ephemCSV_start = i
        elif dataLine.strip() == EPHEM_CSV_END_MARKER:
            ephemCSV_end = i
    assert ephemCSV_start is not None, 'No ephem start'
    assert ephemCSV_end is not None, 'No ephem end'
    assert ephemCSV_start < ephemCSV_end, 'Ephem are a bit odd'
    # The header, the lines *after* the start marker, up to (but not including, because it is a slice) the end marker
    csv_lines = [urlData[ephemCSV_start - 2]] + urlData[ephemCSV_start + 1: ephemCSV_end]
    ephemCSV = ''.join(csv_lines)
    ephemCSVfile = io.BytesIO(ephemCSV)
    ephemerides = pd.DataFrame.from_csv(ephemCSVfile)
    
    return ephemerides
    
def get_mag_radius(familyname, objectname):
    
    if type(objectname) is not str:
        objectname = str(objectname)
    
    # from familyname_images.txt get date range of images for objectname
    date_range = []
    with open('asteroid_families/{}/{}_images.txt'.format(familyname, familyname)) as infile:
        for line in infile.readlines()[1:]:
            if len(line.split()) > 0:
                if objectname == line.split()[0]:
                    date_range.append(float(line.split()[5]))   
                    
    date_range_t = Time(date_range, format='mjd', scale='utc')
    assert  len(date_range_t.iso) > 0
    time_start = ((date_range_t.iso[0]).split())[0] + ' 00:00:00.0'
    time_end = ((date_range_t.iso[-1]).split())[0] + ' 00:00:00.0'

    if time_start == time_end:
        time_end = add_day(time_end, date_range_t)
        
    print " Date range in query: {} -- {}".format(time_start, time_end)

    # change date format from 01-01-2001 00:00 to 01-Jan-2001 00:00
    date_start = change_date(time_start)
    date_end = change_date(time_end)
    
    ephemerides = query_jpl(objectname, date_start, date_end, params=[9, 36])
    
    mag_list =  np.array(ephemerides.icol(2))
    ra_avg = np.mean(np.mean(ephemerides.icol(3)))
    dec_avg = np.mean(np.mean(ephemerides.icol(4)))
        
    if ra_avg > dec_avg:
        r_err = ra_avg / 0.184
    else:
        r_err = dec_avg / 0.184

    if r_err < 25: # 0.003 deg * 3600 "/deg / 0.187 "/pix
        r_err = 25
    
    return mag_list, r_err

def get_ell(familyname, objectname, expnum):
    '''
    Queries the JPL Horizon's ephemeris for rate of change of RA and DEC for a specific day
    '''
    
    if type(objectname) is not str:
        objectname = str(objectname)
    
    # from familyname_images.txt get date range of images for objectname
    image_table = pd.read_table(image_list_path, usecols=[0, 1, 5], header=0, names=['Object', 'Image', 'time'], sep=' ', dtype={'Object':object})
    index = image_table.query(('Image == "{}"'.format(expnum)) and ('Object == "{}"'.format(objectname)))

    date = index['time']          
    date_t = Time(date, format='mjd', scale='utc')
    assert  len(date_t.iso) > 0

    time_start = (date_t.iso)[0].split()[0]+' 00:00:00.0'
    
    time_end_date = (((date_t.iso[0]).split())[0]).split('-')
    day_add_one = int(time_end_date[2])+1
    if day_add_one < 10:
        day = '0{}'.format(day_add_one)
    else:
        day = day_add_one
    if day > 27:
        day = 1
        month = int(time_end_date[1]) + 1
    else:
        month = time_end_date[1]
    time_end = '{}-{}-{} 00:00:00.0'.format(time_end_date[0], month, day)
        
    # change date format from 01-01-2001 00:00 to 01-Jan-2001 00:00
    date_start = change_date(time_start)
    date_end = change_date(time_end)
    
    ephemerides = query_jpl(objectname, date_start, date_end, params=[3])
    
    a = ephemerides['dRA*cosD'][0]
    b = ephemerides['d(DEC)/dt'][0]
    e = (1-(b/a)**2)**0.5
                        
    return e

def change_date(date):
    '''
    Convert time format 01-01-2001 00:00:00.00 to 01-Jan-2001 00:00
    '''
    date_split = date.split('-')
    date_strip = (date_split[2]).split()
    month = int(date_split[1])
    month_name = ['nan', 'Jan', "Feb", 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    for i in range(0,13):
        if i == month:
            month_jpl = month_name[i]
    date_new = date_split[0]+'-'+month_jpl+'-'+date_strip[0]+' 00:00'
    
    return date_new    
    
def sep_phot(data, ap, th):
    ''' 
    Preforms photometry by SEP, similar to source extractor 
    '''

    # Measure a spatially variable background of some image data (numpy array)
    try:    
        bkg = sep.Background(data) #, mask=mask, bw=64, bh=64, fw=3, fh=3) # optional parameters
    except:
        data = data.byteswap(True).newbyteorder()
        bkg = sep.Background(data) #, mask=mask, bw=64, bh=64, fw=3, fh=3) # optional parameters
        
    # Evaluate the spatially variable background and RMS:
    back = bkg.back() # creates an array, same shape and type as data
    rms = bkg.rms()   # creates an array, same shape and type as data
        
    # Directly subtract the background from the data in place
    bkg.subfrom(data)
    bkg.globalback    # Global "average" background level
    bkg.globalrms     # Global "average" RMS of background
        
    # for the background subtracted data, detect objects in data given some threshold
    thresh = th * bkg.globalrms    # ensure the threshold is high enough wrt background        
    objs = sep.extract(data, thresh)

    # calculate the Kron radius for each object, then we perform elliptical aperture photometry within that radius
    kronrad, krflag = sep.kron_radius(data, objs['x'], objs['y'], objs['a'], objs['b'], objs['theta'], ap)
    flux, fluxerr, flag = sep.sum_ellipse(data, objs['x'], objs['y'], objs['a'], objs['b'], objs['theta'], 2.5*kronrad, subpix=1)
    flag |= krflag  # combine flags into 'flag'
   
    r_min = 1.75  # minimum diameter = 3.5
    use_circle = kronrad * np.sqrt(objs['a'] * objs['b']) < r_min
    x = objs['x']
    y = objs['y']
    cflux, cfluxerr, cflag = sep.sum_circle(data, x[use_circle], y[use_circle],
                                            r_min, subpix=1)
    flux[use_circle] = cflux
    fluxerr[use_circle] = cfluxerr
    flag[use_circle] = cflag
   
    # write to ascii table
    table = Table([objs['x'], objs['y'], flux, objs['a'], objs['b']], names=('x', 'y', 'flux', 'a', 'b'))
    return table

def find_neighbours(septable, objectname, expnum, pvwcs, r_err):
    '''
    Computes the nearest neighbours to predicted coordinates within an RA/DEC uncertainty circle
    '''
    tree = cKDTree(zip((np.array(septable['x'])).ravel(), (np.array(septable['y'])).ravel()))    
    
    with open('{}/{}'.format(family_dir, imageinfo)) as infile:
        for line in infile.readlines()[1:]:
            assert len(line.split()) > 0
            expnum_fromfile = line.split()[1]
            object_fromfile = line.split()[0]
            
            # for entries in *_images.txt that correspond to images of the object
            if (expnum_fromfile == expnum) & (object_fromfile == objectname):
                pRA = float(line.split()[3])
                pDEC = float(line.split()[4])
                
                pRA_pix, pDEC_pix = pvwcs.sky2xy(pRA, pDEC) # convert from WCS to pixels
                print " Predicted RA and DEC: {}  {}".format(pRA, pDEC)
                print "  in pixels: {} {}".format(pRA_pix, pDEC_pix)
                
                # parse through table and get RA and DEC closest to predicted coordinates (in pixels)
                coords = np.array([pRA_pix, pDEC_pix])
                i_list = tree.query_ball_point(coords, r_err)
            
                return i_list

    '''table = pd.read_table(image_list_path, usecols=[0, 1, 3, 4], header=0, names=['Object', 'Image', 'RA', 'DEC'], sep=' ', dtype={'Object':object})
    index = table.query(('Object == "{}"'.format(objectname)))
    index2 = index.query(('Image == "{}"'.format(expnum)))
    print index2
    pRA = index2.get_value(0, ['RA'])
    pDEC = index2.get_value(0, ['DEC'])
    print pRA, pDEC, type(pRA), type(pDEC)
    pRA_pix, pDEC_pix = pvwcs.sky2xy(pRA, pDEC) # convert from WCS to pixels
    #print " Predicted RA and DEC: {}  {}".format(pRA, pDEC)
    #print "  in pixels: {} {}".format(pRA_pix, pDEC_pix)
    
    coords = np.array([pRA_pix, pDEC_pix])
    i_list = tree.query_ball_point(coords, r_err)
    print i_list

    return i_list'''
    
    
def iden_good_neighbours(expnum, i_list, septable, zeropt, mag_list_jpl, e_jpl, pvwcs):
    '''
    Selects nearest neighbour object from predicted coordinates as object of interest
    In order:
        Compares measured apparent magnitude to predicted, passes if in range of values
        Calculates eccentricity, passes if greater than minimum value that is inputted
    '''
        
    mRA_pix = None
    mag_sep_list = []
    all_mag = []
    index_list = []
    flux_list = []
    mRA_pix_list = []
    mDEC_pix_list = []
    ecc_list = []
    x_list = []
    y_list = []
    temp = []
    
    mean = np.mean(mag_list_jpl)
    maxmag = np.amax(mag_list_jpl)
    minmag = np.amin(mag_list_jpl)
    
    if ( 2 > maxmag - minmag):
        magrange = 2
    else:
        magrange = maxmag - minmag
    
    for i in i_list:
        flux = septable[i][2]
        x = septable[i][0]
        y = septable[i][1]
        a = septable[i][3]
        b = septable[i][4]
        ecc = (1-(b/a)**2)**0.5

        if flux > 0:
            mag_sep = -2.5*math.log10(flux)+zeropt
            all_mag.append(mag_sep)
            temp.append(x)
        
            # apply both eccentricity and magnitude conditions
            if (abs(mag_sep - mean) < magrange) & (float(e_jpl) - 0.1 < ecc < float(e_jpl) + 0.1):
                mRA_pix = septable[i][0]
                mDEC_pix = septable[i][1]
                
                index_list.append(i)
                flux_list.append(flux)
                mRA_pix_list.append(mRA_pix)
                mDEC_pix_list.append(mDEC_pix)
                mag_sep_list.append(mag_sep)
                ecc_list.append(ecc)
                x_list.append(x)
                y_list.append(y)   
                        
            # if not both, try each
            elif (abs(mag_sep - mean) < magrange):
                mRA_pix = septable[i][0]
                mDEC_pix = septable[i][1]

                index_list.append(i)
                flux_list.append(flux)
                mRA_pix_list.append(mRA_pix)
                mDEC_pix_list.append(mDEC_pix)
                mag_sep_list.append(mag_sep)
                ecc_list.append(ecc)
                x_list.append(x)
                y_list.append(y)
                                
                print "WARNING: Eccentricity is outside range: calculated, {}, measured, {}".format(e_jpl, ecc)       
                    
    good_neighbours = Table([index_list, flux_list, mag_sep_list, x_list, y_list, mRA_pix_list, mDEC_pix_list, ecc_list], 
                names=('index', 'flux', 'mag', 'x', 'y', 'RA', 'DEC', 'ecc'))
                    
    if len(all_mag) == 0:
        print "WARNING: Flux of nearest neighbours measured to be 0.0 00000000000000000000"
        return
    
    if mRA_pix == None:
        print "WARNING: Magnitude condition could not be satisfied 000000000000000000000"
        print '  Nearest neighbour list: {}'.format(i_list)
        print "  Mag mean, accepted error, and fitting mags: {} {} {} {}".format(mean, magrange, all_mag, temp)
        ascii.write(septable, 'asteroid_families/temp_phot_files/{}_phot.txt'.format(expnum))
        return
        
    else: 
        #print "   Flux, mag: {}, {}".format(flux, mag_sep)     
        mRA, mDEC = pvwcs.xy2sky(mRA_pix, mDEC_pix) # convert from pixels to WCS
        #print " Measured RA and DEC: {}  {}".format(mRA, mDEC)
        #print "  Coordinates: {} {}".format(mRA_pix, mDEC_pix)
        #print " Difference: {} {}".format(mRA - pRA, mDEC - pDEC)
        
        return good_neighbours, mean 
        
def check_num_stars(table, size, objectname, expnum_p, username, password, familyname):
        
    # r_old = size * 0.184 / 3600
    # r_new = r_old + 0.05
    
    enough = True
    
    if size < 200:
        r_new = 0.02
        r_old = 0.02
    if 200 < size < 500:
        r_new = 0.03
        r_old = 0.02
    if size > 500:
        r_new = 0.04
        r_old = 0.02
    
    print '  Number of objects in image: {}'.format(len(table))
    if 0 < len(table) < 30:
        enough = False
        print 'Not enough stars () in the stamp ////////////////////////////////'.format(len(table))                    
        #if size < 500:   
            #get_stamps.get_one_stamp(objectname, expnum_p, r_new, username, password, familyname)
        if size > 500:
            print 'SIZE IS AT MAX ({}) and not enough stars <<<<<<<<<<<<<<<<<<<'.format(size)
            #get_stamps.get_one_stamp(objectname, expnum_p, 0.05, username, password, familyname)
    
    return r_new, r_old, enough 
    
def print_output(familyname, objectname, expnum_p, object_data, ap, th):

    if object_data is None:
        print "WARNING: Could not identify object {} in image".format(objectname, expnum_p)
        
    else:
        out_filename = '{}_r{}_t{}_output.txt'.format(familyname, ap, th)
        with open('asteroid_families/{}/{}_stamps/{}'.format(familyname, familyname, out_filename), 'a') as outfile:
            try:
                for i in range(0, len(object_data)):
                    #'Object', "Image", 'flux', 'mag', 'RA', 'DEC', 'ecc', 'index'
                    outfile.write('{} {} {} {} {} {} {} {}\n'.format(
                          objectname, expnum_p, object_data[i][1], object_data[i][2], object_data[i][5], object_data[i][6], 
                          object_data[i][7], object_data[i][0]))
            except:
                print "ERROR: cannot write to outfile <<<<<<<<<<<<<<<<<<<<<<<<<<"
    
        return object_data
                    
def init_dirs(familyname, objectname):
    
    # initiate vos directories 
    global vos_dir
    vos_dir = 'vos:kawebb/postage_stamps/{}'.format(familyname)
    assert exists(vos_dir, force=True)

    
    # initiate local directories
    dir_path_base = '/Users/admin/Desktop/MainBeltComets/getImages/asteroid_families'
    global family_dir
    family_dir = os.path.join(dir_path_base, familyname)
    if os.path.isdir(family_dir) == False:
        print "Invalid family name"
    global stamps_dir
    stamps_dir = os.path.join(family_dir, familyname+'_stamps')
    if os.path.isdir(stamps_dir) == False:
        os.makedirs(stamps_dir)
        
    global imageinfo
    imageinfo = '{}_images.txt'.format(familyname)
    global image_list_path
    image_list_path = '{}/{}'.format(family_dir, imageinfo)
        
    assert os.path.exists('{}/{}'.format(family_dir, imageinfo))
    
    return family_dir, stamps_dir, vos_dir


def add_day(time_end, date_range_t):
    print "WARNING: only searching for one day"
    time_end_date = (((date_range_t.iso[-1]).split())[0]).split('-')
    day_add_one = int(time_end_date[2])+1
    if day_add_one < 10:
        day = '0{}'.format(day_add_one)
    else:
        day = day_add_one
    month = int(time_end_date[1])
    if day > 27:
        day = 1
        if month == 12:
            month = 1
        else:
            month = month + 1

    time_end = '{}-{}-{} 00:00:00.0'.format(time_end_date[0], month, day)        
    return time_end

    
if __name__ == '__main__':
    main()
    