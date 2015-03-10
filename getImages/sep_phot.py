import os
import io
import sep
import re
import vos
import urllib2 as url
import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.io import ascii
from astropy.time import Time
import argparse
from scipy.spatial import cKDTree
import math
import pandas as pd
import sys
from shapely.geometry import Polygon, Point

from ossos_scripts import storage
import ossos_scripts.wcs as wcs
from ossos_scripts.storage import get_astheader, exists
#from ossos_scripts.util import match_lists

from get_images import get_image_info
from find_family import find_family_members
import get_stamps

client = vos.Client()

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
                            
    args = parser.parse_args()
    
    find_objects_by_phot(args.family, args.object, float(args.aperture), float(args.thresh), args.filter, args.type)
    
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
            iterate_thru_images(familyname, imageobject, expnum_list[index], ap, th, filtertype, imagetype)    
    else:  
        for index, imageobject in enumerate(image_list):
            if objectname == imageobject:
                print 'Finding asteroid {} in family {} '.format(objectname, familyname)
                iterate_thru_images(familyname, objectname, expnum_list[index], ap, th, filtertype, imagetype)
        

def iterate_thru_images(familyname, objectname, expnum_p, username, password, ap=10.0, th=5.0, filtertype='r', imagetype='p'):

    success = False            
    # initiate directories
    init_dirs(familyname, objectname) 
    
    try:
        print "-- Performing photometry on image {} ".format(expnum_p)
        septable, exptime, zeropt, size, pvwcs, stamp_found, start, end = get_fits_data(familyname, objectname, expnum_p, username, password, ap, th, filtertype, imagetype)
        if stamp_found == False:
            print "WARNING: no stamps exist"
            get_stamps.get_one_stamp(objectname, expnum_p, 0.02, username, password, familyname)
            return
    except Exception, e:
        print "ERROR: Error while doing photometry, {}".format(e)
        return
    
    try:
        print "-- Querying JPL Horizon's ephemeris"
        mag_list_jpl, r_sig = get_mag_rad(familyname, objectname)
        pRA, pDEC, ra_dot, dec_dot = get_coords(familyname, objectname, expnum_p, start, end)
    except Exception, e:
        print 'ERROR: Error while doing JPL query, {}'.format(e)
        raise
    
    table = append_table(septable, pvwcs, zeropt)
    transients, num_cat_objs = compare_to_catalogue(table, pvwcs)
    
    r_new, r_old, enough = check_num_stars(num_cat_objs, size, objectname, expnum_p, username, password, familyname)
    if enough == False:
        return
        
    print '-- Identifying object from nearest neighbours'
    i_list, found = find_neighbours(transients, pvwcs, r_sig, pRA, pDEC, expnum_p)
    if found == False:
        success = True
        return success
    
    '''try:
        
        good_neighbours, mean = iden_good_neighbours(expnum_p, i_list, transients, zeropt, mag_list_jpl, ra_dot, dec_dot, exptime, pvwcs)
        print good_neighbours

        print_output(familyname, objectname, expnum_p, good_neighbours, ap, th)
        
        if len(good_neighbours) == 1:
            involved = check_involvement(good_neighbours, table)
            if involved == False:
                print '-- Cutting out recentered postage stamp'
                #cut_centered_stamp(familyname, objectname, expnum_p, good_neighbours, r_old, username, password)
        
        success = True
        
    except Exception, e:
        print 'ERROR: {}'.format(e)
        #get_stamps.get_one_stamp(objectname, expnum_p, r_new, username, password, familyname)   
    '''
    good_neighbours, r_err = iden_good_neighbours(expnum_p, i_list, transients, zeropt, mag_list_jpl, ra_dot, dec_dot, exptime, pvwcs)
    print good_neighbours

    print_output(familyname, objectname, expnum_p, good_neighbours, ap, th)
    
    if len(good_neighbours) == 1:
        involved = check_involvement(good_neighbours, table, r_err)
        if involved == False:
            print '-- Cutting out recentered postage stamp'
            #cut_centered_stamp(familyname, objectname, expnum_p, good_neighbours, r_old, username, password)
    
    success = True
                
    return success
            
def get_fits_data(familyname, objectname, expnum_p, username, password, ap, th, filtertype, imagetype):    
    
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
                        #print hdulist.info()
                        
                        if hdulist[0].data is None:
                            print 'IMAGE is mosaic'
                            data1 = fits.getdata(file_path, 1)
                            data2 = fits.getdata(file_path, 2)
                            header = fits.getheader(file_path, 1)
                            header2 = fits.getheader(file_path, 2)
                            size = header['NAXIS1']+header2['NAXIS1']
                            print 'good until here'
                        
                            table1 = sep_phot(data1, ap, th)
                            print 'this time it worked'
                            table2 = sep_phot(data2, ap, th)
                            table = vstack([table1, table2])
                            #ascii.write(table, os.path.join(stamps_dir, '{}_phot.txt'.format(expnum_p)))
                            
                        else:
                            data = fits.getdata(file_path)
                            header = fits.getheader(file_path)
                            size = header['NAXIS1']
                            table = sep_phot(data, ap, th)
                            #ascii.write(table, os.path.join(stamps_dir, '{}_phot.txt'.format(expnum_p)))

                except Exception, e:
                    print 'ERROR: {} xxxxxxxxxxx'.format(e)
                    get_stamps.get_one_stamp(objectname, expnum_p, 0.03, username, password, familyname)
                    raise
               
                os.unlink('{}/{}'.format(stamps_dir, file))
                pvwcs = wcs.WCS(header)
                zeropt = header['PHOTZP']
                exptime = header['EXPTIME']
                start = '{} {}'.format(header['DATE-OBS'], header['UTIME'])
                end = '{} {}'.format(header['DATEEND'], header['UTCEND'])
               
                return table, exptime, zeropt, size, pvwcs, stamp_found, start, end
                                        
def get_mag_rad(familyname, objectname):
    
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
        
    #print " Date range in query: {} -- {}".format(time_start, time_end)
    
    ephemerides = query_jpl(objectname, time_start, time_end, params=[9, 36], step=1)
    
    mag_list =  np.array(ephemerides.icol(2))
    ra_sig = np.mean(np.mean(ephemerides.icol(3)))
    dec_sig = np.mean(np.mean(ephemerides.icol(4)))
        
    if ra_sig > dec_sig:
        r_sig = ra_sig / 0.184
    else:
        r_sig = dec_sig / 0.184

    if r_sig < 10: # 0.003 deg * 3600 "/deg / 0.187 "/pix
        r_sig = 10
    
    return mag_list, r_sig

def get_coords(familyname, objectname, expnum, time_start, time_end):
    '''
    Queries the JPL Horizon's ephemeris for rate of change of RA and DEC for a specific day
    '''
    
    if type(objectname) is not str:
        objectname = str(objectname)
    
    ephemerides = query_jpl(objectname, time_start, time_end, params=[1, 3], step='1', su='m')
            
    mid = int(len(ephemerides)/2)
    ra = ephemerides['R.A._(ICRF/J2000.0)'][mid]
    dec = ephemerides[' DEC_(ICRF/J2000.0)'][mid]
    ra_dot = ephemerides[' dRA*cosD'][mid]
    dec_dot = ephemerides['d(DEC)/dt'][mid]
    
    if dec.split()[0] < 0:
        sign = -1
    else:
        sign = 1
    
    ra_deg = float(HMS2deg(ra=ra))
    dec_deg = float(HMS2deg(dec=dec))
                        
    return ra_deg, dec_deg, ra_dot, dec_dot

def query_jpl(objectname, time_start, time_end, params, step, su='d'):
    
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
    
def sep_phot(data, ap, th):
    ''' 
    Preforms photometry by SEP, similar to source extractor 
    '''

    # Measure a spatially variable background of some image data (np array)
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
    table = Table([objs['x'], objs['y'], flux, objs['a'], objs['b'], objs['theta']], names=('x', 'y', 'flux', 'a', 'b', 'theta'))
    return table
    
def append_table(table, pvwcs, zeropt):
    
    ra_list = []
    dec_list = []
    mag_sep_list = []
    for row in range(len(table)):
        try:
            ra, dec = pvwcs.xy2sky(table['x'][row], table['y'][row])
            ra_list.append(ra)
            dec_list.append(dec)
            mag_sep_list.append(-2.5*math.log10(table['flux'][row])+zeropt)
        except:
            print row
    table['ra'] = ra_list
    table['dec'] = dec_list
    table['mag'] = mag_sep_list
    return table

def compare_to_catalogue(table, pvwcs):
    
    #convert sep table elements from pixels to WCS
    septable = pd.DataFrame(np.array(table))
    
    Eall = pd.read_table('catalogue/Eall.photcat', usecols=[0, 1, 4], header=0, names=['ra', 'dec', 'mag'], sep='      |     |    |   |  ', engine='python')
    Hall = pd.read_table('catalogue/Hall.photcat', usecols=[0, 1, 4], header=0, names=['ra', 'dec', 'mag'], sep='      |     |    |   |  ', engine='python')
    Lall = pd.read_table('catalogue/Lall.photcat', usecols=[0, 1, 4], header=0, names=['ra', 'dec', 'mag'], sep='      |     |    |   |  ', engine='python')
    Oall = pd.read_table('catalogue/Oall.photcat', usecols=[0, 1, 4], header=0, names=['ra', 'dec', 'mag'], sep='      |     |    |   |  ', engine='python')
    catalogue = pd.concat([Eall, Hall, Lall, Oall])
    catalogue.reset_index(drop=True, inplace=True)
    
    #pos1 = septable.as_matrix(columns=['ra', 'dec'])
    #pos2 = catalogue.as_matrix(columns=['ra', 'dec'])
    #match1, match2 = match_lists(pos1, pos2, tolerance=0.005111111)
    #print match1, match2
     
    sep_tol = 5 * 0.184 / 3600 # pixels to degrees
    mag_tol = 0.2
    
    cat_objs = 0
    trans_ra = []
    trans_dec = []
    trans_x = []
    trans_y = []
    trans_mag = []
    trans_a = []
    trans_b = []
    trans_theta = []
    for row in range(len(septable)):
        #index = catalogue[( abs(catalogue.ra - septable['ra'][row]) < 0.0051111 )]
        index = catalogue.query('({} < ra < {}) & ({} < dec < {}) & ({} < mag < {})'.format(septable['ra'][row]-sep_tol, septable['ra'][row]+sep_tol, 
                                                                                            septable['dec'][row]-sep_tol, septable['dec'][row]+sep_tol,
                                                                                            septable['mag'][row]-mag_tol, septable['mag'][row]+mag_tol))
        if len(index) == 0:
            trans_ra.append(septable['ra'][row])
            trans_dec.append(septable['dec'][row])
            trans_x.append(septable['x'][row])
            trans_y.append(septable['y'][row])
            trans_a.append(septable['a'][row])
            trans_b.append(septable['b'][row])
            trans_theta.append(septable['theta'][row])
            trans_mag.append(septable['mag'][row])
        else:
            cat_objs += len(index)
    
    if len(trans_x) == 0:
        print "WARNING: No transients identified"
            
    transients = Table([trans_x, trans_y, trans_a, trans_b, trans_ra ,trans_dec, trans_mag, trans_theta], names=['x', 'y', 'a', 'b', 'ra', 'dec', 'mag', 'theta'])
    return transients, cat_objs

def find_neighbours(transients, pvwcs, r_err, pRA, pDEC, expnum_p):
    
    pX, pY = pvwcs.sky2xy(pRA, pDEC)
    print "  Predicted RA, DEC : {}  {}".format(pRA, pDEC)
    print "  Predicted x, y : {}  {}".format(pX, pY)
    
    found = True
    
    i_list = neighbour_search(transients, pvwcs, r_err, pRA, pDEC)
    if len(i_list) == 0:
        print '-- Expanding radius of nearest neighbour search by 1.5x'
        i_list = neighbour_search(transients, pvwcs, 2*r_err, pRA, pDEC)
        if len(i_list) == 0:
            print 'WARNING: No nearest neighbours were found within {} ++++++++++++++++++++++'.format(r_err*1.5)
            ascii.write(transients, 'asteroid_families/temp_phot_files/{}_phot.txt'.format(expnum_p))
            found = False
                
    return i_list, found
    
def neighbour_search(septable, pvwcs, r_err, pRA, pDEC):    
    '''
    Computes the nearest neighbours to predicted coordinates within an RA/DEC uncertainty circle
    '''
    tree = cKDTree(zip((np.array(septable['x'])).ravel(), (np.array(septable['y'])).ravel()))    
    
    pX, pY = pvwcs.sky2xy(pRA, pDEC) # convert from WCS to pixels
    
    # parse through table and get RA and DEC closest to predicted coordinates (in pixels)
    coords = np.array([pX, pY])
    i_list = tree.query_ball_point(coords, r_err)
    
    return i_list
    
def iden_good_neighbours(expnum, i_list, septable, zeropt, mag_list_jpl, ra_dot, dec_dot, exptime, pvwcs):
    '''
    Selects nearest neighbour object from predicted coordinates as object of interest
    In order:
        Compares measured apparent magnitude to predicted, passes if in range of values
        Calculates eccentricity, passes if greater than minimum value that is inputted
    '''
    
    err = 30.0
    r_pix = exptime * ( (ra_dot)**2 + (dec_dot)**2 )**0.5 / (3600 * 0.184)
    r_pix_err = exptime * (err/100) * ( abs(ra_dot) + abs(dec_dot) ) / (3600 * 0.184)
    assert r_pix_err != 0
    print '  Error allowance is set to {} percent'.format(err)    
        
    mRA_pix = None
    mag_sep_list = []
    index_list = []
    mRA_pix_list = []
    mDEC_pix_list = []
    ra_list = []
    dec_list = []
    theta_list = []
    a_list = []
    b_list = []
    
    mean = np.mean(mag_list_jpl)
    maxmag = np.amax(mag_list_jpl)
    minmag = np.amin(mag_list_jpl)
    
    if ( 2 > maxmag - minmag):
        magrange = 2
    else:
        magrange = maxmag - minmag
    
    print '  Theoretical: {} +/- {}'.format(r_pix, r_pix_err)
    
    something_good = False
    for i in i_list:
        # table format: x, y, a, b, ra, dec, mag
        mag_sep = septable['mag'][i]
        a = septable['a'][i]
        b = septable['b'][i]
        
        r = 2 * ( a**2 - b**2 )**0.5

        if mag_sep > 0:
            something_good = True
            # apply both eccentricity and magnitude conditions
            if (abs(mag_sep - mean) < magrange) & (abs(r-r_pix) < r_pix_err):
                mRA_pix = septable['x'][i]
                mDEC_pix = septable['y'][i]
                
                index_list.append(i)
                mRA_pix_list.append(mRA_pix)
                mDEC_pix_list.append(mDEC_pix)
                mag_sep_list.append(mag_sep) 
                ra_list.append(septable['ra'][i])             
                dec_list.append(septable['dec'][i])
                theta_list.append(septable['theta'][i])
                a_list.append(septable['a'][i])
                b_list.append(septable['b'][i])

                print '  Measured values for index {}: {}, {} {}'.format(i, r, a, b)
            
            elif (abs(r-r_pix) < r_pix_err):
                mRA_pix = septable['x'][i]
                mDEC_pix = septable['y'][i]
                
                index_list.append(i)
                mRA_pix_list.append(mRA_pix)
                mDEC_pix_list.append(mDEC_pix)
                mag_sep_list.append(mag_sep) 
                ra_list.append(septable['ra'][i])             
                dec_list.append(septable['dec'][i])
                theta_list.append(septable['theta'][i])
                a_list.append(septable['a'][i])
                b_list.append(septable['b'][i])
                
                print '  Measured values for index {}: {}, {} {}'.format(i, r, a, b)
                                
                print "WARNING: Magnitude is outside range for index {}: calculated, {} +/- {}, measured, {}".format(i, mean, magrange, mag_sep)
                
            elif (abs(mag_sep - mean) < magrange):
                mRA_pix = septable['x'][i]
                mDEC_pix = septable['y'][i]
                
                index_list.append(i)
                mRA_pix_list.append(mRA_pix)
                mDEC_pix_list.append(mDEC_pix)
                mag_sep_list.append(mag_sep) 
                ra_list.append(septable['ra'][i])             
                dec_list.append(septable['dec'][i])
                theta_list.append(septable['theta'][i])
                a_list.append(septable['a'][i])
                b_list.append(septable['b'][i])
                
                print '  Measured magnitude and predicted: {} {}'.format(mag_sep, mean)
                print "WARNING: Ellipticity is outside range for index {}: calculated, {} +/- {}, measured, {}".format(i, r_pix, r_pix_err, r)
                  
                    
    good_neighbours = Table([index_list, mRA_pix_list, mDEC_pix_list, ra_list, dec_list, mag_sep_list, a_list, b_list, theta_list], 
                names=('index', 'x', 'y', 'ra', 'dec', 'mag', 'a', 'b', 'theta'))
                    
    if something_good == False:
        print "WARNING: Flux of nearest neighbours measured to be 0.0"
        return
    
    if mRA_pix == None:
        print "WARNING: No condition could not be satisfied <<<<<<<<<<<<<<<<<<<<<<<<<<<<"
        print '  Nearest neighbour list: {}'.format(i_list)
        print "  Mag mean, accepted error, and fitting mags: {} {}".format(mean, magrange)
        ascii.write(septable, 'asteroid_families/temp_phot_files/{}_phot.txt'.format(expnum))
        return
        
    else: 
        return good_neighbours, r_pix_err 
        
def check_involvement(objectdata, septable, r_err):
    
    tree = cKDTree(zip((np.array(septable['x'])).ravel(), (np.array(septable['y'])).ravel()))
    search_r = 50
    i_list = tree.query_ball_point([objectdata['x'][0], objectdata['y'][0]], search_r)
    
    print '>> Selects first object in good_neighbours to make polygon from'
    a = objectdata['a'][0]+r_err
    b = objectdata['b'][0]+r_err
    th = objectdata['theta'][0]
    x = objectdata['x'][0]
    y = objectdata['y'][0]
    a_x = a*math.cos(th)
    a_y = a*math.sin(th)
    b_x = b*math.sin(th)
    b_y = b*math.cos(th)
    
    p1 = (x+b_x-a_x, y-a_y+b_y)
    p2 = (x+a_x+b_x, y+a_y+b_y)
    p3 = (x+a_x-b_x, y+a_y-b_y)
    p4 = (x-a_x-b_x, y-a_y-b_y)
    
    print p1, p2, p3, p4
    
    polygon = Polygon([p1, p2, p3, p4])
    
    if len(i_list) > 1:
        print '  Objects within {} pixels of identified asteroid (inc. ast.): {}'.format(search_r, i_list)
        for i in i_list:
            x1 = septable['x'][i]
            y1 = septable['y'][i]
        
            point = Point(x1, y1)
            if polygon.contains(point) == True and x1 != x:
                print '>>> Object is involved with point at x, y: {} {}'.format(x1, y1)
            else:
                print '  Object is not involved.'
    
    else:
        print 'No other objects within {} pixels of identified asteroid'.format(search_r)
        
def check_num_stars(num_objs, size, objectname, expnum_p, username, password, familyname):
        
    # r_old = size * 0.184 / 3600
    # r_new = r_old + 0.05
    
    enough = True
    
    r_old = size * 0.184 / 3600
    r_new = r_old + 0.01
    
    '''if size < 200:
        r_new = 0.02
        r_old = 0.02
    if 200 < size < 500:
        r_new = 0.03
        r_old = 0.02
    if 500 < size < 750:
        r_new = 0.04
        r_old = 0.02
    '''
    
    print '  Number of objects in image: {}'.format(num_objs)
    if 0 < num_objs <= 10:
        r_new = r_old + 0.03
        enough = False
        print '  Not enough stars in the stamp <<<<<<<<<<<<<<<<<<<'              
        get_stamps.get_one_stamp(objectname, expnum_p, r_new, username, password, familyname)
    if 10 < num_objs < 30: 
        r_new = r_old + 0.01
        enough = False
        print '  Not enough stars in the stamp <<<<<<<<<<<<<<<<<<<'                    
        get_stamps.get_one_stamp(objectname, expnum_p, r_new, username, password, familyname)
    
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
                    #index      flux          mag            x             y      
                    outfile.write('{} {} {} {} {} {}\n'.format(
                          objectname, expnum_p, object_data[i][1], object_data[i][2], object_data[i][3], object_data[i][4], 
                          object_data[i][0]))
            except:
                print "ERROR: cannot write to outfile <<<<<<<<<<<<<<<<<<<<<<<<<<"
    
        return object_data
                    
def cut_centered_stamp(familyname, objectname, expnum_p, object_data, r_old, username, password):
        
    if object_data is not None:
        #for i in range(0, len(object_data)):
            # objectname, expnum_p, r_new, RA, DEC, username, password, familyname
        print '-- Cutting recentered stamp'
        get_stamps.centered_stamp(objectname, expnum_p, r_old, object_data[0][5], object_data[0][6], username, password, familyname)                    

def init_dirs(familyname, objectname):
    
    # initiate vos directories 
    global vos_dir
    vos_dir = 'vos:kawebb/postage_stamps/{}'.format(familyname)
    assert exists(vos_dir, force=True)

    
    # initiate local directories
    dir_path = os.path.dirname(os.path.abspath(__file__))
    dir_path_base = '{}/asteroid_families'.format(dir_path)
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

def HMS2deg(ra='', dec=''):
    RA, DEC, rs, ds = '', '', 1, 1
    if dec:
        D, M, S = [float(i) for i in dec.split()]
        if str(D)[0] == '-':
            ds, D = -1, abs(D)
        deg = D + (M/60) + (S/3600)
        DEC = '{0}'.format(deg*ds)
  
    if ra:
        H, M, S = [float(i) for i in ra.split()]
        if str(H)[0] == '-':
            rs, H = -1, abs(H)
        deg = (H*15) + (M/4) + (S/240)
        RA = '{0}'.format(deg*rs)

    if ra and dec:
        return (RA, DEC)
    else:
        return RA or DEC    
    
if __name__ == '__main__':
    main()
    