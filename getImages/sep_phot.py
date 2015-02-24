import os
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
    - read in files from VOSpace - NOT WORKING
    - add in psf.py to get PSF
    - find uncertianty of predicted coordinates and use that as nearest neighbour ball search
    - reformat list appending
    - get rid of images with incorrect WCS
    - go through all objects for a family
'''


def find_objects_by_phot(familyname, objectname, ap=10.0, th=5.0, filtertype='r', imagetype='p', elim=0.3):
    
    print 'Finding asteroid {} in family {} '.format(objectname, familyname)
    
    # From the given input, identify the desired filter and rename appropriately
    if filtertype.lower().__contains__('r'):
        filtertype = 'r.MP9601'  # this is the old (standard) r filter for MegaCam
    if filtertype.lower().__contains__('u'):
        filtertype = 'u.MP9301'
    
    # initiate directories
    init_dirs(familyname, objectname)
    out_filename = '{}_r{}_t{}_output.txt'.format(familyname, ap, th)

    with open('{}/{}'.format(stamps_dir, out_filename), 'w') as outfile:
        outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format('Object', "Image", 'flux', 'mag', 'RA', 'DEC', 'ecc'))        

    # get range of magnitudes of object over time span that the images were taken
    print "----- Querying JPL Horizon's ephemeris for apparent magnitudes -----"
    
    mag_list_jpl = parse_mag_jpl(objectname, step=1)
    
    # preform the photometry and identify object
    print "----- Preforming photometry on all images of {} in family {} -----".format(objectname, familyname)
    
    stamp_list = []
    for file in client.listdir(vos_dir):
    #for file in os.listdir(stamps_dir):
        if file.endswith('.fits') == True:
            objectname_file = file.split('_')[0]
            if objectname_file == objectname:
            # images named with convention: object_expnum_RA_DEC.fits
                expnum_p = file.split('_')[1]
                stamp_list.append(expnum_p)
                storage.copy('{}/{}'.format(vos_dir, file), '{}/{}'.format(stamps_dir, file))
                #assert os.path.exists('{}/{}'.format(vos_dir, file))
                #with fits.open('{}/{}'.format(vos_dir, file)) as hdulist:
                with fits.open('{}/{}'.format(stamps_dir, file)) as hdulist:    
                    print " Preforming photometry on image {} ".format(expnum_p)
                
                    if hdulist[0].data is None: # what if more than 2ccd mosaic? could just be aperture values?
                        try:
                            zeropt = fits.getval('{}/{}'.format(stamps_dir, file), 'PHOTZP', 1)
                            exptime = fits.getval('{}/{}'.format(stamps_dir, file), 'EXPTIME', 1)
                            table1 = sep_phot(hdulist[1].data, ap, th)
                            table2 = sep_phot(hdulist[2].data, ap, th)
                            table = vstack([table1, table2])
                            #ascii.write(table, os.path.join(stamps_dir, '{}_phot.txt'.format(expnum_p))) # write all phot data to file in directory familyname/famlyname_objectname/sep_phot_output
                            astheader = hdulist[0].header
                        except LookupError: # maybe not correct error type?
                            print "ERROR: no PHOTZP in header, skipping image "
                    
                    else:
                        try:
                            zeropt = fits.getval('{}/{}'.format(stamps_dir, file), 'PHOTZP', 0)
                            exptime = fits.getval('{}/{}'.format(stamps_dir, file), 'EXPTIME', 0)
                            table = sep_phot(hdulist[0].data, ap, th)
                            astheader = hdulist[0].header
                            #ascii.write(table, os.path.join(stamps_dir, '{}_phot.txt'.format(expnum_p))) # write all phot data to file in directory familyname/famlyname_objectname/sep_phot_output
                        except LookupError:
                            print "ERROR: no PHOTZP in header, skipping image "
                
                    pvwcs = wcs.WCS(astheader)
                    i_list = find_neighbours(table, expnum_p, pvwcs)
                    good_neighbours, mean = iden_good_neighbours(i_list, table, zeropt, mag_list_jpl, elim, pvwcs)
                    print good_neighbours
                    
                    '''
                    thingy = daophot_revised(hdulist, zeropt, filtertype, exptime, good_neighbours, ap, swidth=10, apcor=0.3, maxcount=30000.0)
                    
                    
                    object_data = iden_object(good_neighbours, mean)
                
                    if object_data == None:
                        print "WARNING: Could not identify object {} in image".format(objectname, expnum_p)
                    else:
                        with open('{}/{}'.format(stamps_dir, out_filename), 'a') as outfile:
                            try:
                                outfile.write('{} {} {} {} {} {}\n'.format(
                                          expnum_p, object_data[1], object_data[2], object_data[5], object_data[6], 
                                          object_data[7]))
                            except:
                                print "ERROR: cannot write to outfile"
                    '''
               
                    os.unlink('{}/{}'.format(stamps_dir, file))
    
    if len(stamp_list) == 0:
        print "WARNING: No stamps found"
                
def query_jpl(objectname, step=1, su='d'):
    '''
    Constructs a URL to query JPL Horizon's for apparent magnitude in a date range
    '''
    # from familyname_images.txt get date range of images for objectname
    date_range = []
    with open('{}/{}'.format(family_dir, imageinfo)) as infile:
        for line in infile.readlines()[1:]:
            assert len(line.split()) > 0
            if objectname == line.split()[0]:
                date_range.append(float(line.split()[5]))
                
    if len(date_range) == 0:
        print "WARNING: No images of this object exist"
        assert len(date_range) != 0
        
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
    
    s = '9' # select parameter for apparent magnitude
    
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
            
    return urlData, date_start, date_end

def parse_mag_jpl(objectname, step=1):
    
    urlData, date_start, date_end = query_jpl(objectname, step=1)
    
    # parse through urlData for indexes of start and end dates    
    index_end = None
    index_start = None
    for idx, line in enumerate(urlData):  #testing
        assert line.split() > 0
        try:
            date_jpl = line.split()[0]+' '+(line.split()[1]).strip(',')
            if date_start == date_jpl:
                index_start = idx
        except:
            None
    if index_start is None:
        print "WARNING: index start could not be obtained, set to index 69"
        index_start = 69
    for idx, line in enumerate(urlData):  #testing
        try:
            date_jpl = line.split()[0]+' '+(line.split()[1]).strip(',')
            if date_end == date_jpl:
                index_end = idx
        except:
            None
    if index_end is None:
        print "WARNING: index end could not be obtained, set to index 70"
        index_end = 70
        
    # for indexes from start to end dates, get apparent magnitude values
    mag_list = []
    for line in urlData[index_start:index_end+1]:
        assert len(line.split()) > 0 
        try:
            mag_list.append(float((line.split()[4]).strip(',')))
        except:
            print "WARNING: could not get magnitude from JPL line"
    
    assert len(mag_list) > 0
    return mag_list
       
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
   
    # write to ascii table
    table = Table([objs['x'], objs['y'], flux, objs['a'], objs['b']], names=('x', 'y', 'flux', 'a', 'b'))
    return table

def find_neighbours(septable, expnum_p, pvwcs):
    '''
    Computes the nearest neighbours to predicted coordinates, should eventually find all those in uncertainty ellipse
    '''
    num_neighbours = 25
    if num_neighbours > len(septable):
        num_neighbours = len(septable)
    
    tree = cKDTree(zip((np.array(septable['x'])).ravel(), (np.array(septable['y'])).ravel()))
    
    with open('{}/{}'.format(family_dir, imageinfo)) as infile:
        for line in infile.readlines()[1:]:
            assert len(line.split()) > 0
            expnum_p_fromfile = line.split()[1]
            
            # for entries in *_images.txt that correspond to images of the object
            if expnum_p_fromfile == expnum_p:
                objectname = line.split()[0]
                pRA = float(line.split()[3])
                pDEC = float(line.split()[4])
                expnum = (line.split()[1]).rstrip('p')
                
                pRA_pix, pDEC_pix = pvwcs.sky2xy(pRA, pDEC) # convert from WCS to pixels
                #print " Predicted RA and DEC: {}  {}".format(pRA, pDEC)
                #print "  in pixels: {} {}".format(pRA_pix, pDEC_pix)
                
                # parse through table and get RA and DEC closest to predicted coordinates (in pixels)
                coords = np.array([pRA_pix, pDEC_pix])
                d_list, i_list = tree.query(coords, num_neighbours)
                
                return i_list

def iden_good_neighbours(i_list, septable, zeropt, mag_list_jpl, elim, pvwcs):
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
        
    for i in i_list:
        flux = septable[i][2]
        x = septable[i][0]
        y = septable[i][1]
        a = septable[i][3]
        b = septable[i][4]
        ecc = (1-(b/a)**2)**0.5

        if flux > 0:
            mag_sep = -2.5*math.log10(flux)+zeropt
            mean = np.mean(mag_list_jpl)
            maxmag = np.amax(mag_list_jpl)
            minmag = np.amin(mag_list_jpl)
            all_mag.append(mag_sep)
        
            if ( 2 > maxmag - minmag):
                magrange = 2
            else:
                magrange = maxmag - minmag
                    
            # apply both eccentricity and magnitude conditions
            try:
                if (abs(mag_sep - mean) < magrange) & (ecc > float(elim)):
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
            except:              
                if (abs(mag_sep - mean) < magrange):
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
                    
                    if ecc < 0.5:
                        print "Eccentricity is low: {}".format(ecc)        
    
    # write to ascii table
    good_neighbours = Table([index_list, flux_list, mag_sep_list, x_list, y_list, mRA_pix_list, mDEC_pix_list, ecc_list], 
                names=('index', 'flux', 'mag', 'x', 'y', 'RA', 'DEC', 'ecc'))
                    
    if mRA_pix == None:
        print "WARNING: Magnitude condition could not be satisfied"
        print '  {}'.format(i_list)
        print "  {} {} {}".format(i, mean, all_mag)
        
    else: 
        #print "   Flux, mag: {}, {}".format(flux, mag_sep)     
        mRA, mDEC = pvwcs.xy2sky(mRA_pix, mDEC_pix) # convert from pixels to WCS
        #print " Measured RA and DEC: {}  {}".format(mRA, mDEC)
        #print "  Coordinates: {} {}".format(mRA_pix, mDEC_pix)
        #print " Difference: {} {}".format(mRA - pRA, mDEC - pDEC)
        
        return good_neighbours, mean 
        
def iden_object(good_neighbours, mean):
    
    '''
    mag_diff = []
    for objects in good_neighbours:
        mag = objects[2]
        mag_diff.append(mean - mag)
        
    min_mag = np.amin(mag_diff)
    
    for objects in good_neighbours:
        mag = objects[2]
        if mag == (mean - min_mag):
            the_object = objects
    #print the_object
    '''
    the_object = good_neighbours[0]
    return the_object     
    
def daophot_revised(fits_filename, zeropt, filtertype, exptime, good_neighbours, ap, swidth, apcor, maxcount):
    """
    Compute the centroids and magnitudes of a bunch sources detected on CFHT-MEGAPRIME images.
    Args:
      fits_filename: str  -- The name of the file containing the image to be processed.
    Returns:
      a MOPfiles data structure.
    """
          
    x_in = good_neighbours['x'][0]
    y_in = good_neighbours['y'][0]

    # setup IRAF to do the magnitude/centroid measurements
    iraf.set(uparm="./")
    iraf.digiphot()
    iraf.apphot()
    iraf.daophot(_doprint=0)

    iraf.photpars.apertures = ap
    iraf.photpars.zmag = zeropt
    iraf.datapars.datamin = 0
    iraf.datapars.datamax = maxcount
    iraf.datapars.exposur = ""
    iraf.datapars.itime = exptime
    iraf.fitskypars.annulus = ap + 5
    iraf.fitskypars.dannulus = swidth
    iraf.fitskypars.salgorithm = "mode"
    iraf.fitskypars.sloclip = 5.0
    iraf.fitskypars.shiclip = 5.0
    iraf.centerpars.calgori = "centroid"
    iraf.centerpars.cbox = 5.
    iraf.centerpars.cthreshold = 0.
    iraf.centerpars.maxshift = 3.
    iraf.centerpars.clean = 'no'
    iraf.phot.update = 'no'
    iraf.phot.verbose = 'no'
    iraf.phot.verify = 'no'
    iraf.phot.interactive = 'no'

    # Used for passing the input coordinates
    coofile = tempfile.NamedTemporaryFile(suffix=".coo", delete=False)
    coofile.write('{} {}\n'.format(x_in, y_in))

    # Used for receiving the results of the task
    # mag_fd, mag_path = tempfile.mkstemp(suffix=".mag")
    magfile = tempfile.NamedTemporaryFile(suffix=".mag", delete=False)

    # Close the temp files before sending to IRAF due to docstring:
        # "Whether the name can be used to open the file a second time, while the named temporary file is still open, varies across platforms"
    coofile.close()
    magfile.close()
    os.remove(magfile.name)

    iraf.phot(fits_filename, coofile.name, magfile.name)

    # TODO: Move this filtering downstream to the user.
    phot_filter = "PIER==0 && CIER==0 && SIER==0"

    pdump_out = iraf.pdump(magfile.name, "XCENTER,YCENTER,MAG,MERR,ID,XSHIFT,YSHIFT,LID",
                           phot_filter, header='no', parameters='yes',
                           Stdout=1)

    if not len(pdump_out) > 0:
        mag_content = open(magfile.name).read()
        raise TaskError("photometry failed. {}".format(mag_content))

    os.remove(coofile.name)
    os.remove(magfile.name)

    # setup the mop output file structure
    hdu = {}
    hdu['header'] = {'image': fits_filename,
                     'aper': ap,
                     's_aper': sky,
                     'd_s_aper': swidth,
                     'aper_cor': apcor,
                     'zeropoint': zeropt}
    hdu['order'] = ['X', 'Y', 'MAG', 'MERR', 'ID', 'XSHIFT', 'YSHIFT', 'LID']
    hdu['format'] = {'X': '%10.2f',
                     'Y': '%10.2f',
                     'MAG': '%10.2f',
                     'MERR': '%10.2f',
                     'ID': '%8d',
                     'XSHIFT': '%10.2f',
                     'YSHIFT': '%10.2f',
                     'LID': '%8d'}
    hdu['data'] = {}
    for col in hdu['order']:
        hdu['data'][col] = []

    for line in pdump_out:
        values = line.split()
        for col in hdu['order']:
            if re.match('\%.*f', hdu['format'][col]):
                if col == 'MAG':
                    values[0] = float(values[0]) - float(apcor)
                hdu['data'][col].append(float(values.pop(0)))
            elif re.match('\%.*d', hdu['format'][col]):
                hdu['data'][col].append(int(values.pop(0)))
            else:
                hdu['data'][col].append(values.pop(0))

    # Clean up temporary files generated by IRAF
    os.remove("datistabe.par")
    os.remove("datpdump.par")

    return hdu


def phot_mag(*args, **kwargs):
    """Wrapper around phot which only returns the computed magnitude directly."""
    hdu = phot(*args, **kwargs)

    try:
        return hdu["data"]["X"][0], hdu["data"]["Y"][0], hdu["data"]["MAG"][0], hdu["data"]["MERR"][0],
    except IndexError:
        raise TaskError("Photometry cannot be performed.  "+str(hdu))

    

def init_dirs(familyname, objectname):
    
    global imageinfo
    imageinfo = familyname+'_images_test.txt' # FOR TESTING, CALLS TEST FILE
    print "WARNING: Using test file"
    
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
        time_end = '{}-{}-{} 00:00:00.0'.format(time_end_date[0], time_end_date[1], day)
        
        return time_end
        
        
def main(): 
    
    parser = argparse.ArgumentParser(
                        description='For a given set of fits images, \
                         preforms photometry on a specified object given an aperture size and threshold, \
                        and then selects the object in the image from the predicted coordinates, magnitude, and eventually shape')
    parser.add_argument("--family", '-f',
                        action="store",
                        default=None,
                        help="Asteroid family name. Usually the asteroid number of the largest member.")
    parser.add_argument("--aperture", '-ap',
                        action='store',
                        default=10.0,
                        help='aperture (degree) of circle for photometry.')
    parser.add_argument("--thresh", '-t',
                        action='store',
                        default=5.0,
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
                        default='0.4',
                        help="restrict type of image (unprocessed, reduced, calibrated)")
                            
    args = parser.parse_args()
    
    if args.family == None:
        familyname = 'all'
        all_images = pd.read_table('asteroid_families/all/all_images_test.txt') # FOR TESTING, CALLS TEST FILE
        object_list = all_images['Object'].values
        for objectname in object_list:
            find_objects_by_phot(familyname, str(objectname), float(args.aperture), float(args.thresh), args.filter, args.type, args.elim)
    else:
        if args.object == None:
            object_list = find_family_members(args.family)
            for objectname in object_list:
                find_objects_by_phot(args.family, objectname, float(args.aperture), float(args.thresh), args.filter, args.type, args.elim)
        else:        
            find_objects_by_phot(args.family, args.object, float(args.aperture), float(args.thresh), args.filter, args.type, args.elim)
    
if __name__ == '__main__':
    main()
    