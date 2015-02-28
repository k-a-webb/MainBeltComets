import argparse
import getpass
import requests
import os
import vos
from astropy.time import Time
import urllib2 as url
import numpy as np
from astropy.table import Table, Column

import sys
sys.path.append('/Users/admin/Desktop/MainBeltComets/getImages/ossos_scripts/')

from ossos_scripts import storage
from ossos_scripts import coding
from ossos_scripts import mpc
from ossos_scripts import util



BASEURL = "http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/vospace/auth/synctrans"

"""
Retrieval of cutouts of the FITS images associated with the CFHT/MegaCam detections.
Takes a table (get_images.py output) as input
An Example URL for cutouts from OSSOS (not CFHT/MegaCam)
http://www.canfar.phys.uvic.ca/vospace/auth/synctrans?
TARGET=vos://cadc.nrc.ca~vospace/OSSOS/dbimages/1625356/1625356p.fits&
DIRECTION=pullFromVoSpace&
PROTOCOL=ivo://ivoa.net/vospace/core%23httpget&
view=cutout&
cutout=CIRCLE+ICRS+242.1318+-12.4747+0.05
"""

def main():

    parser = argparse.ArgumentParser(
        description='Parse an familyname.images.txt input file (get_images.py output) and create links in the postage stamp directory '
                    'that allow retrieval of cutouts of the FITS images associated with the CHT/MegaCam detections. '
                    'Cutouts are defined on the WCS RA/DEC of the object position.')
    parser.add_argument("--family", '-f',
                        action="store",
                        default=None,
                        help="The input .txt files of astrometry/photometry measurements.")
    parser.add_argument("--radius", '-r',
                        action='store',
                        default=0.01,
                        help='Radius (degree) of circle of cutout postage stamp.')
    parser.add_argument("--suffix", '-s',
                        action='store',
                        default=None,
                        help='Suffix of mba without family designation')
    args = parser.parse_args()
    
    # CADC PERMISSIONS
    username = raw_input("CADC username: ")
    password = getpass.getpass("CADC password: ")
    
    get_stamps(args.family, username, password, args.radius, args.suffix)
    
def get_stamps(familyname, username, password, radius=0.01, suffix=None):
    
    print "----- Cutting postage stamps of objects in family {}  from CFHT/MegaCam images -----".format(familyname)	                  
    
    dir_path_base = '/Users/admin/Desktop/MainBeltComets/getImages/asteroid_families'
    family_dir = os.path.join(dir_path_base, familyname)
    if os.path.isdir(family_dir) == False:
        print "Invalid family name or directory does not exist"

    image_list = '{}/{}_images.txt'.format(family_dir, familyname)
    with open(image_list) as infile: 
        for line in infile.readlines()[9360:]: # skip header info
            assert len(line.split()) > 0
            objectname = line.split()[0]
            expnum = line.split()[1]
            RA = float(line.split()[3])   
            DEC = float(line.split()[4])
                        
            vos_dir = 'vos:kawebb/postage_stamps/{}'.format(familyname)
            if not storage.exists(vos_dir, force=True):
                storage.mkdir(vos_dir)
            #assert storage.exists(vos_dir, force=True)

            postage_stamp_filename = "{}_{}_{:8f}_{:8f}.fits".format(objectname, expnum, RA, DEC)
            if storage.exists('{}/{}'.format(vos_dir, postage_stamp_filename)) == True:
                print "  Stamp already exists"
            else:
                
                urlData, date_start, date_end = query_jpl(familyname, objectname, step=1)
                print "----- Querying JPL Horizon's ephemeris for RA and DEC uncertainties -----"
                RA_3sigma, DEC_3sigma = parse_mag_jpl(urlData, date_start, date_end) # in arcseconds
            
                RA_3sigma_avg = np.mean(RA_3sigma) / 3600 # convert to degrees
                DEC_3sigma_avg = np.mean(DEC_3sigma) / 3600
            
                if RA_3sigma_avg > DEC_3sigma_avg:
                    r_temp = RA_3sigma_avg
                else:
                    r_temp = DEC_3sigma_avg
                
                if r_temp > radius:
                    radius = r_temp
                
                cutout(objectname, expnum, RA, DEC, radius, username, password, familyname)
                
def get_one_stamp(objectname, expnum, radius, username, password, familyname):
    
    print "-- Cutting postage stamps of {} {}".format(objectname, expnum)	                  
    
    dir_path_base = '/Users/admin/Desktop/MainBeltComets/getImages/asteroid_families'
    family_dir = os.path.join(dir_path_base, familyname)
    if os.path.isdir(family_dir) == False:
        print "Invalid family name or directory does not exist"

    image_list = '{}/{}_images.txt'.format(family_dir, familyname)
    with open(image_list) as infile: 
        for line in infile.readlines()[1:]: # skip header info
            assert len(line.split()) > 0
            objectname = line.split()[0]
            expnum_file = line.split()[1]
            RA = float(line.split()[3])   
            DEC = float(line.split()[4])
                        
            vos_dir = 'vos:kawebb/postage_stamps/{}'.format(familyname)
            if not storage.exists(vos_dir, force=True):
                storage.mkdir(vos_dir)
            #assert storage.exists(vos_dir, force=True)
            
            if expnum == expnum_file:
                
                postage_stamp_filename = "{}_{}_{:8f}_{:8f}.fits".format(objectname, expnum, RA, DEC)
                storage.remove('vos:kawebb/postage_stamps/{}/{}'.format(familyname, postage_stamp_filename))
                
                cutout(objectname, expnum, RA, DEC, radius, username, password, familyname)
                return                
	
def centered_stamp(objectname, expnum,radius, RA, DEC, username, password, familyname):
    
    vos_dir = 'vos:kawebb/postage_stamps/{}'.format(familyname)
    postage_stamp_filename = "{}_{}_{:8f}_{:8f}_centered.fits".format(objectname, expnum, RA, DEC)
    if storage.exists('{}/{}'.format(vos_dir, postage_stamp_filename)) == True:
        print "  Stamp already exists"
    else:
        cutout(objectname, expnum, RA, DEC, radius, username, password, familyname)
    return
    
def cutout(objectname, image, RA, DEC, radius, username, password, familyname):
    
    ''' 
    Test for image known to work   
    image = '1667879p'
    RA = 21.1236333333
    DEC = 11.8697277778
    '''
    vos_dir = 'vos:kawebb/postage_stamps/{}'.format(familyname)
    output_dir = 'asteroid_families/{}/{}_stamps'.format(familyname, familyname)
    if os.path.isdir(output_dir) == False:
        os.makedirs(output_dir)
        
    this_cutout = "CIRCLE ICRS {} {} {}".format(RA, DEC, radius)                                 
    print "cut out: {} {} {} {} {}".format(objectname, image, RA, DEC, radius)

    expnum = image.split('p')[0] # only want calibrated images    
    target = storage.vospace.fixURI(storage.get_uri(expnum)) 
    direction = "pullFromVoSpace"
    protocol = "ivo://ivoa.net/vospace/core#httpget"
    view = "cutout"
    params = {"TARGET": target,
                  "PROTOCOL": protocol,
                  "DIRECTION": direction,
                  "cutout": this_cutout,
                  "view": view}
    
    r = requests.get(BASEURL, params=params, auth=(username, password))
    
    try:
          r.raise_for_status()  # confirm the connection worked as hoped
    
          postage_stamp_filename = "{}_{}_{:8f}_{:8f}.fits".format(objectname, image, RA, DEC)

          with open('{}/{}'.format(output_dir, postage_stamp_filename), 'w') as tmp_file:
              object_dir = 'asteroid_families/{}/{}_stamps/{}'.format(familyname, familyname, postage_stamp_filename)
              assert os.path.exists(object_dir)
              tmp_file.write(r.content)
              storage.copy(object_dir, '{}/{}'.format(vos_dir, postage_stamp_filename))
          #os.unlink(object_dir)  # easier not to have them hanging around    
    
    except: 
        print 'Connection Failed'
        return
    

def query_jpl(familyname, objectname, step=1, su='d'):
    '''
    Constructs a URL to query JPL Horizon's for RA and DEC uncertainties
    '''
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

    s = '36' # select parameter for RA and DEC uncertainty

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

def parse_mag_jpl(urlData, date_start, date_end):
        
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
    RA_3sigma = []
    DEC_3sigma = []
    for line in urlData[index_start:index_end+1]:
        assert len(line.split()) > 0 
        try:
            RA_3sigma.append(float((line.split()[4]).strip(',')))
            DEC_3sigma.append(float((line.split()[5]).strip(',')))
        except:
            print "WARNING: could not get magnitude from JPL line"
    
    assert len(RA_3sigma) > 0
    return RA_3sigma, DEC_3sigma

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
    	
if __name__ == '__main__':
    main()	
		
		
		
		
		