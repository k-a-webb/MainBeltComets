""" Retrieval of cutouts of the FITS images associated with the OSSOS detections.
    Takes a directory of .ast file (in dbase format) as input """
    
import argparse
import logging
import getpass
import requests
import os

from ossos import mpc
from ossos import storage

BASEURL = "http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/vospace/auth/synctrans"



def cutout(obj, obj_dir, radius, username, password):
    for obs in obj.mpc_observations:  # FIXME: TESTING ONLY #######################################################################
        if obs.null_observation:
            continue
        expnum = obs.comment.frame.split('p')[0]  # only want calibrated images
        # Using the WCS rather than the X/Y, as the X/Y can be unreliable on a long-term basis
        this_cutout = "CIRCLE ICRS {} {} {}".format(obs.coordinate.ra.degree,
                                               obs.coordinate.dec.degree,
                                               radius)
        print this_cutout
        # FIXME: should be able to use line below, but bug in VOSpace requires direct-access workaround, for now.
        # postage_stamp = storage.get_image(expnum, cutout=cutout)

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
        r.raise_for_status()  # confirm the connection worked as hoped
        postage_stamp_filename = "{}_{:11.5f}_{:09.5f}_{:+09.5f}.fits".format(obj.provisional_name,
                                                                              obs.date.mjd,
                                                                              obs.coordinate.ra.degree,
                                                                              obs.coordinate.dec.degree)
        logging.info("{}".format(postage_stamp_filename))
        with open(postage_stamp_filename, 'w') as tmp_file:
            tmp_file.write(r.content)
            storage.copy(postage_stamp_filename, obj_dir + "/" + postage_stamp_filename)
        os.unlink(postage_stamp_filename)  # easier not to have them hanging around



def main():
    
    # IDENTIFY PARAMETERS FOR QUERY OF SSOIS FROM INPUT

    parser = argparse.ArgumentParser(
        description='Parse a directory of TNO .ast files and create links in the postage stamp directory '
                    'that allow retrieval of cutouts of the FITS images associated with the OSSOS detections. '
                    'Cutouts are defined on the WCS RA/DEC of the object position.')

    parser.add_argument("version",
                        help="The OSSOS data release version these stamps should be assigned to.")
    parser.add_argument("--ossin",
                        action="store",
                        default="lixImages.txt",
                        help="The vospace containerNode that clones ossin dbaseclone"
                             "holding the .ast files of astrometry/photometry measurements.")
    parser.add_argument("--blocks", "-b",
                        action="store",
                        default=["o3e", "o3o"],
                        choices=["o3e", "o3o", "O13BL", "Col3N"],
                        help="Prefixes of object designations to include.")
    parser.add_argument("--radius", '-r',
                        action='store',
                        default=0.02,
                        help='Radius (degree) of circle of cutout postage stamp.')
    parser.add_argument("--debug", "-d",
                        action="store_true")
    parser.add_argument("--verbose", "-v",
                        action="store_true")

    args = parser.parse_args()
    
    # NEED PERMISSIONS STILL
    username = raw_input("CADC username: ")
    password = getpass.getpass("CADC password: ")
    
    # From the given input, make list of objects to query
    
    in_file = args.infile
    with open(in_file) as infile: 
        filestr = infile.read()
    input_lines = filestr.split('\n') # array of objects to query
    
    # CONFIRM that the input is the proper format to search for the appropriate ephemeris
    
    for fn in input_lines:
        obj = mpc.MPCReader(args.ossin + fn)  # let MPCReader's logic determine the provisional name
        for block in args.blocks:
            if obj.provisional_name.startswith(block):
                obj_dir = '{}/{}/{}'.format(storage.POSTAGE_STAMPS, args.version, obj.provisional_name)
                if not storage.exists(obj_dir, force=True):
                    storage.mkdir(obj_dir)
                # assert storage.exists(obj_dir, force=True)
                cutout(obj, obj_dir, args.radius, username, password)


if __name__ == '__main__':
    main()