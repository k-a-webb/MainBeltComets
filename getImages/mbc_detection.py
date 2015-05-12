__author__ = 'kw'

"""
Detects the deviation of an asteroid's PSF profile from a stellar model PSF profile

- Reads in the SEP ouput for the asteroid from a text file in _PHOT_DIR
- Pulls the postage stamp of the asteroid from VOSpace and reads in fits header and data
- Pulls the stellar model PSF data from VOSpace for the exposure
    - if no PSF, stop and write object name and expnum to a file _OUTPUT_NO_MKPSF
- If the object is too bright, <18.6 mag, stop and write object name and expnum to a file _OUTPUT_TOO_BRIGHT
    - this is becase the PSF is saturated
- Pulls the file with the fwhm value for the exposure from VOSpace
- Preform photometry on the postage stamp to get background value and flux/fluxerr of the asteroid
- Calculates the asteroid PSF profile in get_asteroid_data
- Calculates the stellar model PSF profile in build_star_psf
- Creates a fake comet psf from the stellar model PSF profile
    - to use instead of asteroid profile in the comparison, specify test = True
- Compares the asteroid PSF profile to the stellar model PSF profile and checks for data points which deviate
  from the mean by more than 3,4,5 sigma
"""

activate_this = '/Users/admin/Desktop/MainBeltComets/bin/activate_this.py'
execfile(activate_this, dict(__file__=activate_this))
import pandas as pd
import numpy as np
import vos
import argparse
import os
import sep
from astropy.io import fits
from scipy.ndimage.interpolation import rotate
import math
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import sys
from pyraf import iraf

pd.set_option('display.max_columns', 30)
client = vos.Client()
sys.path.append('User/admin/Desktop/OSSOS/MOP/src/ossos-pipeline/ossos')
from ossos import storage

_VOS_DIR = 'vos:kawebb/postage_stamps'
_DIR_PATH_BASE = os.path.dirname(os.path.abspath(__file__))
_STAMPS_DIR = '{}/postage_stamps'.format(_DIR_PATH_BASE)
if not os.path.exists(_STAMPS_DIR):
    os.mkdir(_STAMPS_DIR)
_OSSOS_PATH_BASE = 'vos:OSSOS/dbimages'
_PHOT_DIR = '{}/phot_output'.format(_DIR_PATH_BASE)
_OUTPUT_DIR = '{}/psf_output'.format(_DIR_PATH_BASE)
if not os.path.exists(_OUTPUT_DIR):
    os.mkdir(_OUTPUT_DIR)

_CCD = 'EXTNAME'
_ZMAG = 'PHOTZP'
_TARGET = 'OBJECT'

_OBJECT_HEADER = 'object'
_EXPNUM_HEADER = 'expnum'
_MAG_HEADER = 'mag'
_XMIN_HEADER = 'xmin'
_XMAX_HEADER = 'xmax'
_YMIN_HEADER = 'ymin'
_YMAX_HEADER = 'ymax'
_XMID_HEADER = 'x_mid'
_YMID_HEADER = 'y_mid'
_THETA_HEADER = 'theta'
_A_HEADER = 'a'
_B_HEADER = 'b'

_INPUT_HEADERS = [_OBJECT_HEADER, _EXPNUM_HEADER, 'p_ra', 'p_dec', 'ra', 'dec', 'p_x', 'p_y', 'x', 'y', _XMID_HEADER,
                  _YMID_HEADER, _XMIN_HEADER, _XMAX_HEADER, _YMIN_HEADER, _YMAX_HEADER, _A_HEADER, _B_HEADER,
                  _THETA_HEADER, 'p_f', 'f', 'consistent_f', 'p_mag', _MAG_HEADER, 'flux', 'consistent_mag',
                  'date_time']

_BUFFER1 = 1.5  # aperture of asteroid is this * the fwhm.
_BUFFER2 = 20.  # the size of the cutout of the asteroid before rotation.
SPINE_COLOR = 'gray'  # for latexifcation

_OUTPUT_NO_MKPSF = 'no_image_psf.txt'.format(_DIR_PATH_BASE)
_OUTPUT_TOO_BRIGHT = 'too_bright.txt'.format(_DIR_PATH_BASE)
_INPUT_FILE = 'output.txt'
if not os.path.exists(_OUTPUT_NO_MKPSF):
    open(_OUTPUT_NO_MKPSF, 'w').close()
if not os.path.exists(_OUTPUT_TOO_BRIGHT):
    open(_OUTPUT_TOO_BRIGHT, 'w').close()

'''
headers in {}_output:
object expnum p_ra p_dec ra dec p_x p_y x y x_mid y_mid xmin xmax ymin ymax a b theta p_f f consistent_f p_mag mag
flux consistent_mag date_time
'''


def main():
    parser = argparse.ArgumentParser(
        description='For a given set of fits images: preforms photometry, identifies a specific object, and returns \
                        the orbital elements of the object')
    parser.add_argument("--family", '-f',
                        action="store",
                        default='all',
                        help="Asteroid family name. Usually the asteroid number of the largest member.")
    parser.add_argument("--object", '-o',
                        action='store',
                        default=None,
                        help='The object to be identified.')
    parser.add_argument("--expnum", '-x',
                        action='store',
                        default=None,
                        help='The expsure in which the object is being identified.')
    parser.add_argument("--test",
                        action='store',
                        default=False,
                        help='Use a test comet profile.')

    args = parser.parse_args()

    phot_file_path = '{}/{}/{}_{}'.format(_PHOT_DIR, args.family, args.family, _INPUT_FILE)
    assert os.path.exists(phot_file_path), 'Phot output file ({}) does not exist'.format(phot_file_path)

    phot_table = pd.read_table(phot_file_path, dtype={_OBJECT_HEADER: object}, index_col=False, header=None, sep=' ',
                               names=_INPUT_HEADERS)

    if args.object is None:
        for i in range(len(phot_table)):
            detect_mbc(args.family, phot_table[_OBJECT_HEADER][i], phot_table[_EXPNUM_HEADER][i], phot_table, args.test)
            print '\n'

    elif args.expnum is None:
        object_table = phot_table.query('object == "{}"'.format(args.object))
        object_table.reset_index(drop=True, inplace=True)
        for i in range(len(object_table)):
            detect_mbc(args.family, args.object, object_table[_EXPNUM_HEADER][i], object_table, args.test)
            print '\n'
    else:
        detect_mbc(args.family, args.object, args.expnum, phot_table, args.test)


def detect_mbc(family_name, object_name, expnum, phot_table, test):
    """
    Compare psf of asteroid with mean of stars to detect possible activity
    """

    # read in asteroid identification values from the photometry output
    asteroid_id = phot_table.query(
        '{} == "{}" & {} == "{}"'.format(_OBJECT_HEADER, object_name, _EXPNUM_HEADER, expnum))
    print asteroid_id
    assert len(asteroid_id) == 1, 'No object or multiple objects identified'

    # read in postage stamp header and data, do photometry to measure background (needed for saturation check)
    header, exp_data, fits_file = fits_data(object_name, expnum, family_name)

    # make sure that a mean star psf has been created form the OSSOS pipeline
    if not storage.get_status(expnum.strip('p'), header[_CCD].split('d')[1], 'mkpsf'):
        print '>> PSF does not exist'
        ast_sky_psf = build_ast_profile(asteroid_id, exp_data, fwhm, family_name) # to build cutout of object
        write_no_mkpsf(family_name, '{} {}'.format(expnum.strip('p'), header[_CCD].split('d')[1]))
        return

    # reject any object too bright that will definetly be saturated
    mag = asteroid_id[_MAG_HEADER].values[0]
    if mag < 18.5:
        print '>> Object is too bright for accurate photometry'
        ast_sky_psf = build_ast_profile(asteroid_id, exp_data, fwhm, family_name) # to build cutout of object
        write_too_bright(family_name, asteroid_id)
        return

    # get fwhm from OSSOS VOSpace file
    fwhm = storage.get_fwhm(expnum.strip('p'), header[_CCD].split('d')[1])
    bkg, flux, fluxerr = sep_phot(exp_data, asteroid_id)

    ast_sky_psf = build_ast_profile(asteroid_id, exp_data, fwhm, family_name)
    ast_psf = np.subtract(ast_sky_psf, bkg)
    try:
        star_psf = build_star_profile(ast_psf, expnum, header, asteroid_id, fwhm, flux)
    except Exception, e:
        print 'Error calculating star psf: {}'.format(e)
        return

    if test is not False:
        comet_psf = build_comet_psf(star_psf, fwhm)
        detection, sig = compare_psf(star_psf, comet_psf, comet_psf, family_name, test)
        if detection:
            print '>> Detect possible comae <<'

    else:
        print '-- Comparing PSFs'
        detection, sig = compare_psf(star_psf, ast_psf, ast_sky_psf, family_name, expnum)
        if detection:
            print '>> Detect possible comae <<'
            write_to_file(asteroid_id, family_name, sig)
        else:
            write_to_no_detection_file(asteroid_id, family_name)


def fits_data(object_name, expnum, family_name):
    """
    Creates local copy of fits file from VOSpace
    """

    if family_name == 'none':
        vos_dir = '{}/none'.format(_VOS_DIR)
    else:
        vos_dir = '{}/all'.format(_VOS_DIR)

    assert storage.exists(vos_dir), 'Vos directory does not exist, or permissions have expired'
    for fits_file in client.listdir(vos_dir):  # images named with convention: object_expnum_RA_DEC.fits

        if fits_file.endswith('.fits'):
            objectname_file = fits_file.split('_')[0]
            expnum_file = fits_file.split('_')[1]
            if (expnum_file == expnum) and (objectname_file == object_name):
                file_path = '{}/{}'.format(_STAMPS_DIR, fits_file)
                storage.copy('{}/{}'.format(vos_dir, fits_file), file_path)

                with fits.open(file_path) as hdulist:
                    data = hdulist[0].data
                    header = hdulist[0].header

                os.unlink(file_path)
                return header, data, fits_file


def sep_phot(exp_data, asteroid_id, ap=10.0):
    """
    Measure background of postage stamp and the flux_err of the asteroid
    """

    data2 = np.ones(exp_data.shape) * exp_data
    # np.copyto(data2, exp_data)
    try:
        bkg = sep.Background(data2)
    except ValueError:
        data3 = data2.byteswap(True).newbyteorder()
        bkg = sep.Background(data3)

    # Directly subtract the background from the data in place
    bkg.subfrom(data2)

    # calculate the Kron radius for each object, then we perform elliptical aperture photometry within that radius
    kronrad, krflag = sep.kron_radius(data2, asteroid_id[_XMID_HEADER], asteroid_id[_YMID_HEADER],
                                      asteroid_id[_A_HEADER],
                                      asteroid_id['b'], asteroid_id[_THETA_HEADER], ap)
    flux, fluxerr, flag = sep.sum_ellipse(data2, asteroid_id[_XMID_HEADER], asteroid_id[_YMID_HEADER],
                                          asteroid_id[_A_HEADER], asteroid_id[_B_HEADER],
                                          asteroid_id[_THETA_HEADER], 2.5 * kronrad, subpix=1, err=bkg.globalrms)

    return bkg.globalback, flux, fluxerr


def build_ast_profile(object_data, data, fwhm, family_name):
    """
    Calculate psf of asteroid, taking into acount trailing effect
    Cut a square around the polygon of the object, remove saturated rows, rotate so that the ellipse is parallel to the
    horizontal axis, then mask shape of ellipse around the object
    """

    # define parameters of the square to cut out from the polygon around the elongated object
    x_min = int(object_data[_XMIN_HEADER].values[0] - _BUFFER2)
    x_max = int(object_data[_XMAX_HEADER].values[0] + _BUFFER2)
    y_min = int(object_data[_YMIN_HEADER].values[0] - _BUFFER2)
    y_max = int(object_data[_YMAX_HEADER].values[0] + _BUFFER2)

    # ensure that you are selecting data that is on the CCD
    data_y_max, data_x_max = data.shape
    if x_min < 0:
        x_min = 0
    if y_min < 0:
        y_min = 0
    if x_max > data_x_max:
        x_max = data_x_max
    if y_max > data_y_max:
        y_max = data_y_max

    data_obj = np.ones((len(range(y_min, y_max)), len(range(x_min, x_max))))
    np.copyto(data_obj, data[y_min:y_max, x_min:x_max])

    # rotate the data about the angle of elongation, semi major axis is along x as reading out rows (not columns)
    data_rot = rotate(data_obj, math.degrees(object_data[_THETA_HEADER].values[0]))

    # cut out a circular aperture around the object
    a = 0.5 * math.sqrt((object_data[_XMAX_HEADER].values[0] - object_data[_XMIN_HEADER].values[0]) ** 2 +
                        (object_data[_YMAX_HEADER].values[0] - object_data[_YMIN_HEADER].values[0]) ** 2)
    ell_buffer = _BUFFER1 * fwhm
    cy, cx = np.divide(data_rot.shape, 2)

    mask = np.ones(data_rot.shape, dtype=bool)
    mask2 = np.zeros(data_rot.shape)

    x_range = range(int(cx - a - ell_buffer + 1), int(cx + a + ell_buffer - 1))

    if x_range[0] < 0 and x_range[-1] > data_rot.shape[1]:
        x_range = range(0, data_rot.shape[1])
    elif x_range[0] < 0:
        x_range = range(0, int(cx + a + ell_buffer - 1))
    elif x_range[-1] > data_rot.shape[1]:
        x_range = range(int(cx - a - ell_buffer + 1), data_rot.shape[1])

    for x in x_range:
        y1 = int(cy - math.sqrt((a + ell_buffer) ** 2 - (x - cx) ** 2))
        y2 = int(cy + math.sqrt((a + ell_buffer) ** 2 - (x - cx) ** 2))
        for y in range(y1, y2):
            mask[y][x] = False
            mask2[y][x] = 1

    data_cutout = np.ma.array(data_rot, mask=mask)

    hdu = fits.PrimaryHDU()
    hdu.data = np.multiply(data_rot, mask2)
    hdu.writeto('psf_output/{}/cutout_{}_{}.fits'.format(family_name, object_data[_OBJECT_HEADER].values[0],
                                                         object_data[_EXPNUM_HEADER].values[0]), clobber=True)

    # take the mean of all the values in each ROW
    data_mean = np.ma.mean(data_cutout, axis=1)
    data_nonzero = data_mean[np.nonzero(np.ma.fix_invalid(data_mean, fill_value=0))[0]]
    assert len(data_nonzero) > 0, "All data in cutout is zero"

    return data_nonzero


def build_star_profile(data_ast, expnum, header, asteroid_id, fwhm, flux):
    """
    Center the star data on the peak of the asteroid psf, and interpolate data points
    """

    mag = (-2.5 * np.log10(flux) + header[_ZMAG])[0]

    try:
        data_str = get_star_data(asteroid_id, mag, expnum, header)
    except:
        print '>> ERROR: Iraf seepsf error'
        raise Exception

    x_ast = range(len(data_ast))
    x_str = range(len(data_str))
    p_str = [np.amax(data_str), float(len(x_str)) / 2, fwhm, 0.]
    p_ast = [np.amax(data_ast), float(len(x_ast)) / 2., fwhm, 0.]

    # shift the center of the star psf to the center of the asteroid psf
    try:
        # fit a gaussian to the asteroid psf and then the star psf, will shift from fit peak x-value
        fitp_ast, fitco_ast = curve_fit(gauss, x_ast, data_ast, p_ast)
        fitp_str, fitco_str = curve_fit(gauss, x_str, data_str, p_str)
        gauss_shift = fitp_str[1] - fitp_ast[1]
    except RuntimeError:  # if gaussian fit cannot be obtained in 100 iterations
        gauss_shift = np.argmax(data_str) - np.argmax(data_ast)

    # interpolate the psf values of the star at the x values that there is data for the asteroid psf
    y_str = lambda x: np.interp(x, np.subtract(x_str, gauss_shift), data_str)
    data_str_at_astpts = y_str(x_ast)

    return data_str_at_astpts


def get_star_data(asteroid_id, mag, expnum, header):
    """
    From ossos psf fitted image, calculate mean of the flux of each row of the rotated PSF
    """

    # calculate mean psf
    uri = storage.get_uri(expnum.strip('p'), header[_CCD].split('d')[1])
    ossos_psf = '{}.psf.fits'.format(uri.strip('.fits'))
    local_psf = '{}{}.psf.fits'.format(expnum, header[_CCD].split('d')[1])
    local_file_path = '{}/{}'.format(_STAMPS_DIR, local_psf)
    storage.copy(ossos_psf, local_file_path)

    # pvwcs = wcs.WCS(header)
    # x, y = pvwcs.sky2xy(asteroid_id['ra'].values, asteroid_id['dec'].values)
    x = asteroid_id[_XMID_HEADER].values[0]
    y = asteroid_id[_YMID_HEADER].values[0]

    # run seepsf on the mean psf image
    iraf.set(uparm="./")
    iraf.digiphot(_doprint=0)
    iraf.apphot(_doprint=0)
    iraf.daophot(_doprint=0)
    iraf.seepsf(local_file_path, local_psf, xpsf=x, ypsf=y, magnitude=mag)

    with fits.open(local_psf) as hdulist:
        data = hdulist[0].data

    th = math.degrees(asteroid_id[_THETA_HEADER].values[0])
    data_rot = rotate(data, th)
    data_rot = np.ma.masked_where(data_rot == 0, data_rot)

    data_mean = np.ma.mean(data_rot, axis=1)

    os.unlink(local_psf)
    os.unlink(local_file_path)

    return data_mean[np.nonzero(np.ma.fix_invalid(data_mean, fill_value=0))[0]]


def compare_psf(star_psf, ast_psf, ast_sky_psf):
    """
    Compare psf of asteroid against mean of stars, check if anomaly in wings
    >> changed data to not normalized nor baseline subtracted
    """

    x_ast = range(len(ast_psf))

    # star_psf_sig = np.sqrt(np.absolute(star_psf))
    star_psf_sig = 0  # Assume uncertainty on the stellar model is zero
    ast_psf_sig = np.sqrt(np.absolute(ast_sky_psf))
    r = np.divide(ast_psf, star_psf)
    r_sig = r * np.sqrt((np.divide(ast_psf_sig, ast_psf) ** 2 + np.divide(star_psf_sig, star_psf) ** 2))
    r_mean = np.ma.mean(np.sort(r)[4:-5])  # don't include outliers in calculation of the mean
    ratio = (r - r_mean) / r_sig

    #peak = np.argmax(ast_psf)
    #r_mean_peak = np.mean(np.divide(ast_psf[peak - 2: peak + 1], star_psf[peak - 2: peak + 1]))
    #ratio_peak = (r - r_mean_peak) / r_sig

    # print r
    # print r_sig
    # print r_mean
    print '>> (r - r_mean) / r_sig'
    print ratio
    #print ratio_peak

    with sns.axes_style('ticks'):
        # latexify()
        plt.plot(x_ast, ast_psf, label='Asteroid', ls='-')
        plt.plot(x_ast, star_psf, label='Stellar Model', ls='--')
        plt.ylabel('Mean Flux')
        plt.legend()
        plt.show()

    if len(ratio[np.greater_equal(np.absolute(ratio), 5)]) > 2:
        return True, 5
    elif len(ratio[np.greater_equal(np.absolute(ratio), 4)]) > 2:
        return True, 4
    elif len(ratio[np.greater_equal(np.absolute(ratio), 3)]) > 2:
        return True, 3
    else:
        return False, 0


def write_to_file(asteroid_id, family_name, sig):
    with open('{}/{}/abv_{}_sig.txt'.format(_OUTPUT_DIR, family_name, sig), 'a') as outfile:
        outfile.write('{} {}\n'.format(asteroid_id[_OBJECT_HEADER].values[0], asteroid_id[_EXPNUM_HEADER].values[0]))


def write_to_no_detection_file(asteroid_id, family_name):
    with open('{}/{}/no_detection.txt'.format(_OUTPUT_DIR, family_name), 'a') as outfile:
        outfile.write('{} {}\n'.format(asteroid_id[_OBJECT_HEADER].values[0], asteroid_id[_EXPNUM_HEADER].values[0]))


def write_no_mkpsf(family_name, expnum_ccd):
    with open('{}/{}/{}'.format(_OUTPUT_DIR, family_name, _OUTPUT_NO_MKPSF), 'r') as outfile:
        contents = []
        for line in outfile:
            contents.append(line.strip('\n'))
    if expnum_ccd not in contents:
        with open('{}/{}/{}'.format(_OUTPUT_DIR, family_name, _OUTPUT_NO_MKPSF), 'a') as outfile:
            outfile.write('{}\n'.format(expnum_ccd))


def write_too_bright(family_name, asteroid_id):
    object_expnum = []
    asteroid_expnum = '{} {}'.format(asteroid_id[_OBJECT_HEADER].values[0], asteroid_id[_EXPNUM_HEADER].values[0])
    with open('{}/{}/{}'.format(_OUTPUT_DIR, family_name, _OUTPUT_TOO_BRIGHT), 'r') as outfile:
        for line in outfile:
            object_expnum.append('{} {}'.format(line.split(' ')[0], line.split(' ')[0]))
    if asteroid_expnum not in object_expnum:
        with open('{}/{}/{}'.format(_OUTPUT_DIR, family_name, _OUTPUT_TOO_BRIGHT), 'a') as outfile:
            outfile.write(
                '{} {} {}\n'.format(asteroid_id[_OBJECT_HEADER].values[0], asteroid_id[_EXPNUM_HEADER].values[0],
                                    asteroid_id[_MAG_HEADER].values[0]))


def latexify(fig_width=None, fig_height=None, columns=1):
    """Set up matplotlib's RC params for LaTeX plotting.
    Call this before plotting a figure.

    Parameters
    ----------
    fig_width : float, optional, inches
    fig_height : float,  optional, inches
    columns : {1, 2}
    """

    # code adapted from http://www.scipy.org/Cookbook/Matplotlib/LaTeX_Examples

    # Width and max height in inches for IEEE journals taken from
    # computer.org/cms/Computer.org/Journal%20templates/transactions_art_guide.pdf

    assert (columns in [1, 2])

    if fig_width is None:
        fig_width = 5 if columns == 1 else 6.9  # width in inches

    if fig_height is None:
        golden_mean = (math.sqrt(5) - 1.0) / 2.0  # Aesthetic ratio
        fig_height = fig_width * golden_mean  # height in inches

    max_height_inches = 8.0
    if fig_height > max_height_inches:
        print("WARNING: fig_height too large:" + fig_height +
              "so will reduce to" + max_height_inches + "inches.")
        fig_height = max_height_inches

    params = {'backend': 'ps',
              'text.latex.preamble': ['\usepackage{gensymb}'],
              'axes.labelsize': 8,  # fontsize for x and y labels (was 10)
              'axes.titlesize': 8,
              'font.size': 10,  # was 10
              'legend.fontsize': 8,  # was 10
              'xtick.labelsize': 8,
              'ytick.labelsize': 8,
              'text.usetex': True,
              'figure.figsize': [fig_width, fig_height],
              'font.family': 'serif'
              }

    matplotlib.rcParams.update(params)


def gauss(x, *p):
    amp, mu, sigma, b = p
    return amp * np.exp(-(x - mu) ** 2 / (2. * sigma ** 2)) + b


def build_comet_psf(data_str, fwhm):
    c = np.argmax(data_str)
    x_min = int(c - 1.5 * fwhm)
    x_max = int(c + 1.5 * fwhm + 1)

    length = 2 * 2 * fwhm
    b = np.mean(np.sort(data_str)[:5])
    amp = np.amax(data_str) - b

    coma = []
    for x in range(len(data_str)):
        if x_min <= x <= x_max:
            y = 0.9 * amp * np.sin(np.divide(x - x_min, length) * np.pi)
            if y > (data_str[x] - b):
                coma.append(y - data_str[x] + b)
            else:
                coma.append(0)
        else:
            coma.append(0)

    return np.add(data_str, coma)


if __name__ == '__main__':
    main()