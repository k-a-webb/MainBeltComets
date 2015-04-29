__author__ = 'kw'

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
from scipy.optimize import curve_fit
from scipy.stats import ttest_ind
from scipy.stats import ttest_1samp
import sys
from pyraf import iraf

pd.set_option('display.max_columns', 30)
client = vos.Client()
sys.path.append('User/admin/Desktop/OSSOS/MOP/src/ossos-pipeline/ossos')
from ossos import storage

_VOS_DIR = 'vos:kawebb/postage_stamps'
_DIR_PATH_BASE = os.path.dirname(os.path.abspath(__file__))
_STAMPS_DIR = '{}/postage_stamps'.format(_DIR_PATH_BASE)
_OSSOS_PATH_BASE = 'vos:OSSOS/dbimages'
_PHOT_DIR = '{}/phot_output'.format(_DIR_PATH_BASE)
_OUTPUT_DIR = '{}/psf_output'.format(_DIR_PATH_BASE)

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

_INPUT_HEADERS = [_OBJECT_HEADER, _EXPNUM_HEADER, 'p_ra', 'p_dec', 'ra', 'dec', 'p_x', 'p_y', 'x', 'y', _XMID_HEADER,
                  _YMID_HEADER, _XMIN_HEADER, _XMAX_HEADER, _YMIN_HEADER, _YMAX_HEADER, 'a', 'b', _THETA_HEADER, 'p_f',
                  'f', 'consistent_f', 'p_mag', _MAG_HEADER, 'flux', 'consistent_mag', 'date_time']

_BUFFER1 = 2.5  # aperture of asteroid is this * the fwhm.
_BUFFER2 = 30.  # the size of the cutout of the asteroid before rotation.

_OUTPUT_NO_MKPSF = 'no_image_psf.txt'.format(_DIR_PATH_BASE)
_OUTPUT_TOO_BRIGHT = 'too_bright.txt'.format(_DIR_PATH_BASE)
_INPUT_FILE = 'output.txt'

SPINE_COLOR = 'gray'

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

    args = parser.parse_args()

    phot_table = pd.read_table('{}/{}/{}_{}'.format(_PHOT_DIR, args.family, args.family, _INPUT_FILE),
                               dtype={'object': object}, index_col=False, header=None, sep=' ',
                               names=_INPUT_HEADERS)

    if args.object is None:
        for i in range(3, len(phot_table)):
            detect_mbc(args.family, phot_table['object'][i], phot_table['expnum'][i], phot_table)
            print '\n'

    elif args.expnum is None:
        object_table = phot_table.query('object == "{}"'.format(args.object))
        object_table.reset_index(drop=True, inplace=True)
        for i in range(len(object_table)):
            detect_mbc(args.family, args.object, object_table['expnum'][i], object_table)
            print '\n'
    else:
        detect_mbc(args.family, args.object, args.expnum, phot_table)


def detect_mbc(family_name, object_name, expnum, phot_table):
    """
    Compare psf of asteroid with mean of stars to detect possible activity
    """

    # read in asteroid identification values from the photometry output
    asteroid_id = phot_table.query('object == "{}" & expnum == "{}"'.format(object_name, expnum))
    print asteroid_id
    assert len(asteroid_id) == 1, 'Multiple object identified'

    # read in postage stamp header and data, do photometry to measure background (needed for saturation check)
    header, data, fits_file = fits_data(object_name, expnum, family_name)

    # make sure that a mean star psf has been created form the OSSOS pipeline
    if not storage.get_status(expnum.strip('p'), header[_CCD].split('d')[1], 'mkpsf'):
        print 'PSF does not exist'
        write_no_mkpsf(family_name, '{} {}'.format(expnum.strip('p'), header[_CCD].split('d')[1]))
        return

    '''
    target = header[_TARGET]
    if 'BH' not in target:
    print '-- NOT in H Block'
    return
    '''

    # reject any object too bright that will definetly be saturated
    mag = asteroid_id['mag'].values[0]
    if mag < 18.5:
        print '>> Object is too bright for accurate photometry'
        write_too_bright(family_name, asteroid_id)
        return

    # get fwhm from OSSOS VOSpace file
    fwhm = storage.get_fwhm(expnum.strip('p'), header[_CCD].split('d')[1])

    ast_psf = get_asteroid_data(asteroid_id, data, fwhm, family_name)
    star_psf = build_star_psf(ast_psf, expnum, header, asteroid_id, data, fwhm)

    # comet_psf = build_comet_psf(star_psf)

    print '-- Comparing PSFs'
    detection, sig = compare_psf(star_psf, ast_psf)
    if detection:
        print '>> Detect possible comae <<'
        write_to_file(asteroid_id, family_name, sig)


def fits_data(object_name, expnum, family_name):
    """
    Creates local copy of fits file from VOSpace
    """

    if family_name == 'none':
        vos_dir = '{}/none'.format(_VOS_DIR)
    else:
        vos_dir = '{}/all'.format(_VOS_DIR)

    assert storage.exists(vos_dir)
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


def get_asteroid_data(object_data, data, fwhm, family_name):
    """
    Calculate psf of asteroid, taking into acount trailing effect
    Cut a square around the polygon of the object, remove saturated rows, rotate so that the ellipse is parallel to the
    horizontal axis, then mask shape of ellipse around the object
    """

    # define parameters of the square to cut out from the polygon around the elongated object
    x_min = int(object_data['xmin'].values[0] - _BUFFER2)
    x_max = int(object_data['xmax'].values[0] + _BUFFER2)
    y_min = int(object_data['ymin'].values[0] - _BUFFER2)
    y_max = int(object_data['ymax'].values[0] + _BUFFER2)

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
    data_rot = rotate(data_obj, math.degrees(object_data['theta'].values[0]))

    # cut out a circular aperture around the object
    a = 0.5 * math.sqrt((object_data['xmax'].values[0] - object_data['xmin'].values[0]) ** 2 +
                        (object_data['ymax'].values[0] - object_data['ymin'].values[0]) ** 2)
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
    '''
    hdu = fits.PrimaryHDU()
    hdu.data = np.multiply(data_rot, mask2)
    hdu.writeto('psf_output/{}/cutout_{}_{}.fits'.format(family_name, object_data['object'].values[0],
                                                         object_data['expnum'].values[0]), clobber=True)
    '''

    # take the mean of all the values in each ROW
    data_mean = np.ma.mean(data_cutout, axis=1)
    data_nonzero = data_mean[np.nonzero(np.ma.fix_invalid(data_mean, fill_value=0))[0]]
    assert len(data_nonzero) > 0, "All data in cutout is zero"

    return data_nonzero


def build_star_psf(data_ast, expnum, header, asteroid_id, exp_data, fwhm):
    """
    Center the star data on the peak of the asteroid psf, and interpolate data points
    """

    bkg, flux, fluxerr = sep_phot(exp_data, asteroid_id)

    mag_max = -2.5 * np.log10(flux + fluxerr) + header[_ZMAG]
    mag_min = -2.5 * np.log10(flux - fluxerr) + header[_ZMAG]
    step = (mag_max - mag_min) / 2

    # calculate mean psf
    uri = storage.get_uri(expnum.strip('p'), header[_CCD].split('d')[1])
    ossos_psf = '{}.psf.fits'.format(uri.strip('.fits'))
    local_psf = '{}{}.psf.fits'.format(expnum, header[_CCD].split('d')[1])
    local_file_path = '{}/{}'.format(_STAMPS_DIR, local_psf)
    storage.copy(ossos_psf, local_file_path)

    ratio = 0
    while ratio < 0.97:
        data_str, ratio = find_best_mag(asteroid_id, bkg, mag_min, data_ast, fwhm, local_file_path, local_psf)
        mag_min += step
        os.unlink(local_psf)

    os.unlink(local_file_path)

    return data_str


def find_best_mag(asteroid_id, bkg, mag, data_ast, fwhm, local_file_path, local_psf):
    """
    Iterate through magnitude values to find best amplitude fit for star psf proflie
    """

    data_str = get_star_data(asteroid_id, bkg, mag, local_file_path, local_psf)

    x_ast = range(len(data_ast))
    x_str = range(len(data_str))
    p_str = [np.amax(data_str), float(len(x_str)) / 2, fwhm, 0.]
    p_ast = [np.amax(data_ast), float(len(x_ast)) / 2., fwhm, 0.]

    # fit a gaussian to the asteroid psf and then the star psf, will shift from fit peak x-value
    fitp_ast, fitco_ast = curve_fit(gauss, x_ast, data_ast, p_ast)
    fitp_str, fitco_str = curve_fit(gauss, x_str, data_str, p_str)

    # shift the center of the star psf to the center of the asteroid psf
    gauss_shift = fitp_str[1] - fitp_ast[1]
    x_str_shift = np.subtract(x_str, gauss_shift)

    # interpolate the psf values of the star at the x values that there is data for the asteroid psf
    y_str = lambda x: np.interp(x, x_str_shift, data_str)
    data_str_at_astpts = y_str(x_ast)

    ast_baseline = np.mean(np.sort(data_ast)[:6])
    str_baseline = np.mean(np.sort(data_str_at_astpts)[:6])
    data_str_at_astpts += (ast_baseline - str_baseline)

    ast_amp = np.amax(data_ast) - ast_baseline
    str_amp = np.amax(data_str) - str_baseline
    ratio = str_amp / ast_amp
    print ratio, mag

    return data_str_at_astpts, ratio


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
    kronrad, krflag = sep.kron_radius(data2, asteroid_id['x_mid'], asteroid_id['y_mid'], asteroid_id['a'],
                                      asteroid_id['b'], asteroid_id['theta'], ap)
    flux, fluxerr, flag = sep.sum_ellipse(data2, asteroid_id['x'], asteroid_id['y'], asteroid_id['a'], asteroid_id['b'],
                                          asteroid_id['theta'], 2.5 * kronrad, subpix=1, err=bkg.globalrms)

    return bkg.globalback, flux, fluxerr


def get_star_data(asteroid_id, bkg, mag, local_file_path, local_psf):
    """
    From ossos psf fitted image, calculate line profile
    """

    # pvwcs = wcs.WCS(header)
    # x, y = pvwcs.sky2xy(asteroid_id['ra'].values, asteroid_id['dec'].values)
    x = asteroid_id['x_mid'].values[0]
    y = asteroid_id['y_mid'].values[0]



    # run seepsf on the mean psf image
    iraf.set(uparm="./")
    iraf.digiphot(_doprint=0)
    iraf.apphot(_doprint=0)
    iraf.daophot(_doprint=0)
    iraf.seepsf(local_file_path, local_psf, xpsf=x, ypsf=y, magnitude=mag)

    with fits.open(local_psf) as hdulist:
        data = hdulist[0].data

    th = math.degrees(asteroid_id['theta'].values[0])
    data_rot = rotate(data, th)
    data_rot = np.ma.masked_where(data_rot == 0, data_rot)

    # add in background value of exposure instead of subtracting it from the asteroid PSF
    data_rot += bkg

    data_mean = np.ma.mean(data_rot, axis=1)

    return data_mean[np.nonzero(np.ma.fix_invalid(data_mean, fill_value=0))[0]]


def gauss(x, *p):
    amp, mu, sigma, b = p
    return amp * np.exp(-(x - mu) ** 2 / (2. * sigma ** 2)) + b


def build_comet_psf(data_str):
    x_str = range(len(data_str))

    multiple = []
    step = 0.25
    for x in range(0, int(len(x_str) / 4)):
        multiple.append(1 + step * x)

    mult_arr = []
    for item in multiple:
        mult_arr.append(item)
    for item in multiple[::-1]:
        mult_arr.append(item)
    for item in multiple:
        mult_arr.append(item)
    for item in multiple[::-1]:
        mult_arr.append(item)

    if len(mult_arr) != len(data_str):
        mult_arr.append(1)

    amp = np.amax(data_str) - np.amin(data_str)

    return np.add(data_str, np.multiply(amp / 2, mult_arr))


def compare_psf(star_psf, ast_psf):
    """
    Compare psf of asteroid against mean of stars, check if anomaly in wings
    >> changed data to not normalized nor baseline subtracted
    """

    x_ast = range(len(ast_psf))

    two_sample = ttest_ind(ast_psf, star_psf)
    print two_sample

    print ">> Ratio test ast/str"
    star_psf_sig = np.sqrt(np.absolute(star_psf))
    ast_psf_sig = np.sqrt(np.absolute(ast_psf))
    r = np.divide(ast_psf, star_psf)
    # r_sig = r * (np.divide(ast_psf_sig, ast_psf) + np.divide(star_psf_sig, star_psf))
    r_sig2 = r * np.sqrt((np.divide(ast_psf_sig, ast_psf) ** 2 + np.divide(star_psf_sig, star_psf) ** 2))

    print r
    print r_sig2
    print np.ma.mean(r), np.ma.mean(np.sort(r)[2:-3])
    print '>> (r - r_mean) / r_sig'
    ratio = (r - np.ma.mean(np.sort(r)[2:-3])) / r_sig2
    print ratio

    with sns.axes_style('ticks'):
        # latexify()
        plt.plot(x_ast, ast_psf, label='Asteroid', ls='-')
        plt.plot(x_ast, star_psf, label='Stellar Model', ls='--')
        plt.xlabel('Pixels')
        plt.ylabel('Mean Flux')
        plt.legend()
        plt.show()

    with sns.axes_style('ticks'):
        # latexify()
        plt.plot(x_ast, r, label='r', ls='-')
        plt.plot(x_ast, np.ones(len(x_ast)) * np.ma.mean(np.sort(r)[2:-3]) + r_sig2, label='rmean+rsig2 Model', ls='--')
        plt.plot(x_ast, np.ones(len(x_ast)) * np.ma.mean(np.sort(r)[2:-3]), label='rmean', ls=':')
        plt.legend()
        plt.show()

    if len(ratio[np.greater_equal(ratio, 5)]) > 2:
        return True, 5
    elif len(ratio[np.greater_equal(ratio, 4)]) > 2:
        return True, 4
    elif len(ratio[np.greater_equal(ratio, 3)]) > 2:
        return True, 3
    else:
        return False, 0


def write_to_file(asteroid_id, family_name, sig):
    with open('{}/{}/abv_{}_sig.txt'.format(_OUTPUT_DIR, family_name, sig), 'a') as outfile:
        outfile.write('{} {}\n'.format(asteroid_id['object'].values[0], asteroid_id['expnum'].values[0]))


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
    asteroid_expnum = '{} {}'.format(asteroid_id['object'].values[0], asteroid_id['expnum'].values[0])
    with open('{}/{}/{}'.format(_OUTPUT_DIR, family_name, _OUTPUT_TOO_BRIGHT), 'r') as outfile:
        for line in outfile:
            object_expnum.append('{} {}'.format(line.split(' ')[0], line.split(' ')[0]))
    if asteroid_expnum not in object_expnum:
        with open('{}/{}/{}'.format(_OUTPUT_DIR, family_name, _OUTPUT_TOO_BRIGHT), 'a') as outfile:
            outfile.write(
                '{} {} {}\n'.format(asteroid_id['object'].values[0], asteroid_id['expnum'].values[0],
                                    asteroid_id['mag'].values[0]))


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


if __name__ == '__main__':
    main()