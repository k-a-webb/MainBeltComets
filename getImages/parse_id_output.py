import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import numpy as np

_DIR_PATH_BASE = os.path.dirname(os.path.abspath(__file__))
_DIR_PATH = '{}/asteroid_families'.format(_DIR_PATH_BASE)
_INPUT = '{}/asteroid_families/all/phot_output/test_output.txt'.format(_DIR_PATH_BASE)


def main():
    parser = argparse.ArgumentParser(
        description='For a given set of fits images: preforms photometry, identifies a specific object, and returns \
                        the orbital elements of the object')
    parser.add_argument("--object", '-obj',
                        action="store",
                        default='74024',
                        help="Asteroid name.")

    args = parser.parse_args()

    light_curve(args.object)


def light_curve(objectname):
    """
    Produce a graph of magnitude as a function of time
    """

    input_table = pd.read_table(_INPUT, sep='\t', usecols=[0, 4, 5, 9, 10, 11], dtype={'object': object})

    object_table = input_table.query('object == "{}"'.format(objectname))
    object_table.drop_duplicates()
    print object_table.sort('time')

    time_start = np.amin(object_table['time'])
    time_stop = np.amax(object_table['time'])

    with sns.axes_style("ticks"):
        pg = sns.PairGrid(object_table, x_vars=["time"], y_vars=['mag'], aspect=4)
        # plt.xlim(time_start, time_stop)
        # plt.ylim(13, 22)
        pg.map(plt.scatter)
        plt.show()

def parse_ids():
    """
    Select IDs with good number of background stars and satisfy both magnitude and ellipticity conditions
    """

    input_table = pd.read_table(_INPUT, sep='\t', usecols=[0, 1, 4, 5, 8, 9, 10, 11], dtype={'object': object})

    good_objects = input_table.query('(stars > 30) & (consistent_f == "yes") & (consistent_mag == "yes")')
    print good_objects

def remove_duplicates():

    print '-- Removing duplicates'

    infile = 'no_object_found.txt'
    input = '{}/all/phot_output/{}'.format(_DIR_PATH, infile)

    lines = []
    with open(input, 'rw+') as infile:
        for line in infile:
            if line not in lines:
                lines.append(line)
            else:
                del line

if __name__ == '__main__':
    main()