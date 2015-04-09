activate_this = '/Users/admin/Desktop/MainBeltComets/bin/activate_this.py'
execfile(activate_this, dict(__file__=activate_this))

import os
import requests
from astropy.table import Table
import pandas as pd
import numpy as np
from collections import Counter
import argparse
import matplotlib.pyplot as plt
import seaborn as sns

_DIR_PATH_BASE = os.path.dirname(os.path.abspath(__file__))
_FAMILY_LISTS = '{}/family_lists'.format(_DIR_PATH_BASE)
_STAMPS_DIR = '{}/postage_stamps'.format(_DIR_PATH_BASE)
_IMAGE_LISTS = '{}/image_lists'.format(_DIR_PATH_BASE)
_OUTPUT_DIR = '{}/phot_output'.format(_DIR_PATH_BASE)

_BASE_URL = 'http://hamilton.dm.unipi.it/~astdys2/propsynth/numb.syn'


def main():
    parser = argparse.ArgumentParser(description='Get metadata')
    parser.add_argument("--family", '-f',
                        action="store",
                        default='all',
                        help="Asteroid family name. Usually the asteroid number of the largest member.")

    args = parser.parse_args()

    # parse_metadata(args.family)
    get_metadata_for_images(args.family)
    # get_metadata_for_family(args.family)


def get_metadata_for_images(familyname):
    """
    Get orbital information from image list
    """

    print '----- Searching for orbital data for objects with family designations -----'

    try:
        image_table = pd.read_table('{}/{}_images.txt'.format(_IMAGE_LISTS, familyname), sep='\t',
                                    dtype={'Object': object})
        object_list = image_table['Object'].values
        print 'success'
    except Exception, e:
        print e
        return

    counts = Counter(object_list)
    object_list_no_repeats = []
    for item in counts:
        object_list_no_repeats.append(item)

    try:
        astdys_table = pd.read_table('{}/astdys_table.txt'.format(_FAMILY_LISTS), sep='\t',
                                     dtype={'objectname': object})
        print 'success'
    except:
        print '-- Querying AtsDys Database'
        r = requests.get(_BASE_URL)
        r.raise_for_status()

        tableobject_list = []
        a_list = []
        e_list = []
        sini_list = []

        table = r.content
        table_lines = table.split('\n')
        for line in table_lines[2:]:
            if len(line.split()) > 0:
                tableobject_list.append(line.split()[0])
                a_list.append(float(line.split()[2]))
                e_list.append(float(line.split()[3]))
                sini_list.append(float(line.split()[4]))
        i_list = np.degrees(np.arcsin(sini_list))
        table_arrays = {'objectname': tableobject_list, 'semimajor_axis': a_list, 'eccentricity': e_list,
                        'inclination': i_list}
        astdys_table = pd.DataFrame(data=table_arrays)
        # print all_objects_table

    a_list2 = []
    e_list2 = []
    i_list2 = []
    name_list = []
    occurance = []

    for objectname in object_list_no_repeats[1:]:
        occur = counts['{}'.format(objectname)]
        try:
            index = astdys_table.query('objectname == "{}"'.format(objectname))
            name_list.append(objectname)
            a_list2.append(float(index['semimajor_axis']))
            e_list2.append(float(index['eccentricity']))
            i_list2.append(float(index['inclination']))
            occurance.append(occur)
        except Exception, e:
            print 'Could not find object {}'.format(objectname)

    object_data = pd.DataFrame({'objectname': name_list,
                                'occurance': occurance,
                                'semimajor_axis': a_list2,
                                'eccentricity': e_list2,
                                'inclination': i_list2})
    object_data.to_csv('{}/{}_metadata.txt'.format(_IMAGE_LISTS, familyname), sep='\t', encoding='utf-8', index=False)

    return object_data


def get_metadata_for_family(family_name):
    """
    Searches AstDys for all orbital elements for every member of a given family
    """

    object_list = []
    with open('asteroid_families/{}/{}_family.txt'.format(family_name, family_name), 'r') as infile:
        for line in infile:
            object_list.append(line.strip('\n'))

    base_url = 'http://hamilton.dm.unipi.it/~astdys2/propsynth/numb.syn'
    r = requests.get(base_url)
    r.raise_for_status()

    tableobject_list = []
    a_list = []
    e_list = []
    sini_list = []

    table = r.content
    table_lines = table.split('\n')
    for line in table_lines[2:]:
        if len(line.split()) > 0:
            tableobject_list.append(line.split()[0])
            a_list.append(float(line.split()[2]))
            e_list.append(float(line.split()[3]))
            sini_list.append(float(line.split()[4]))
    i_list = np.degrees(np.arcsin(sini_list))
    astdys_table = pd.DataFrame({'object_name': tableobject_list,
                                 'semimajor_axis': a_list,
                                 'eccentricity': e_list,
                                 'inclination': i_list})

    index_list = []
    for object_name in object_list:
        index = astdys_table.query('object_name == "{}"'.format(object_name))
        index_list.append(index.index[0])

    family_table = astdys_table.iloc[index_list, :]
    family_table.to_csv('asteroid_families/{}/{}_all_metadata.txt'.format(family_name, family_name), sep='\t',
                        encoding='utf-8', index=False)


def parse_metadata(familyname):
    if os.path.exists('{}/{}_metadata.txt'.format(_IMAGE_LISTS, familyname)):
        data_table = pd.read_table('{}/{}_metadata.txt'.format(_IMAGE_LISTS, familyname))
    else:
        data_table = get_metadata_for_images(familyname)

    print data_table.sort('occurance')

    frequency = []
    for i in range(1, 34):
        occ = data_table.query('occurance == {}'.format(i))
        frequency.append(float(len(occ)))
    
    print frequency
    
    with sns.axes_style('ticks'):
        plt.hist(data_table['occurance'], 34, histtype="stepfilled", alpha=1)
        plt.ylabel('Number of asteroids')
        plt.xlabel('Observation occurance')
        plt.show()


if __name__ == '__main__':
    main()                
