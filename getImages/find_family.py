import numpy as np
import requests
import argparse
import os
import pandas as pd

_DIR_PATH_BASE = os.path.dirname(os.path.abspath(__file__))
_OUTPUT_DIR = '{}/family_lists'.format(_DIR_PATH_BASE)
if not os.path.isdir(_OUTPUT_DIR):
    os.makedirs(_OUTPUT_DIR)

_BASEURL_NUMB_FAM = 'http://hamilton.dm.unipi.it/~astdys2/propsynth/numb.famtab'
_BASE_URL_NUMB_MEM = 'http://hamilton.dm.unipi.it/~astdys2/propsynth/numb.members'
_BASE_URL_FAMREC = 'http://hamilton.dm.unipi.it/~astdys2/propsynth/numb.famrec'
_BASE_URL_SYN = 'http://hamilton.dm.unipi.it/~astdys2/propsynth/numb.syn'

"""
Creates a list of members in a designated family - 'all' for all objects, or a specific family eg '3330' Lixiaohua
Also parses astdys database for some orbital parameters

"""


def main():
    # From the given input, make list of MBCs to query

    parser = argparse.ArgumentParser(
        description='Queries the AstDys database for members of specified family, family name is name of largest member.')
    parser.add_argument("--family", '-f',
                        action="store",
                        default="all",
                        help='Whether in designated family (all) or not (none)')
    parser.add_argument("--output", "-o",
                        action="store",
                        default="family.txt",
                        help='Location and name of output file containing object names')
    parser.add_argument("--status",
                        action='store',
                        default=4,
                        help='Family status, options: 0 not in family, 1 inner MB, 2 MB, 3 outter MB')
    parser.add_argument('--fromfile',
                        action='store',
                        default=True,
                        help='Whether to find family members from search of from input file')

    args = parser.parse_args()

    # find_family_members(args.family, args.output)
    mba_list = find_by_status(args.status)
    name_list = parse_for_all(mba_list, args.status, fromfile=True)


def find_family_members(familyname, output=None):
    """
    Queries the AstDys database for members of specified family, family name is name of largest member
    """

    if output is None:
        output = '{}/{}_family.txt'.format(_OUTPUT_DIR, familyname)
    else:
        output = '{}/{}_family.txt'.format(_OUTPUT_DIR, familyname)

    r = requests.get(_BASE_URL_NUMB_MEM)
    r.raise_for_status()
    table = r.content
    table_lines = table.split('\n')

    print "----- Searching for obects of in family {} -----".format(familyname)

    asteroid_list = []
    for line in table_lines[1:len(table_lines) - 1]:
        assert len(line.split()) > 0
        familyname_infile = line.split()[3]
        if familyname_infile == '{}'.format(familyname):
            asteroid_list.append(line.split()[0])

    assert len(asteroid_list) > 0
    print " Number of members in family {}: {}".format(familyname, len(asteroid_list))

    with open(output, 'w') as outfile:
        for item in asteroid_list:
            outfile.write("{}\n".format(item))

    return asteroid_list


def get_all_families_list(output=None):
    """
    Queries the AstDys database for members of all known families, family name is name of largest member
    """

    if output is None:
        output = '{}/all_family_new.txt'.format(_OUTPUT_DIR)

    r = requests.get(_BASEURL_NUMB_FAM)
    r.raise_for_status()
    table = r.content
    table_lines = table.split('\n')

    print "----- Searching for list of all family names -----"

    familyname_list = []
    for line in table_lines[2:len(table_lines) - 1]:  # skip header lines
        assert len(line.split()) > 0
        familyname_list.append(line.split()[0])

    assert len(familyname_list) > 0
    print " Number of asteroid families: {}".format(len(familyname_list))

    with open(output, 'w') as outfile:
        for item in familyname_list:
            outfile.write("{}\n".format(item))

    return familyname_list


def find_by_status(status=4):
    """
    Queries the AstDys database for members of the main asteroid belt
    with a particular family designation status [0, 1, 2, 3]

    At the moment it takes:
        0 as not in family
        1 in inner main belt
        2 in middle main belt
        3 outer main belt
        4 either 1,2,3
    i.e. it searches either for members of families (3) or not (0)
    """

    output = '{}/status_{}.txt'.format(_OUTPUT_DIR, status)

    # create duplicate output to satisfy naming conventions
    if status == '0':
        output2 = '{}/none_family.txt'.format(_OUTPUT_DIR)
    else:
        output2 = '{}/all_family.txt'.format(_OUTPUT_DIR)

    r = requests.get(_BASE_URL_FAMREC)
    r.raise_for_status()
    table = r.content
    table_lines = table.split('\n')

    print "----- Searching for objects of in the Main Asteroid Belt family designation status {} -----".format(status)

    object_list = []
    for line in table_lines[1:]:
        if len(line.split()) > 0:
            table_name = line.split()[0]
            table_status = line.split()[2]
            if status == 4:
                if str(status) != '0':
                    object_list.append(table_name)
            else:
                if str(status) == str(table_status):
                    object_list.append(table_name)

    assert len(object_list) > 0
    print " Number of objects in the M.A.B. family designation status {}: {}".format(status, len(object_list))

    with open(output, 'w') as outfile:
        for item in object_list:
            outfile.write("{}\n".format(item))
    with open(output2, 'w') as outfile:
        for item in object_list:
            outfile.write("{}\n".format(item))

    return object_list


def get_all_astrometry():
    """
    Parse through the orbital information from AstDys for semimajor axis, eccentricity, inclination, write to new table
    I hope I don't actually use this anywhere but I'm afraid to delete it
    """

    print '----- Getting astrometry -----'

    r = requests.get(_BASE_URL_SYN)
    r.raise_for_status()
    table = r.content
    table_lines = table.split('\n')

    tableobject_list = []
    a_list = []
    e_list = []
    sini_list = []
    for line in table_lines[2:]:
        if len(line.split()) > 0:
            tableobject_list.append(line.split()[0])
            a_list.append(float(line.split()[2]))
            e_list.append(float(line.split()[3]))
            sini_list.append(float(line.split()[4]))
    i_list = np.degrees(np.arcsin(sini_list))
    table_arrays = {'objectname': tableobject_list, 'semimajor_axis': a_list, 'eccentricity': e_list,
                    'inclination': i_list}
    all_objects_table = pd.DataFrame(data=table_arrays)
    all_objects_table.to_csv('{}/astdys_table.txt'.format(_OUTPUT_DIR), sep='\t', encoding='utf-8')
    # print all_objects_table

    return all_objects_table


def parse_for_all(mba_list, status, fromfile=True):
    """
    Queries the AstDys database for all objects of a specified status and gets astrometry
    """

    if fromfile is True:
        all_objects_table = pd.read_table('{}/astdys_table.txt'.format(_OUTPUT_DIR), sep='\t', dtype={'objectname': object})
    else:
        all_objects_table = get_all_astrometry()

    print "----- Searching for obects of in AstDys -----"

    a_list2 = []
    e_list2 = []
    i_list2 = []
    name_list = []

    for objectname in mba_list:
        index = all_objects_table.query('objectname == "{}"'.format(objectname))
        try:
            name_list.append(objectname)
            a_list2.append(float(index['semimajor_axis']))
            e_list2.append(float(index['eccentricity']))
            i_list2.append(float(index['inclination']))
        except:
            print objectname
            print index
    print " Number of objects found: {}".format(len(name_list))

    table2_arrays = {'objectname': name_list, 'semimajor_axis': a_list2, 'eccentricity': e_list2,
                     'inclination': i_list2}
    objects_table = pd.DataFrame(data=table2_arrays)

    objects_table.to_csv('{}/all_data_status_{}.txt'.format(_OUTPUT_DIR, status), sep='\t', encoding='utf-8',
                         index=False)

    '''
    # optional, same output as find_by_status()
    with open('{}/all_list_status_{}.txt'.format(_OUTPUT_DIR, status), 'w') as outfile:
        for item in name_list:
            outfile.write('{}\n'.format(item))
    '''

    return name_list


def parse_for_mba(mba_list, status=3):
    """
    Queries the AstDys database for members of the main asteroid belt and astrometry
    a > 2.067 AU, e < 0.45, i < 40 deg, sini < 0.642

    WARNING: this should not be used, we don't want to exclude objects that we have observations of

    """

    try:
        all_objects_table = pd.read_table('{}/astdys_table.txt'.format(_OUTPUT_DIR), sep='\t')
    except:
        all_objects_table = get_all_astrometry()

    print "----- Searching for obects of in the Main Asteroid Belt with a > 2.064 AU, e < 0.45, i < 40 deg -----"

    a_list2 = []
    e_list2 = []
    sini_list2 = []
    name_list = []
    x = 0
    for objectname in mba_list:
        index = all_objects_table.query('objectname == "{}"'.format(objectname))
        try:
            if (float(index['semimajor_axis']) > 2.064) & (float(index['eccentricity']) < 0.45) & (
                        float(index['sin_inclination']) < 0.642):
                name_list.append(objectname)
                a_list2.append(float(index['semimajor_axis']))
                e_list2.append(float(index['eccentricity']))
                sini_list2.append(float(index['sin_inclination']))
        except:
            x += 1
    i_list2 = np.degrees(np.arcsin(sini_list2))
    table2_arrays = {'objectname': name_list, 'semimajor_axis': a_list2, 'eccentricity': e_list2,
                     'inclination': i_list2}
    objects_in_mb_table = pd.DataFrame(data=table2_arrays)

    print objects_in_mb_table

    print " Number of objects found: {}".format(len(name_list))

    objects_in_mb_table.to_csv('{}/main_belt_status_{}.txt'.format(_OUTPUT_DIR, status), sep='\t', encoding='utf-8',
                               index=False)

    return name_list


if __name__ == '__main__':
    main()