import os
import requests
from astropy.table import Table
import pandas as pd
import numpy as np
from collections import Counter
import argparse
import matplotlib.pyplot as plt
import seaborn as sns


def main(): 
    
    parser = argparse.ArgumentParser(
                        description='Get metadata')
    parser.add_argument("--family", '-f',
                        action="store",
                        default='all',
                        help="Asteroid family name. Usually the asteroid number of the largest member.")
                            
    args = parser.parse_args()
    
    parse_metadata(args.family)

def get_metadata(familyname):

    print '----- Searching for orbital data for objects with family designations -----'
    
    if os.path.exists('asteroid_families/{}/{}_images.txt'.format(familyname, familyname)) == False:
        objects, all_images = get_all_image_list(familyname)
    else:
        objects = []
        with open('asteroid_families/{}/{}_images.txt'.format(familyname, familyname)) as infile:
            next(infile)
            for line in infile:
                objects.append(line.split(' ')[0])
    
    counts = Counter(objects)
    new_objects = []
    for item in counts:
        new_objects.append(item)        

    BASEURL = 'http://hamilton.dm.unipi.it/~astdys2/propsynth/numb.syn'
    r = requests.get(BASEURL)
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
    table_arrays = {'objectname': tableobject_list, 'semimajor_axis': a_list, 'eccentricity': e_list, 'sin_inclination': sini_list}
    all_objects_table = pd.DataFrame(data=table_arrays)
    #print all_objects_table

    a_list2 = []
    e_list2 = []
    sini_list2 = []
    name_list = []
    occurance = []

    for objectname in new_objects[1:]:
        occur = objects.count(objectname)
        try:
            index = all_objects_table.query('objectname == "{}"'.format(objectname))
            name_list.append(objectname)
            a_list2.append(float(index['semimajor_axis']))
            e_list2.append(float(index['eccentricity']))
            sini_list2.append(float(index['sin_inclination']))
            occurance.append(occur)
        except:
            print 'Could not find object {}'.format(objectname)

    i_list2 = np.degrees(np.arcsin(sini_list2))

    objects_in_fam_table = pd.DataFrame({'objectname': name_list, 
                                            'occurance': occurance, 
                                            'semimajor_axis': a_list2, 
                                            'eccentricity': e_list2, 
                                            'inclination': i_list2})
    objects_in_fam_table.to_csv('asteroid_families/{}/{}_metadata.txt'.format(familyname, familyname), sep='\t', encoding='utf-8', index=False)

    return objects_in_fam_table

def get_all_image_list(familyname):
    
    family_list = []
    objects = []
    all_images = []
    
    with open('asteroid_families/families_with_images.txt') as infile:
        for line in infile:
            family_list.append(line.strip('\n'))
    
    for familyname in family_list:
        with open('asteroid_families/{}/{}_images.txt'.format(familyname, familyname), 'r+') as infile:
            next(infile)
            for line in infile:
                objects.append(line.split(' ')[0])
                all_images.append(line)

    with open('asteroid_families/all_images.txt', 'w') as outfile:
        for item in all_images:
              outfile.write("{}".format(item))
              
    return objects

def parse_metadata(familyname):
    
    if os.path.exists('asteroid_families/{}/{}_metadata.txt'.format(familyname, familyname)) == True:
        data_table = pd.read_table('asteroid_families/{}/{}_metadata.txt'.format(familyname, familyname))
    else:
        data_table = get_metadata(familyname)
      
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
        
        
        
        
        
        
        
        
        
        
        
        