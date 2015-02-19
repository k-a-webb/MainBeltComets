import os
from collections import Counter
import requests
from astropy.table import Table
import pandas as pd
import numpy as np

print '---------- \nSearching for orbital data for objects with family designations \n----------'


family_list = []
object_list = []
objects = []

with open('asteroid_families/families_with_images.txt') as infile:
    for line in infile:
        family_list.append(line.strip('\n'))
        
for familyname in family_list:
    with open('asteroid_families/{}/{}_images.txt'.format(familyname, familyname), 'r+') as infile:
        next(infile)
        for line in infile:
            objects.append(line.split(' ')[0])
            #object_list.append(line.strip('\n'))
            
'''
with open('asteroid_families/all_images.txt', 'w') as outfile:
    for item in object_list:
          outfile.write("{}\n".format(item))
'''
counts = Counter(objects)

new_objects = []
for item in counts:
    new_objects.append(item)
    
num_images = sum(counts.values())
values =  counts.values()

#print num_images
#print counts.most_common(10)

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

#print counts.values(new_objects[2])
#print objects.count(new_objects[2])

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

table2_arrays = {'objectname': name_list, 'num_of_observations': occurance, 'semimajor_axis': a_list2, 'eccentricity': e_list2, 'sin_inclination': sini_list2}
objects_in_fam_table = pd.DataFrame(data=table2_arrays)

#print objects_in_fam_table

objects_in_fam_table.to_csv('asteroid_families/mba_fam_data.txt', sep='\t', encoding='utf-8', index=False)

#with open('asteroid_families/mba_fam_data.txt', 'w') as outfile:
#    outfile.write('{}'.format(objects_in_fam_table))



                
        
        
        
        
        
        
        
        
        
        
        
        