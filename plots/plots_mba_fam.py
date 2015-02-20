import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

#input_table = pd.DataFrame.from_csv('asteroid_families/mba_fam_data.csv')

name = []
occurance = []
a = []
e = []
sini = []

with open('../getImages/asteroid_families/mba_fam_data.txt') as infile:
    next(infile) # skip header 
    for line in infile:
        name.append(line.split()[2])
        occurance.append(float(line.split()[1]))
        a.append(float(line.split()[3]))
        e.append(float(line.split()[0]))
        sini.append(float(line.split()[4]))

i_rad = np.arcsin(sini)
i = np.degrees(i_rad)

input_arrays = {'name': name, 'occurance': occurance, 'semimajor_axis': a, 'eccentricity': e, 'inclination': i}
input_table = pd.DataFrame(data=input_arrays)
#print input_table

pg = sns.PairGrid(input_table, x_vars=["semimajor_axis"], y_vars=["eccentricity", "inclination"], 
                  aspect=4, hue="occurance", palette="coolwarm")
plt.xlim(2.15, 4)
pg.axes[0][0].set_ylim([0,0.35])
plt.ylim(0,40)  
pg.map(plt.scatter, alpha=0.7) 
pg.add_legend(title='occurance')
plt.show()

pg = sns.PairGrid(input_table, x_vars=["semimajor_axis"], y_vars=['occurance'], aspect=4)
plt.xlim(2.15, 4)
plt.ylim(0,35)  
pg.map(plt.scatter) 
plt.show()

pg = sns.PairGrid(input_table, x_vars=["inclination"], y_vars=['occurance'], aspect=4)
plt.xlim(0, 35)
plt.ylim(0, 35)  
pg.map(plt.scatter) 
plt.show()

pg = sns.PairGrid(input_table, x_vars=["eccentricity"], y_vars=['occurance'], aspect=4)
plt.xlim(0, 0.3)
plt.ylim(0, 35)  
pg.map(plt.scatter) 
plt.show()