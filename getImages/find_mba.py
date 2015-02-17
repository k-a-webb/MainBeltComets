import numpy as np
import requests
import os

def main():
    
    find_mbas()
    
def find_mbas(output=None):    
    '''
    Queries the AstDys database for members of the main asteroid belt
    2.064 AU < a < 3.277 AU, e < 0.45, i < 40 deg
    '''
    
    if output == None:
        output = 'mba_list.txt'
        
    output_dir = '/Users/admin/Desktop/MainBeltComets/getImages/mba'
    
    BASEURL = 'http://hamilton.dm.unipi.it/~astdys2/propsynth/numb.syn'
    
    r = requests.get(BASEURL)
    r.raise_for_status()
    
    table = r.content
    table_lines = table.split('\n')
        
    print "----- Searching for obects of in the Main Asteroid Belt -----"
        
    mba_list = []
    mba_wofam_list = []
    
    with open('mba/mba_wofam_list.txt') as infile:
        for line in infile:
            mba_wofam_list.append(line)
            
        
    for line in table_lines[2:]:
        if len(line.split()) > 0:
            mba_name = line.split()[0]
            mba_a = float(line.split()[2])
            mba_e = float(line.split()[3])
        
            for item in mba_wofam_list:
                if (item == mba_name) & (mba_a > 2.064) & (mba_e < 0.45):
                    mba_list.append(mba_name)

    #assert len(mba_list) > 0
    print " Number of objects in the M.A.B.: {}".format(len(mba_list))            
    
    with open('{}/{}'.format(output_dir, output), 'w') as outfile:
        for item in mba_list:
              outfile.write("{}\n".format(item))
              
    return mba_list


if __name__ == '__main__':
    main()