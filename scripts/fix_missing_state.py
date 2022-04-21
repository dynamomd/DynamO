#!/usr/bin/env python3

import os,pickle

root='HSTetherNVTPhiTWD'

def tryconvert(a):
    try:
        return int(a)
    except ValueError:
        try:
            return float(a)
        except ValueError:
            return a

for subdir, dirs, files in os.walk(root):
    #skip the first entry, its just the root folder
    if subdir == root:
        continue
    
    if 'state.pkl' not in files:
        #Try to regenerate the state
        tags = subdir.split('_')
        #The first entry contains the subdir, remove it
        tags[0] = tags[0].split('/')[1]

        #The last element is the run, skip it
        tags = tags[:-1]

        #Check if the strings can be turned into int or float
        tags = [tryconvert(a) for a in tags]
        
        tags = tuple((a,b) for a,b in zip(tags[0::2], tags[1::2]))

        with open(subdir+'/state.pkl', 'wb') as f:
            pickle.dump(tags, f)
        
        print(subdir, tags)
        
