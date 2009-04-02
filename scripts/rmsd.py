#!/usr/bin/env python

import math
import numpy
import sys
import os
import pickle

xmlstarlet='xml'

def rmsd(crds1, crds2):
  """Returns RMSD between 2 sets of [nx3] numpy array"""
  v, s, w_tr = numpy.linalg.svd(numpy.dot(numpy.transpose(crds1), crds2))
  return numpy.sqrt(max([(sum(sum(crds1 * crds1)) + sum(sum(crds2 * crds2)) - 2.0*sum(s)) / float(numpy.shape(crds1)[0]), 0.0]))

def get_structlist(outputfile):
  cmd = 'bzcat '+outputfile+' | '+xmlstarlet+' sel -t -m \'/OutputData/StructureImages\' -m \'Image\' -m \'Atom\' -v \'@x\' -o "," -v \'@y\' -o \',\' -v \'@z\' -o \':\' -b -n | gawk \'{if (NF) print $0}\''

  structlist = []
  for line in os.popen(cmd).readlines():
    line = line.rstrip(':')
    line = line.rstrip(':\n')
    structlist.append(numpy.array([[float(x) for x in atom.split(',')] for atom in line.split(':')]))

  print "Number of structures is "+str(len(structlist))
  if (len(structlist) > 100):
    print "\nTruncating to 100"
    del structlist[100:(len(structlist)-1)]

  return structlist


from operator import add

def get_best_crds(structlist):
  #initialise the best value using the first entry
  bestcrds = structlist[0]
  minsum = 1e308

  #Now test if any have a lower sum
  for atomlist1 in structlist:
    locsum = sum((min(rmsd(atomlist1, atomlist2), rmsd(atomlist1,atomlist2[::-1])) for atomlist2 in structlist))
    
    if (locsum < minsum):
      minsum = locsum
      bestcrds = atomlist1
  
  return bestcrds, minsum / len(structlist)

def get_temperature(file):
    cmd = 'bzcat '+file+' | '+xmlstarlet+' sel -t -v \'/OutputData/KEnergy/T/@val\''
    return float(os.popen(cmd).read())
  
filedata = []
for file in  sys.argv[1:]:
  print "Processing file "+file
  if (os.path.exists(file+'.rmspickle')):
    print "Found pickled data, skipping minimisation"
    f = open(file+'.rmspickle', 'r')
    filedata.append(pickle.load(f))
    f.close()
  else:
    print "Could not find pickled minimisation data"
    print "Processing file "+file
    bestcrds, avgmsd = get_best_crds(get_structlist(file))
    tempdata = [get_temperature(file), bestcrds, avgmsd]
    filedata.append(tempdata)
    f = open(file+'.rmspickle', 'w')
    pickle.dump(tempdata,f)
    f.close()
      
filedata.sort()
    


f = open('avgrmsd.dat', 'w')
try:
  for data in filedata:
    print >>f, data[0], data[2]
finally:
  f.close()


f = open('rmsddiff.dat', 'w')
try:
  for val in range(len(filedata)-1):
    print >>f, filedata[val][0], min(rmsd(filedata[val][1],filedata[val+1][1]), rmsd(filedata[val][1],filedata[val+1][1][::-1]))
finally:
  f.close()
