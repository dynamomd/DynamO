#!/usr/bin/env python

import math
import numpy

def rmsd(crds1, crds2):
  """Returns RMSD between 2 sets of [nx3] numpy array"""
  assert(crds1.shape[1] == 3)
  assert(crds1.shape == crds2.shape)
  n_vec = numpy.shape(crds1)[0]
  correlation_matrix = numpy.dot(numpy.transpose(crds1), crds2)
  v, s, w_tr = numpy.linalg.svd(correlation_matrix)
  is_reflection = (numpy.linalg.det(v) * numpy.linalg.det(w_tr)) < 0.0
  if is_reflection:
    s[-1] = - s[-1]
  E0 = sum(sum(crds1 * crds1)) + \
       sum(sum(crds2 * crds2))
  rmsd_sq = (E0 - 2.0*sum(s)) / float(n_vec)
  rmsd_sq = max([rmsd_sq, 0.0])
  return numpy.sqrt(rmsd_sq)


def get_crds(atomlist):
  return numpy.array(atomlist)

def mirror_crds(atomlist):
  crds = atomlist

  for i in range(len(atomlist)):
    crds[i,0] = -crds[i,0]

  return crds
  
import os

def get_structlist(outputfile):
  cmd = 'bzcat '+outputfile+' | xmlstarlet sel -t -m \'/OutputData/StructureImages\' -m \'Image\' -m \'Atom\' -v \'@x\' -o "," -v \'@y\' -o \',\' -v \'@z\' -o \':\' -b -n | gawk \'{if (NF) print $0}\''

  structlist = []
  for line in os.popen(cmd).readlines():
    line = line.rstrip(':')
    line = line.rstrip(':\n')
    
    atomlist = [[float(x) for x in atom.split(',')] for atom in line.split(':')]
        
    structlist.append(atomlist)

  print "Number of structures is "+str(len(structlist))
  if (len(structlist) > 100):
    print "\nTruncating to 100"
    del structlist[100:len(structlist)]

  return structlist


def get_best_crds(structlist):
  
  #initialise the best value using the first entry
  bestcrds = get_crds(structlist[0])
  minsum = 0

  for atomlist2 in structlist:
    crds2 = get_crds(atomlist2)
    minsum += rmsd(bestcrds, crds2)

  #Now test if any have a lower sum
  for atomlist1 in structlist:
    crds1 = get_crds(atomlist1)
    sum = 0

    for atomlist2 in structlist:
      crds2 = get_crds(atomlist2)
      sum += min(rmsd(crds1, mirror_crds(crds2)), rmsd(crds1, crds2))
    
    if (sum < minsum):
      minsum = sum
      bestcrds = crds1
  
  return bestcrds, minsum / len(structlist)

def get_temperature(file):
    cmd = 'bzcat '+file+' | xmlstarlet sel -t -v \'/OutputData/KEnergy/T/@val\''
    return float(os.popen(cmd).read())
  
import sys
 
filedata = []
for file in  sys.argv[1:]:
  print "Processing file "+file
  bestcrds, avgmsd = get_best_crds(get_structlist(file))
  filedata.append([get_temperature(file), bestcrds, avgmsd])

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
    print >>f, filedata[val][0], min(rmsd(filedata[val][1],filedata[val+1][1]),rmsd(filedata[val][1], mirror_crds(filedata[val+1][1])))
finally:
  f.close()
