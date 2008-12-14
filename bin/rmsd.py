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
    crds = numpy.zeros((len(atomlist), 3), float)
    for i, a in enumerate(atomlist):
        crds[i,0] = a[0]
        crds[i,1] = a[1]
        crds[i,2] = a[2]
    return crds

def get_crds_mirror(atomlist):
    crds = numpy.zeros((len(atomlist), 3), float)
    for i, a in enumerate(atomlist):
        crds[i,0] = -a[0]
        crds[i,1] = a[1]
        crds[i,2] = a[2]
    return crds

import os

def get_structlist(outputfile):
  structlist = []
  cmd = 'bzcat '+outputfile+' | xmlstarlet sel -t -m \'/OutputData/StructureImages\' -m \'Image\' -m \'Atom\' -v \'@x\' -o "," -v \'@y\' -o \',\' -v \'@z\' -o \':\' -b -n | gawk \'{if (NF) print $0}\''

  for line in os.popen(cmd).readlines():
    line = line.rstrip(':')
    line = line.rstrip(':\n')
    
    atomlist = []
    for atom in line.split(':'):
      atomcoords = []
      
      for x in atom.split(','):
        atomcoords.append(float(x))
        
      atomlist.append(atomcoords)
        
    structlist.append(atomlist)

  print "Number of structures is "+len(structlist)
  if (len(structlist) > 100):
    print "\nTruncating to 200"
    del structlist[200:len(structlist)]

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
      sum += min(rmsd(crds1, get_crds_mirror(atomlist2)), rmsd(crds1, get_crds(atomlist2)))
    
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

for data in filedata:
  print data[0], data[2]
