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

def get_temperature(file):
    cmd = 'bzcat '+file+' | '+xmlstarlet+' sel -t -v \'/OutputData/KEnergy/T/@val\''
    return float(os.popen(cmd).read())  

from operator import add

def get_rmsd_data(structlist):
  arrayList = numpy.zeros( (len(structlist),len(structlist)), float)
  for i in range(len(structlist)):    
    for j in range(i+1,len(structlist)):
      arrayList[i][j] = min(rmsd(structlist[i], structlist[j]), rmsd(structlist[i], (structlist[j])[::-1]))
      arrayList[j][i] = arrayList[i][j]

  minimum = 0
  minsum = 1e308
  for i in range(len(structlist)):  
    localsum = sum(arrayList[i])
    if (localsum < minsum):
      minimum = i
      minsum = localsum

  return minimum, minsum / len(structlist), arrayList

import copy

def cluster_analysis(arrayList, length, threshold):
  localCluster = []
  for i in range(len(arrayList)):
    local = []
    count = 0
    for j in range(len(arrayList)):
      if (arrayList[i][j] < length):
        local.append(j)
        count += 1
    localCluster.append([count, local, i])
  
  clusters = []
  maxval = copy.deepcopy(max(localCluster))
  while (maxval[0] > threshold):
    clusters.append(maxval)
    for j in range(len(arrayList)):
      for id in maxval[1]:
        if (id in localCluster[j][1]):
          localCluster[j][1].remove(id)
          localCluster[j][0] -= 1
    maxval = copy.deepcopy(max(localCluster))
    
  return clusters

filedata = []
for file in  sys.argv[1:]:
  print "Processing file "+file
  if (os.path.exists(file+'.rmspickle2')):
    print "Found pickled data, skipping minimisation"
    f = open(file+'.rmspickle2', 'r')
    filedata.append(pickle.load(f))
    f.close()
  else:
    print "Could not find pickled minimisation data"
    tempdata = [ get_temperature(file), get_rmsd_data(get_structlist(file)) ]
    filedata.append(tempdata)
    f = open(file+'.rmspickle2', 'w')
    pickle.dump(tempdata,f)
    f.close()

filedata.sort()


f = open('rmsdhist.dat', 'w')
try:
  for i in range(len(filedata[0][1][2])):
    for j in range(i+1,len(filedata[0][1][2])):
      print >>f, filedata[0][1][2][i][j]
finally: 
  f.close()
  
#Outputs the min rmsd and the structure number
f = open('clusters.dat', 'w')
g = open('clusterfraction.dat', 'w')
h = open('avgrmsd.dat', 'w')
try:
  for data in filedata:
    clusters = cluster_analysis(data[1][2], 0.4, 10)
    
    sum = 0
    for cluster in clusters:
      sum += cluster[0]
      
    print >>f, data[0], len(clusters)
    print >>g, data[0], float(sum) / float(len(data[1][2]))
    print >>h, data[0], data[1][1], data[1][0]
finally:
  f.close()
  g.close()
  h.close()
  
