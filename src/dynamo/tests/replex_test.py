#!/usr/bin/env python3
#   dynamo:- Event driven molecular dynamics simulator 
#   http://www.dynamomd.org
#   Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>
#
#   This program is free software: you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   version 3 as published by the Free Software Foundation.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
import dynamo
import os
import math
import sys
import getopt
import xml.etree.ElementTree as ET

e_swaps=2500
p_swaps=25000
swap_time=0.5
run=True

#This may not exactly correspond to the real finish time, it may be +-
#swaptime
e_finish_time=float(e_swaps)*swap_time 
p_finish_time=float(p_swaps)*swap_time 


error_count = 0

shortargs=""
longargs=["dynarun=", "dynamod=", "dynahist_rw="]
try:
    options, args = getopt.gnu_getopt(sys.argv[1:], shortargs, longargs)
except getopt.GetoptError as err:
    # print help information and exit:
    print(str(err)) # will print something like "option -a not recognized"
    sys.exit(2)

dynarun_cmd="NOT SET"
dynamod_cmd="NOT SET"
dynahist_rw_cmd="NOT SET"

for o,a in options:
    if o == "--dynarun":
        dynarun_cmd = a
    if o == "--dynamod":
        dynamod_cmd = a
    if o == "--dynahist_rw":
        dynahist_rw_cmd = a


for name,exe in [("dynahist_rw", dynahist_rw_cmd), ("dynamod", dynamod_cmd), ("dynarun", dynarun_cmd)]:
    if not(os.path.isfile(exe) and os.access(exe, os.X_OK)):
        raise RuntimeError("Failed to find "+name+" executabe at "+exe)
        
import subprocess

###### INITIALISATION
Temperatures=[1.0, 2.0, 5.0, 10.0]
for i,T in enumerate(Temperatures):
    cmd=[dynamod_cmd, "-m2", "-T"+str(T), "-oc"+str(i)+".xml"]
    print(" ".join(cmd))
    if run:
        subprocess.call(cmd)

###### EQUILIBRATION
cmd=[dynarun_cmd, "--engine=2", "-oc%ID.xml", "--out-data-file=o%ID.xml", "-N4", "-i"+str(swap_time), "-f"+str(e_finish_time)]+["c"+str(i)+".xml" for i in range(len(Temperatures))]
print(" ".join(cmd))
if run:
    subprocess.call(cmd)


###### PRODUCTION
cmd=[dynarun_cmd, "--engine=2", "-oc%ID.xml", "--out-data-file=o%ID.xml", "-N4", "-i"+str(swap_time), "-LIntEnergyHist", "-f"+str(p_finish_time)]+["c"+str(i)+".xml" for i in range(len(Temperatures))]
print(" ".join(cmd))
if run:
    subprocess.call(cmd)

###### OUTPUT VALIDATION

import re
replex_calls=None
with open("replex.stats") as f:
    for line in f:
        match = re.search("Number_of_replex_cycles ([0-9]+)", line)
        if match:
            replex_calls = int(match.group(1))
            break

if replex_calls is None:
    raise RuntimeError("Could not parse the number of replex calls performed")

if abs(replex_calls - p_swaps) > 1:
    raise RuntimeError("The number of actual swaps ("+str(replex_calls)+") is different to the number of requested swaps +("+str(p_swaps)+")")

expectedCvs={10.0: 0.304967, 5.0: 2.52568, 2.0: 20.5657, 1.0: 23.1892}

for i,T in enumerate(Temperatures):
    outputfile="o"+str(i)+".xml"
    xmldoc=ET.parse(outputfile)

    measured_T=float(xmldoc.getroot().find(".//Temperature").attrib["Mean"])
    if not dynamo.isclose(measured_T, T, 4e-2):
        error_count = error_count + 1
        print("Simulation Temperature is different to what is expected:"+str(measured_T)+"!="+str(T))
    
    simtime=float(xmldoc.getroot().find(".//Duration").attrib["Time"])
    expected_simtime = swap_time * replex_calls * math.sqrt(min(Temperatures)/T)
    if not dynamo.isclose(simtime, expected_simtime, 1e-4):
        error_count = error_count + 1
        print("Simulation duration is different to what is expected:"+str(simtime)+"!="+str(expected_simtime))

    measured_Cv=float(xmldoc.getroot().find(".//ResidualHeatCapacity").attrib["Value"])
    expected_Cv=expectedCvs[T]

    if not dynamo.isclose(measured_Cv, expected_Cv, 0.1):
        error_count = error_count + 1
        print("Simulation heat capacity is different to what is expected:"+str(measured_Cv)+"!="+str(expected_Cv))

###### dynahist_rw VALIDATION
cmd=[dynahist_rw_cmd]+["o"+str(i)+".xml" for i in range(len(Temperatures))]
print(" ".join(cmd))
if run:
    subprocess.call(cmd)

for line,ref in zip(open('logZ.out', 'r'), [0, -37.1325363028306097, -50.8797611158942076, -53.720182897102171]):
    logZval = float(line.split()[1])
    if ref == 0:
        if logZval != 0:
            raise RuntimeError("First logZ value is not zero")
    else:
        if not dynamo.isclose(logZval, ref, 1e-2):
            error_count = error_count + 1
            print("Calculated logZ value is incorrect:"+str(logZval)+"!="+str(ref))


def interpolate(yin, xin, xout):
  i1 = None
  for i in range(len(xin)):
      if xin[i] == xout:
          return yin[i]
      if xin[i] < xout:
          i1 = i
          break
  if i1 is None:
      raise ValueError("xout out of range")

  x0 = xin[i1-1]  
  x1 = xin[i1]  
  y0 = yin[i1-1]  
  y1 = yin[i1] 
  
  return (xout - x0) / (x1 - x0) * (y1 - y0) + y0  

Cvdata=[list(map(float, line.split())) for line in open("Cv.out")]

for T in expectedCvs:
    measuredCv = interpolate([data[1] for data in Cvdata],  [data[0] for data in Cvdata], T)
    expectedCv = expectedCvs[T]
    if not dynamo.isclose(measuredCv, expectedCv, 1e-1):
        error_count = error_count + 1
        print("Histogram reweighted Cv value is incorrect:"+str(measuredCv)+"!="+str(expectedCv))

print("Total errors:", error_count)
print("Only reporting a fatal error if there are multiple errors (one failures is \"normal\")")
sys.exit(error_count > 1)
