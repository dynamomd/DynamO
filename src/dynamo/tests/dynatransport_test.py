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
import os
import math
import sys
import getopt
import xml.etree.ElementTree as ET

run=True

shortargs=""
longargs=["dynarun=", "dynamod=", "dynahist_rw=", "dynatransport=", "python="]
try:
    options, args = getopt.gnu_getopt(sys.argv[1:], shortargs, longargs)
except getopt.GetoptError as err:
    # print help information and exit:
    print(str(err)) # will print something like "option -a not recognized"
    sys.exit(2)

dynarun_cmd="NOT SET"
dynamod_cmd="NOT SET"
dynahist_rw_cmd="NOT SET"
dynatransport_cmd="NOT SET"
python_cmd="NOT SET"

for o,a in options:
    if o == "--dynarun":
        dynarun_cmd = a
    if o == "--dynamod":
        dynamod_cmd = a
    if o == "--dynahist_rw":
        dynahist_rw_cmd = a
    if o == "--dynatransport":
        dynatransport_cmd = a
    if o == "--python":
        python_cmd = a

for name,exe in [("dynahist_rw", dynahist_rw_cmd), ("dynamod", dynamod_cmd), ("dynarun", dynarun_cmd), ("dynatransport", dynatransport_cmd), ("python", python_cmd)]:
    if not(os.path.isfile(exe) and os.access(exe, os.X_OK)):
        raise RuntimeError("Failed to find "+name+" executabe at "+exe)
        
import subprocess
cmd=[dynamod_cmd, "-m0", "-d0.5", "-oc.xml"]
print(" ".join(cmd))
if run:
    subprocess.call(cmd)

###### EQUILIBRATION
cmd=[dynarun_cmd, "c.xml", "-oc.xml", "--out-data-file=o.xml", "-c300000"]
if run:
    subprocess.call(cmd)

###### PRODUCTION
cmd=[dynarun_cmd, "c.xml", "-oc.xml", "--out-data-file=o.xml", "-c2000000"]
if run:
    subprocess.call(cmd)

cmd=[python_cmd, dynatransport_cmd, "o.xml", "-c1.0", "-s0.3"]
out = subprocess.check_output(cmd)
visc=float(out.decode().split("\n")[0].split()[1])
thermal=float(out.decode().split("\n")[4].split()[1])

def isclose(a,b,tol):
    return (abs(a-b) <= abs(tol*a)) or (abs(a-b) <= abs(tol*b))

if not isclose(visc, 0.58, 2*0.06/0.58):
    raise RuntimeError("Simulation viscosity is different to what is expected:"+str(visc)+"!="+str(0.58))

if not isclose(thermal, 2.3, 2 * 0.2/2.3):
    raise RuntimeError("Simulation thermal conductivity is different to what is expected:"+str(thermal)+"!="+str(2.3))

