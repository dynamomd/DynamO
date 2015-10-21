#!/usr/bin/env python
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

run=True

shortargs=""
longargs=["dynarun=", "dynamod=", "dynahist_rw="]
try:
    options, args = getopt.gnu_getopt(sys.argv[1:], shortargs, longargs)
except getopt.GetoptError as err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
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
        raise RuntimeError("Failed to find "+name+" executabe at "+exe+"\,"+str(sys.argv))
        
doc = dynamo.basicDoc

dynamo.addParticle(doc, pos=[0,0,0], vel=[1,0,0])
dynamo.addParticle(doc, pos=[0.99,0,0], vel=[-1,0,0])


interactions = doc.find('./Simulation/Interactions')
interaction = ET.SubElement(interactions, 'Interaction', {'Type':'Stepped', 'LengthScale':'1', 'EnergyScale':'1', 'Name':'Bulk'})
ET.SubElement(interaction, 'IDPairRange', {'Type':'All'})
potential = ET.SubElement(interaction, 'Potential', {'Type':'Stepped', 'Direction':'Left'})
ET.SubElement(potential, 'Step', {'R':'1.0', 'E':'-0.5'})
ET.SubElement(potential, 'Step', {'R':'0.5', 'E':'-0.75'})

systemevents = doc.find('./Simulation/SystemEvents')
system = ET.SubElement(systemevents, 'System', {'Type':'Andersen', 'MFT':'0.6', 'Temperature':'1', 'Name':'Thermostat', 'SetPoint':'0.05', 'SetFrequency':'100'})
ET.SubElement(system, 'IDRange', {'Type':'All'})

dynamics = doc.find('./Simulation/Dynamics')
dynamics.attrib['Type'] = 'NewtonianMCCMap'
dynamics.attrib['Interaction'] = 'Bulk'
potential = ET.SubElement(dynamics, 'Potential')
maptag = ET.SubElement(potential, 'Map', {'W':'-1000000', 'Distance':'0'})
ET.SubElement(maptag, 'Contact', {'ID1':'0', 'ID2':'1', 'State':'1'})
maptag = ET.SubElement(potential, 'Map', {'W':'-1000000', 'Distance':'0'})
ET.SubElement(maptag, 'Contact', {'ID1':'0', 'ID2':'1', 'State':'2'})

open('test.xml', 'w').write(dynamo.prettyprint(doc))

import subprocess

cmd=[dynarun_cmd, "test.xml", '-c1000000' ,'--out-data-file=output.xml']
print " ".join(cmd)
subprocess.call(cmd)

output = ET.parse('output.xml')
datatag=output.find('./Misc/UConfigurational')
datatag.attrib['Mean']
datatag.attrib['Min']
datatag.attrib['Max']
