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

import xml.etree.ElementTree as ET

basicDoc = ET.fromstring('<?xml version="1.0"?>'\
                         '<DynamOconfig version="1.5.0">'\
                         '<Simulation>'\
                         '<Scheduler Type="NeighbourList">'\
                         '<Sorter Type="BoundedPQMinMax3"/>'\
                         '</Scheduler>'\
                         '<SimulationSize x="10" y="10" z="10"/>'\
                         '<Genus>'\
                         '<Species Type="Point" Name="Bulk" Mass="1">'\
                         '<IDRange Type="All"/>'\
                         '</Species>'\
                         '</Genus>'\
                         '<BC Type="PBC"/>'\
                         '<Topology/>'\
                         '<Interactions/>'\
                         '<Locals/>'\
                         '<Globals/>'\
                         '<SystemEvents/>'\
                         '<Dynamics Type="Newtonian"/>'\
                         '</Simulation>'\
                         '<Properties/>'\
                         '<ParticleData/>'\
                         '</DynamOconfig>'\
)

def addParticle(doc, pos, vel):
    pdata = doc.find('./ParticleData')
    Pt = ET.SubElement(pdata, 'Pt')
    pos = list(map(str, pos))
    vel = list(map(str, vel))
    P = ET.SubElement(Pt, 'P', {'x':pos[0], 'y':pos[1], 'z':pos[2]})
    V = ET.SubElement(Pt, 'V', {'x':vel[0], 'y':vel[1], 'z':vel[2]})


def prettyprint(doc):
    import xml.dom.minidom
    xml = xml.dom.minidom.parseString(ET.tostring(doc))
    return xml.toprettyxml(indent=' ')

def isclose(a,b,tol):
    return (abs(a-b) <= abs(tol*a)) or (abs(a-b) <= abs(tol*b))
