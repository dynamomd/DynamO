#!/usr/bin/env python
import os
import xml.etree.ElementTree as ET
import sys

#A helpful function to load compressed or uncompressed XML files
def loadXMLFile(filename):
    #Check if the file is compressed or not, and 
    if (os.path.splitext(filename)[1][1:].strip() == "bz2"):
        import bz2
        f = bz2.BZ2File(filename)
        doc = ET.parse(f)
        f.close()
        return doc
    else:
        return ET.parse(filename)

if __name__ == "__main__":
    """A program to take in a bulk fluid configuration file and add a hard wall to the center
    while removing all fluid particles which overlap the wall."""

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', required=True, type=str, help="configuration file")
    parser.add_argument('-r', '--radius', type=float, help="radius of the wall to be inserted")
    parser.add_argument('-N', '--n_particles', type=int, help="number of particles to keep in system")
    parser.add_argument('-L', '--Length', type=float, help="side length of simulation box")
    parser.add_argument('-o', '--output', default="config.fixed.xml", type=str, 
                       help="name of the output config file to use")
    args = parser.parse_args()
    
    radius = args.radius
    filename = args.file
    particles = args.n_particles
    length = args.Length
    ofile = args.output

    XMLDoc=loadXMLFile(filename)
    RootElement=XMLDoc.getroot()
    deleteme = []
    
    oldlength = RootElement.find('.//SimulationSize').get('x')
    if float(oldlength) > length:
        sys.exit('original configuration size larger than length' \
                 + ' specified with -L flag. Either increase -L argument '\
                 + 'or use configuration with smaller size')

    #Calculate distance and delete overlapping particles 
    for particle in RootElement.iter('Pt'):
        xyz = particle.find('P')
        x = float(xyz.attrib['x'])
        y = float(xyz.attrib['y'])
        z = float(xyz.attrib['z'])

        r2 = x**2 + y**2 + z**2
        if (r2 <= (radius+0.5)**2):
            deleteme.append(particle)

    for particle in deleteme:
        RootElement.find("ParticleData").remove(particle)

    #Reindex fluid particles from 1 to N. Wall has index 0
    index = 0
    for particle in RootElement.iter('Pt'):
        index += 1
        particle.set("ID",str(index))

    if particles and particles > index:
            sys.exit('The number of particles requested is larger than'\
                 + ' the number of particles remaining in the system.'\
                 + ' Use an initial configuration with a higher density'\
                 + ' or request a smaller number of particles.')
    
    #Fix fluid species def and add species def for wall
    genus = RootElement.find(".//Genus")
    bulkspecies = genus.find('./Species[@Name="Bulk"]')
    bulkspecies.find('./IDRange').set("Type","Ranged")
    bulkspecies.find('./IDRange').set("Start","1")
    if particles:
        bulkspecies.find('./IDRange').set("End",str(particles))
    else:
        bulkspecies.find('./IDRange').set("End",str(index))
    ET.SubElement(genus, "Species", {"Mass":"100", "Name":"Wall", "Type":"Point"})
    wallspecies = genus.find('*[@Name="Wall"]')
    ET.SubElement(wallspecies, "IDRange", {"Type":"List"})
    ET.SubElement(wallspecies.find("IDRange"),"ID",{"val":"0"})

    #Add wall particle to ParticleData
    Pdata = RootElement.find('ParticleData')
    ET.SubElement(Pdata,"Pt",{"ID":"0"})
    ET.SubElement(Pdata.find('./Pt[@ID="0"]'),"P",{"x":"0","y":"0","z":"0"})
    ET.SubElement(Pdata.find('./Pt[@ID="0"]'),"V",{"x":"0","y":"0","z":"0"})
    wallelement = Pdata.find('./Pt[@ID="0"]')
    Pdata.remove(Pdata.find('./Pt[@ID="0"]'))
    Pdata.insert(0,wallelement) #################################
    #Add wall/fluid interaction and wall self interaction, fix range on fluid/fluid
    # interaction
    interactions = RootElement.find(".//Interactions")
    bulkinteraction = interactions.find('./Interaction[@Name="Bulk"]')
    bulkinteraction.find('IDPairRange').set("Type","Single")
    if particles:
        ET.SubElement(bulkinteraction.find('IDPairRange'),"IDRange",{"Type":"Ranged","Start":"1",
                     "End":str(particles)})
    else:
        ET.SubElement(bulkinteraction.find('IDPairRange'),"IDRange",{"Type":"Ranged","Start":"1",
    "End":str(index)})
    
    ET.SubElement(interactions,"Interaction",{"Type":"HardSphere","Diameter":str(radius+0.5),
    "Name":"fluidwall"})
    fluidwall = interactions.find('./Interaction[@Name="fluidwall"]')
    ET.SubElement(fluidwall,"IDPairRange",{"Type":"Pair"})
    if particles:
        ET.SubElement(fluidwall.find('IDPairRange'),"IDRange",{"Type":"Ranged","Start":"1",
                     "End":str(particles)})
    else:
        ET.SubElement(fluidwall.find('IDPairRange'),"IDRange",{"Type":"Ranged","Start":"1",
                     "End":str(index)})
    ET.SubElement(fluidwall.find('IDPairRange'),"IDRange",{"Type":"List"})
    ET.SubElement(fluidwall.find('IDPairRange/IDRange[@Type="List"]'),"ID",{"val":"0"})
    
    ET.SubElement(interactions,"Interaction",{"Type":"HardSphere","Diameter":"1",
    "Name":"wall"})
    wallself = interactions.find('./Interaction[@Name="wall"]')
    ET.SubElement(wallself,"IDPairRange",{"Type":"Self"})
    ET.SubElement(wallself.find('IDPairRange'),"IDRange",{"Type":"List"})
    ET.SubElement(wallself.find('IDPairRange/IDRange'),"ID",{"val":"0"})
    
    #Remove excess particles if -n flag is used
    if particles:
        for id in range(particles+1,index+1):
            deletep = RootElement.find('.//Pt[@ID="'+str(id)+'"]')
            RootElement.find("ParticleData").remove(deletep)
    
    #change x and y dimensions if requested with -L flag
    if length:
       RootElement.find('.//SimulationSize').set('x',str(length))
       RootElement.find('.//SimulationSize').set('y',str(length))
    RootElement.find('.//SimulationSize').set('z',str(4.0*(radius+0.5)))
    #ET.dump(RootElement.find("./Simulation"))
    
    XMLDoc.write(ofile)
