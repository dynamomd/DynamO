#!/usr/bin/python2.7

from lxml import etree as ET
import sys
import os
import time
import random
import numpy as np
import subprocess
import prototype_positions as pp
import re

#Turn on for full output from DynamO, and a .xyz file output
debug = False

def gauss():
    return str( random.gauss(0.0, 1.0) )

#############################################
###  Read in the command line parameters  ###
#############################################

nonglycine_SC    = list('ACDEFHIKLMNPQRSTVWY')
nonglycine_sites = list(nonglycine_SC) + ['NH', 'CH', 'CO']
sites            = list(nonglycine_sites) + ['G']

HB_strength = 1.0

try:
    filename = str(sys.argv[3])
except:
    filename = 'PRIME_peptide'

psf_fn = filename + '.psf'
xml_fn = filename + '.xml'

try:
    temperature = sys.argv[2]
except:
    temperature = '1.0'

try:
    sys.argv[1] = sys.argv[1].upper()
    if ( sys.argv[1] == 'S1' ):
        sequence = list('PPPWLPYMPPWS')
    elif ( sys.argv[1] == 'N16N' ):
        sequence = list('AYHKKCGRYSYCWIPYDIERDRYDNGDKKC')
    elif ( sys.argv[1][:5] == 'ALPHA' ):
        sequence = list('ACDEFHIKLMNPQRSTVWY')
    else:
        sequence = list(sys.argv[1])
        valid_input = sites + [ str(number) for number in range(0,10) ]
        assert ( [ residue in valid_input for residue in list(sys.argv[1])] )

        matches = re.finditer("[0-9]+", sys.argv[1])

        for match in matches:

            letter_index = match.start()-1
            number = int( match.group() )

            letter = sys.argv[1][letter_index]
            insertion = letter*number
            sequence[letter_index] = insertion

        sequence = filter( lambda x: x.isalpha(), sequence)
        sequence = list(''.join(sequence))

        if len(sequence) > 33:
            sys.exit('33 residues is currently the maximum length for this generator.')

except (ValueError, IndexError):
    raise
    sys.exit('Run as ./peptide_maker.py (sequence) [temperature kT = 1.0] [xml_fn = PRIME_peptide].')

date               = time.strftime('%X %x %Z')
box_size           = 200
n_residues         = len(sequence)
n_bb_sites         = 3*n_residues
n_sc_sites         = n_residues - sequence.count('G')
n_sites            = n_bb_sites + n_sc_sites

expanded_sequence = []
for letter in sequence:
    expanded_sequence += ['NH','CH','CO']
    if letter != 'G':
        expanded_sequence.append(letter)

print 'Sequence:' , ''.join(sequence)
print 'File name:' , xml_fn , '\n'

##########################
###      Geometry      ###
##########################

coords = np.zeros([ len(expanded_sequence), 3 ], dtype=float)
j=0
for i_res, res in enumerate(sequence):

    for i_atom in range(3):
        bb_only_index = i_res*3 + i_atom
        coords[j] = pp.BB_33[bb_only_index]
        j += 1

    if res != 'G':
        res_number = nonglycine_sites.index(res)

        coords[j] = pp.SC_33[ res_number ][ i_res ]
        j += 1

#############################################
###               Set up XML              ###
#############################################

DynamOconfig = ET.Element    ( 'DynamOconfig', attrib = {'version' : '1.5.0'} )

Simulation   = ET.SubElement ( DynamOconfig, 'Simulation', attrib = {'lastMFT':"-nan"} )
Properties   = ET.SubElement ( DynamOconfig, 'Properties' )
ParticleData = ET.SubElement ( DynamOconfig, 'ParticleData' )

####ParticleData section####
for ID in range( n_sites ):
    ET.SubElement( ParticleData, 'Pt', attrib = {'ID' : str(ID) } )

    ET.SubElement( ParticleData[ID], 'P', attrib = {"x":str( coords[ID][0] ), "y":str( coords[ID][1] ), "z":str( coords[ID][2] )} )
    ET.SubElement( ParticleData[ID], 'V', attrib = {"x":gauss(), "y":gauss(), "z":gauss()} )

####Simulation section####

Scheduler      = ET.SubElement ( Simulation, 'Scheduler',      attrib = {'Type':'NeighbourList'} )
SimulationSize = ET.SubElement ( Simulation, 'SimulationSize', attrib = dict(zip(['x','y','z'], [str(box_size)]*3)) )
BC             = ET.SubElement ( Simulation, 'BC',             attrib = {'Type':'PBC'} )
Genus          = ET.SubElement ( Simulation, 'Genus')
Topology       = ET.SubElement ( Simulation, 'Topology' )
Interactions   = ET.SubElement ( Simulation, 'Interactions')
Locals         = ET.SubElement ( Simulation, 'Locals' )
Globals        = ET.SubElement ( Simulation, 'Globals' )
SystemEvents   = ET.SubElement ( Simulation, 'SystemEvents' )
Dynamics       = ET.SubElement ( Simulation, 'Dynamics', attrib = {'Type':'Newtonian'} )

Sorter = ET.SubElement ( Scheduler, 'Sorter', attrib = {'Type':'BoundedPQMinMax3'} )

#Add each species and its list of IDs to the XML tree
PRIME_species = ET.SubElement( Genus, 'Species', attrib = {'Mass':'PRIMEData', 'Name':'PRIMEGroups', 'Type':'Point'} )
ET.SubElement( PRIME_species, 'IDRange', attrib = {'Type':'All'} )

temp = ET.SubElement( Globals, 'Global', attrib = {'Type':'Cells','Name':'SchedulerNBList','NeighbourhoodRange':'7.400000000000e+00'})
ET.SubElement( temp, 'IDRange', attrib = {'Type':'All'})

#Interactions section
PRIME = ET.SubElement( Interactions, 'Interaction', attrib = {'Type':'PRIME', 'Name':'Backbone', 'Topology':"PRIMEData", 'HBStrength':str(HB_strength)} )
ET.SubElement( PRIME, 'IDPairRange', attrib = {'Type':'All'} )

#Topology section
Structure = ET.SubElement ( Topology, 'Structure', attrib = {'Type':'PRIME', 'Name':'PRIMEData'} )
Molecule  = ET.SubElement ( Structure, 'Molecule', attrib = {'StartID':'0', 'Sequence':''.join(sequence)} )

######################
#  Create PSF files  #
######################

psf_atoms_section = ""
psf_bonds_section = ""

i_res=-1
for i_atom, atom in enumerate(expanded_sequence):

    if atom == 'NH':
        i_res+=1
        if sequence[i_res-1] == 'G':
            bond_partner = i_atom - 1
        else:
            bond_partner = i_atom - 2

        #Set equal to itself to signal no bond partner:
        if i_res == 0:
            bond_partner = i_atom

    elif atom == 'CH' or atom == 'CO':
        bond_partner = i_atom-1

    else:
        bond_partner = i_atom-2

    psf_atoms_section += "{0: >8d} {1: <4} {2: <4d} {3: <4} {4: <4} {4: <4} {5: >10} {6: >13} {7: >11}\n".format(i_atom+1, str(0), i_res, sequence[i_res], atom, "0.000000", "0.0000", "0")

    if bond_partner != i_atom:
        psf_bonds_section += "{0: >8d}{1: >8d}".format(bond_partner+1,i_atom+1)

        if len(psf_bonds_section) - psf_bonds_section.rfind("\n") > 63:
            psf_bonds_section += "\n"

with open(psf_fn, 'w') as psf_file:
    psf_file.write("PSF\n\n\t1 !NTITLE\n REMARKS " + ''.join(sequence) + " STRUCTURE FILE\n REMARKS DATE: " + date + "\n\n")
    psf_file.write("{0: >8d}".format(n_sites) + " !NATOM\n" + psf_atoms_section + "\n")
    psf_file.write("{0: >8d}".format(n_sites-1) + " !NBOND\n" + psf_bonds_section + "\n\n")

####################
#  Write XML file  #
####################

input_file = open(xml_fn, 'w')
input_file.write('<!-- DynamO input file contains the PRIME20 model of the sequence: ' + ''.join(sequence) + '. -->\n')
input_file.write('<!-- Created on ' +date + '. -->\n')
[ input_file.write(ET.tostring(DynamOconfig, pretty_print=True)) ]
input_file.close()

#Add thermostat and rescale via dynamod:
thermostat_command = [ 'dynamod',  '-T', temperature, '-r', temperature, '-o', xml_fn, '-Z', xml_fn ]
print "Running this command:", " ".join(thermostat_command)
if debug:
    print subprocess.check_output(thermostat_command)
else:
    silent_stdout = subprocess.check_output(thermostat_command)

#Check config is valid with dynamod:
run_command = ['dynamod', xml_fn, "--check"]
print "Running this command:", " ".join(run_command)
if debug:
    print subprocess.check_output(run_command)
else:
    silent_stdout = subprocess.check_output(run_command)

if debug:
    #Create trajectory file
    traj_command = ['dynamo2xyz', xml_fn]
    print "Running this command:", " ".join(traj_command)
    with open('traj.xyz', 'w') as trajfile:
        xyz = subprocess.check_output(traj_command)
        trajfile.write(xyz)
