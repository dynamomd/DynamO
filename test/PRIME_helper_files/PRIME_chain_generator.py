#!/usr/bin/python2

from lxml import etree as ET
import sys
import os
import time
import random
import numpy as np
import subprocess
import re
from config import config as src_config

##############
#  Settings  #
##############

#Turn on for full output from LAMMPS and DynamO
#including LAMMPS .dcd files.
debug = False

###############
#  Functions  #
###############

def gauss():
    return str( random.gauss(0.0, 1.0) )

def mkdir_p(path):
    if not os.path.exists(path):
        os.makedirs(path)

def joinStr(list_like):
    return " ".join(str(item) for item in list_like)

def rotate_and_centre(coords_in):

    xrot = ( random.random()-0.5 )*np.pi
    yrot = ( random.random()-0.5 )*np.pi

    xrotarray = np.array([ [1.0, 0.0, 0.0], [0.0, np.cos(xrot), -np.sin(xrot)], [0.0, np.sin(xrot), np.cos(xrot)] ])
    yrotarray = np.array([ [np.cos(yrot), 0.0, np.sin(yrot)], [0.0, 1.0, 0.0], [-np.sin(yrot), 0.0, np.cos(yrot)] ])

    coords_out = np.dot(coords_in, xrotarray)
    coords_out = np.dot(coords_out, yrotarray)

    return coords_out - coords_out.mean(axis=0)

#############################################
###  Read in the command line parameters  ###
#############################################

nonglycine_SC    = list('ACDEFHIKLMNPQRSTVWY')
nonglycine_sites = list(nonglycine_SC) + ['NH', 'CH', 'CO']
sites            = list(nonglycine_sites) + ['G']
nonglycine_names = ['Alanine', 'Cysteine', 'Aspartic Acid', 'Glutamic Acid', 'Phenylalanine', 'Histidine', 'Isoleucine', 'Lysine', 'Leucine', 'Methionine', 'Asparagine', 'Proline', 'Glutamine', 'Arginine', 'Serine', 'Threonine', 'Valine', 'Tryptophan', 'Tyrosine' ] + ['Nitrogen+Hydrogen', 'Carbon+Hydrogen', 'Carbon+Oxygen']

try:
    filename = str(sys.argv[4])
except:
    filename = 'Hall_peptide'

psf_fn = filename + '.psf'
xml_fn = filename + '.xml'

try:
    temperature = sys.argv[3]
except:
    temperature = '1.0'

try:
    n_chains = int(sys.argv[2])
except:
    n_chains = 1

try:
    sys.argv[1] = sys.argv[1].upper()
    if ( sys.argv[1] == 'S1' ):
        sequence = list('PPPWLPYMPPWS')
    elif ( sys.argv[1] == 'N16N' ):
        sequence = list('AYHKKCGRYSYCWIPYDIERDRYDNGDKKC')
    elif ( sys.argv[1] == 'N16NN' ):
        sequence = list('AYHKKCGRYSYCWIPYNIQRNRYNNGNKKC')
    elif ( sys.argv[1][:5] == 'ALPHA' ):
        sequence = list('ACDEFGHIKLMNPQRSTVWY')
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

except (ValueError, IndexError) as e:
    print 'Run as ./peptide_maker.py (sequence) [n_chains = 1] [temperature kT = 1.0] [xml_fn = PRIME_peptide].'
    print ''
    raise

date               = time.strftime('%X %x %Z')
inner_box_size     = 40.0
box_pad            = 30.0
box_size           = 2*box_pad + inner_box_size
n_residues         = len(sequence)
n_bb_sites         = 3*n_residues
n_sc_sites         = n_residues - sequence.count('G')
n_sites            = n_bb_sites + n_sc_sites
atom_spacing       = 4.6

expanded_sequence = []
sidechain_IDs = []
for letter in sequence:
    expanded_sequence += ['NH','CH','CO']
    if letter != 'G':
        sidechain_IDs.append(len(expanded_sequence))
        expanded_sequence.append(letter)

nonglycine_expanded_sequence = filter(lambda a: a != 'G', expanded_sequence) #'Real' list, e.g. AGA gives NH,CH,CO,NH,CH,CO,NH,CH,CO,A,A

print 'Sequence:' , ''.join(sequence)
print 'File name:' , xml_fn
print 'N chains:', n_chains, '\n'

#######################################################################
#                     Populate atom coords array                      #
#######################################################################


source_res_count = len(src_config)/4

final_coords = np.empty( [n_chains*n_sites, 3] )
coords = np.empty( [n_sites, 3] )

for i_chain in range(n_chains):
    #start at a random value between 0 and source_res_count-1-len(sequence)
    startres = random.randint( 0, source_res_count-1-len(sequence) )
    prototype_coords = src_config[startres*4 : startres*4+len(sequence)*4]

    #Extract coords without coords for nonexistent glycine sidechain
    i_atom = 0
    for i_res, res in enumerate(sequence):

        if res != 'G':
            coords[i_atom:i_atom+4] = prototype_coords[i_res*4 : i_res*4+4]
            i_atom += 4
        else:
            coords[i_atom:i_atom+3] = prototype_coords[i_res*4 : i_res*4+3]
            i_atom += 3

    #Randomise position of whole chain
    min_distance = 0.0
    while min_distance < atom_spacing:
        #centre the atoms and add a random rotation
        rotated_and_centred = rotate_and_centre( coords )

        #add random displacement if multiple chains exist
        if n_chains > 1:
            final_chain_coords = rotated_and_centred + [ (random.random()-0.5)*inner_box_size, (random.random()-0.5)*inner_box_size, (random.random()-0.5)*inner_box_size, ]
        else:
            final_chain_coords = rotated_and_centred

        #check min_distance
        min_distance = 10.0
        for atom in final_chain_coords:
            for atom2 in final_coords[:i_chain*n_sites]:
                min_distance = min( min_distance, np.linalg.norm(atom - atom2) )

    final_coords[i_chain*n_sites:(i_chain+1)*n_sites] = final_chain_coords

#############################################
###               Set up XML              ###
#############################################

HB_strength = 1.0

DynamOconfig = ET.Element    ( 'DynamOconfig', attrib = {'version' : '1.5.0'} )

Simulation   = ET.SubElement ( DynamOconfig, 'Simulation', attrib = {'lastMFT':"-nan"} )
Properties   = ET.SubElement ( DynamOconfig, 'Properties' )
ParticleData = ET.SubElement ( DynamOconfig, 'ParticleData' )

####ParticleData section####
for ID in range( n_sites*n_chains ):
    ET.SubElement( ParticleData, 'Pt', attrib = {'ID' : str(ID) } )

    ET.SubElement( ParticleData[ID], 'P', attrib = {"x":str( final_coords[ID][0] ), "y":str( final_coords[ID][1] ), "z":str( final_coords[ID][2] )} )
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
for i_chain in range(n_chains):
    Molecule  = ET.SubElement ( Structure, 'Molecule', attrib = {'StartID':str(i_chain*n_sites), 'Sequence':''.join(sequence)} )

######################
#  Create PSF files  #
######################

psf_atoms_section = ""
psf_bonds_section = ""

for i_chain in range(n_chains):
    i_res=-1
    for i_atom, atom in enumerate(expanded_sequence):
        i_atom_global = i_atom+(i_chain*n_sites)

        if atom == 'NH':

            atom_label = 'N'

            i_res+=1
            if sequence[i_res-1] == 'G':
                bond_partner = i_atom_global - 1
            else:
                bond_partner = i_atom_global - 2

            #Set equal to itself to signal no bond partner:
            if i_res == 0:
                bond_partner = i_atom_global

        elif atom == 'CH':
            bond_partner = i_atom_global-1
            atom_label = 'CA'

        elif atom == 'CO':
            bond_partner = i_atom_global-1
            atom_label = 'C'

        else:
            bond_partner = i_atom_global-2
            atom_label = atom+'SC'

        psf_atoms_section += "{0: >8d} {1: <4} {2: <4d} {3: <4} {4: <4} {4: <4} {5: >10} {6: >13} {7: >11}\n".format(i_atom_global+1, str(i_chain), i_res, sequence[i_res], atom_label, "0.000000", "0.0000", "0")

        if bond_partner != i_atom_global:
            psf_bonds_section += "{0: >8d}{1: >8d}".format(bond_partner+1,i_atom_global+1)

            if len(psf_bonds_section) - psf_bonds_section.rfind("\n") > 63:
                psf_bonds_section += "\n"

with open(psf_fn, 'w') as psf_file:
    psf_file.write("PSF\n\n\t1 !NTITLE\n REMARKS " + ''.join(sequence) + " STRUCTURE FILE\n REMARKS DATE: " + date + "\n\n")
    psf_file.write("{0: >8d}".format(n_sites*n_chains) + " !NATOM\n" + psf_atoms_section + "\n")
    psf_file.write("{0: >8d}".format(n_sites-1) + " !NBOND\n" + psf_bonds_section + "\n\n")

#############################################
###              Write file               ###
#############################################

input_file = open(xml_fn, 'w')
input_file.write('<!-- DynamO input file contains ' + str(n_chains) + ' chain(s) of the sequence: ' + ''.join(sequence) + '. -->\n')
input_file.write('<!-- Created on ' +date + '. -->\n')
[ input_file.write(ET.tostring(DynamOconfig, pretty_print=True)) ]
input_file.close()

#Add thermostat and rescale via dynamod:
thermostat_command = [ 'dynamod',  '-T', temperature, '-r', temperature, '-o', xml_fn, '-Z', xml_fn ]
print "Running this command:", " ".join(thermostat_command)
if debug:
    print subprocess.Popen(thermostat_command, stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
else:
    silent_stdout = subprocess.Popen(thermostat_command, stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()

#Check config is valid with dynamod:
check_command = ['dynamod', xml_fn, "--check", '-o', xml_fn]
print "Running this command:", " ".join(check_command)
if debug:
    print subprocess.Popen(check_command, stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
else:
    silent_stdout = subprocess.Popen(check_command, stdout=subprocess.PIPE).communicate()

if debug:
    #Create trajectory file
    traj_command = ['dynamo2xyz', xml_fn]
    print "Running this command:", " ".join(traj_command)
    with open('traj.xyz', 'w') as trajfile:
        xyz = subprocess.check_output(traj_command)
        trajfile.write(xyz)
