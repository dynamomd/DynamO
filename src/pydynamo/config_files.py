import math

import numpy

from pydynamo.file_types import XMLFile, validate_xmlfile
from pydynamo.math import conv_to_14sf


# A XMLFile/ElementTree but specialised for DynamO configuration files
class ConfigFile(XMLFile):
    def __init__(self, filename):
        super().__init__(filename)

    # Number of particles in the config
    def N(self):
        return int(len(self.tree.findall('.//Pt')))

    # Primary image volume
    def V(self):
        dims = self.image_dimensions()
        return dims[0] * dims[1] * dims[2]

    # number density
    def n(self):
        return self.N() / self.V()
    
    def __getitem__(self, key):
        return ConfigFile.config_props[key]['recalc'](self)
    
    def __str__(self):
        return "ConfigFile("+self._filename+")"

    def histogramTether(self, limits=[None, None, None]):
        tethers = self.tree.findall('.//Global/CellOrigins/Origin')
        particles = self.tree.findall('.//Pt/P')
        data = numpy.ndarray((len(tethers), 3))
        
        for idx, tether, particle in zip(range(len(tethers)), tethers, particles):
            p = numpy.array(list(map(float, [particle.attrib['x'], particle.attrib['y'], particle.attrib['z']])))
            t = numpy.array(list(map(float, [tether.attrib['x'], tether.attrib['y'], tether.attrib['z']])))
            data[idx] = p-t
        return numpy.histogramdd(data, range=[limits, limits, limits], bins=11)

    def histogramTether1D(self, limits=[None]):
        tethers = self.tree.findall('.//Global/CellOrigins/Origin')
        particles = self.tree.findall('.//Pt/P')
        data = numpy.ndarray((len(tethers), 1))
        for idx, tether, particle in zip(range(len(tethers)), tethers, particles):
            p = numpy.array(list(map(float, [particle.attrib['x'], particle.attrib['y'], particle.attrib['z']])))
            t = numpy.array(list(map(float, [tether.attrib['x'], tether.attrib['y'], tether.attrib['z']])))
            data[idx] = math.sqrt((p-t).dot(p-t))
        return numpy.histogramdd(data, bins=11)

    def image_dimensions(self):
        V = self.tree.find('.//SimulationSize')
        return [float(V.attrib['x']), float(V.attrib['y']), float(V.attrib['z'])]
    
    def to_freud(self):
        import freud
        N = self.N()
        frame = numpy.zeros((N, 3), dtype=numpy.float32)
        box = self.image_dimensions()
        box = freud.box.Box(Lx = box[0], Ly = box[1], Lz = box[2])
        particles = self.tree.findall('.//Pt/P')
        
        for idx, particle in enumerate(particles):
            frame[idx, 0] = numpy.float32(particle.attrib['x'])
            frame[idx, 1] = numpy.float32(particle.attrib['y'])
            frame[idx, 2] = numpy.float32(particle.attrib['z'])
        
        return box, frame
    
    config_props = {}

def validate_configfile(filename):
    return validate_xmlfile(filename)


# ###############################################
# #          Definition of state vars           #
# ###############################################
ConfigFile.config_props["N"] = {'recalculable':True, 'recalc': lambda config: config.N()}
ConfigFile.config_props["ndensity"] = {'recalculable':True, 'recalc': lambda config: conv_to_14sf(config.n())}
ConfigFile.config_props["InitState"] = {'recalculable':False, 'recalc': lambda config: "FCC"}

def Lambda_config(XMLconfig):
    tag = XMLconfig.tree.find('.//Interaction[@Type="SquareWell"]')
    if tag is None:
        return float('inf')
    else:
        return conv_to_14sf(float(tag.attrib['Lambda']))
ConfigFile.config_props["Lambda"] = {'recalculable':True, 'recalc': lambda config: Lambda_config}

def Rso_config(XMLconfig):
    tag = XMLconfig.tree.find('.//Global[@Type="SOCells"]')
    if tag is None:
        return float('inf')
    else:
        return conv_to_14sf(float(tag.attrib['Diameter']) / 2)    
ConfigFile.config_props["Rso"] = {'recalculable':True, 'recalc':Rso_config}


def PhiT_config(XMLconfig):
    Rso = Rso_config(XMLconfig)
    if Rso == float('inf'):
        return float('inf')
    density = XMLconfig.n()
    phiT = density * math.pi * (4.0/3.0) * (Rso**3)
    return conv_to_14sf(phiT)
def PhiT_gen(state):
    density = state["ndensity"]
    phiT = state["PhiT"]
    if phiT == float('inf'):
        Rso = float('inf')
        if 'Rso' in state:
            if state['Rso'] != float('inf'):
                raise RuntimeError("State contains conflicting Rso and PhiT")
        else:
            state['Rso'] = float('inf')
    else:
        Rso = (phiT / (density * math.pi *(4.0/3.0))) ** (1/3.0)
        if 'Rso' in state:
            if state['Rso'] != conv_to_14sf(Rso):
                raise RuntimeError("State contains conflicting Rso and PhiT")
        else:
            state['Rso'] = conv_to_14sf(Rso)
    return state
ConfigFile.config_props["PhiT"] = {'recalculable':True, 'recalc':PhiT_config, 'gen_state':PhiT_gen}


def kT_config(XMLconfig):
    tag = XMLconfig.tree.find('.//System[@Name="Thermostat"]')
    if tag is None:
        return float('inf')
    else:
        return conv_to_14sf(float(tag.attrib['Temperature']))
ConfigFile.config_props["kT"] = {'recalculable':True, 'recalc':kT_config}
